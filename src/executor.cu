#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <unistd.h>
#include <chrono>

#include <nvToolsExtCudaRt.h>

#include <type_traits>

#include "swift.h"
#include "util.cuh"
#include "util.h"
#include "types.h"
#include "executor.cuh"
#include "wh.cuh"
#include "convert.h"

#define VEC_EL_SIZE(x) sizeof(std::remove_reference<decltype(x)>::type::value_type)

namespace sr
{
	namespace exec
	{
		using hrclock = std::chrono::high_resolution_clock;

		double dt_ms(hrclock::time_point t1, hrclock::time_point t2)
		{
			std::chrono::duration<double, std::milli> millis = t2 - t1;
			return millis.count();
		}

		using namespace sr::wh;
		using namespace sr::util;
		using namespace sr::convert;
		using namespace sr::data;

		struct DeviceParticleUnflaggedPredicate
		{
			template <typename Tuple>
			__host__ __device__ bool operator()(const Tuple &args)
			{
				uint16_t flag = thrust::get<2>(thrust::get<0>(args));
				return flag == 0;
			}
		};

		struct DeviceParticleAlivePredicate
		{
			template <typename Tuple>
			__host__ __device__ bool operator()(const Tuple &args)
			{
				// 0x00FE = 0000 0000 1111 1110
				uint16_t flag = thrust::get<2>(thrust::get<0>(args));
				return ((flag & 0xFE) == 0) && (flag != 0x01);
			}
		};

		struct KillEncounterKernel
		{
			KillEncounterKernel() {}

			template <typename Tuple>
			__host__ __device__ void operator()(Tuple args) const
			{
				uint16_t flags = thrust::get<2>(args);
				// If a particle is inside any big bodies:
				if ((flags & 0x01) == 0x01 && (flags != 0x01))
				{
					// 0xFF00 = 1111 1111 0000 0000
					flags = static_cast<uint16_t>((flags & 0xFF00) | 0x80);
				}
				thrust::get<2>(args) = flags;
			}
		};

		Executor::Executor(HostData &_hd, DeviceData &_dd, const Configuration &_config, std::ostream &out)
			: hd(_hd), dd(_dd), output(out), config(_config), resync_counter(0) {}

		// this should be called after hd is populated
		void Executor::init()
		{
			alloc_packed();

			// temp_log.open("temp_log_executor.txt");

			// glisse only supports helio
			// sr::convert::to_helio(hd);

			if (config.interp_planets)
			{
				// setup interpolator
				interpolator = sr::interp::Interpolator(config, hd.planets, config.planet_history_file, config.read_binary_hist, config.read_single_hist);

				// step interp forward until we're in the correct interval
				while (interpolator.t1 <= t)
				{
					interpolator.next(hd.planets);
				}

				// make sure we call next at least once, otherwise it means that the lookup file started past
				// the current time

				// interpolator.t0 is initialized to +inf, but calling next once makes it finite
				// so if t0 is not infinity, next was called at least once
				ASSERT(!std::isinf(interpolator.t0), "lookup file starts past the given Initial-Time")

				// we might have started on a point that's in the middle of a lookup interval:
				// the stuff below handles that case

				// cur_ts represents the current timestep number in the current interval
				// even if we start in the middle of the interval, set cur_ts to 0
				interpolator.cur_ts = 0;

				// n_ts is the total number of timesteps in the current interval
				// the number of timesteps is the remaining time in the interval divided by dt
				interpolator.n_ts = std::max<size_t>(1, static_cast<size_t>(std::round((interpolator.t1 - t) / config.dt)));

				// set the effective dt appropriately
				interpolator.eff_dt = (interpolator.t1 - t) / static_cast<double>(interpolator.n_ts);

				// relative time since the beginning of the interval
				interpolator.rel_t = t - interpolator.t0;

				// fill the planet location at t=0
				interpolator.fill_one(hd.planets, interpolator.rel_t);
			}

			// setup integrator - this gives all particles an acceleration and also detects initial encounters
			// needs to happen after interpolator init to get initial planet positions
			integrator = sr::wh::WHCudaIntegrator(hd.planets, hd.particles, config, htd_stream);

			// setup encounter integrator
			swift = sr::swift::SwiftEncounterIntegrator(config, hd.particles.n());

			// calculate initial energy
			calculate_planet_metrics(hd.planets, &e_0, nullptr);

			output << std::setprecision(7);
			output << "e_0 (planets) = " << e_0 << std::endl;
			output << "n_particle = " << hd.particles.n() << std::endl;
			output << "n_particle_alive = " << hd.particles.n_alive() << std::endl;
			output << "==================================" << std::endl;
			output << "Sending initial conditions to GPU." << std::endl;

			// create cuda streams
			cudaStreamCreate(&main_stream);
			cudaStreamCreate(&sort_stream);
			cudaStreamCreate(&htd_stream);
			cudaStreamCreate(&dth_stream);

			nvtxNameCudaStream(main_stream, "main");
			nvtxNameCudaStream(sort_stream, "sort");
			nvtxNameCudaStream(htd_stream, "htd");
			nvtxNameCudaStream(dth_stream, "dth");

			cudaEventCreate(&start_event);
			cudaEventCreate(&gpu_finish_event);

			// more initialization
			rollback_state = DeviceParticlePhaseSpace(hd.particles.n());
			dd.particles = DeviceParticlePhaseSpace(hd.particles.n());

			// having two planet data arrays allows us to update planet data while the other one is
			// being used by the integrator, so there is always one being used to integrate while
			// the other one is used to upload data for the next timechunk
			dd.planets0 = DevicePlanetPhaseSpace(hd.planets.n(), config.tbsize);
			dd.planets1 = DevicePlanetPhaseSpace(hd.planets.n(), config.tbsize);

			// planet_data_id chooses whether to use planets0 or planets1, arbitrarily choose 0 to start with
			dd.planet_data_id = 0;

			// ** INITIALIZE PARTICLES

			// partition unflagged (any particles with zero low byte, or deathflag & 0x00FF == 0)
			hd.particles.stable_partition_unflagged(0, hd.particles.n());

			// now, n_alive = n_unflagged
			size_t n_unflagged = hd.particles.n_alive();

			// particles having close encounters with any planets (sun excluded)
			hd.particles.stable_partition_alive(n_unflagged, hd.particles.n() - n_unflagged);

			// now, n_alive = total alive (enc + unflagged)
			size_t n_alive = hd.particles.n_alive();

			// don't need to upload if no particles
			if (hd.particles.n() > 0)
			{
				upload_data(0, hd.particles.n());
				dd.particle_phase_space().n_alive = n_alive;
			}

			// download data right after uploading -
			// i forgot why we need this, but I think it's just so the GPU data can be debugged
			download_data(0, hd.particles.n());

			starttime = hrclock::now();

			output << "       Starting simulation.       " << std::endl
				   << std::endl;

			discard_output = std::ofstream(sr::util::joinpath(config.outfolder, "discard.out"));
			if (config.write_encounter_log)
				*encounter_output << std::setprecision(17);
			discard_output << std::setprecision(17);

			hd.planets_snapshot = hd.planets.base;

			// upload planet data before the first timechunk
			update_planets();

			// write all dead particles before the integration starts.
			write_encounter(n_alive, hd.particles.n(), t);

			// out_timing << "t,      n_enc,  enc,    wswift, iswift, swift,  rswift, dswift, io,     planet, planetu,encup,  sort,   dl,     endenc, resync, cpu,    gpu,    total" << std::endl;
		}

		void Executor::swap_logs()
		{
			// swap new and old logs
			hd.planets.swap_logs();
			integrator.swap_logs();
		}

		void Executor::update_planets()
		{
			prev_dt = cur_dt;
			prev_tbsize = cur_tbsize;

			hrclock::time_point clock = hrclock::now();

			// std::cout << t << std::endl;
			integrator.recalculate_rh(hd.planets);

			if (config.interp_planets)
			{
				// cur_ts = n_ts means that we've reached the end of the lookup interval
				// need to step to the next one
				if (interpolator.cur_ts == interpolator.n_ts)
				{
					// sanity check - make sure that rel_t = t1 - t0
					double diff = interpolator.rel_t - (interpolator.t1 - interpolator.t0);
					ASSERT(std::abs(diff) < 1e-8, "sanity check failed - adjusting t by too much: " + std::to_string(diff) +)

					// step
					interpolator.next(hd.planets);

					// reset cur_ts
					interpolator.cur_ts = 0;

					// force the integration time to be equal to the interval begin time
					t = interpolator.t0;

					// now we mark that we just stepped forward in the interpolator - resync uses this
					// to figure out whether to do the "ending encounter resolution" on swift
					starting_lookup_interval = true;
				}

				// this is for the singular case at the very start of the integration, cur_ts is set to 0
				else if (interpolator.cur_ts == 0)
				{
					starting_lookup_interval = true;
				}
				else
				{
					starting_lookup_interval = false;
				}

				// select dt
				cur_dt = interpolator.eff_dt;

				// select the size of the next timestep
				// it can be no more than tbsize, and cannot go past the end of the interval
				// The used tbsize is the minimum of config.tbsize and interpolator interval size.
				cur_tbsize = std::min(config.tbsize, static_cast<uint32_t>(interpolator.n_ts - interpolator.cur_ts));

				// fill the planet logs

				interpolator.fill(hd.planets, cur_tbsize, interpolator.rel_t, cur_dt);

				// advance rel_t
				interpolator.rel_t += static_cast<double>(cur_tbsize) * cur_dt;

				// advance cur_ts
				interpolator.cur_ts += cur_tbsize;

				// since we interpolated planet positions, the planet accelerations were not calculated
				// planet accelerations are needed to calculate h0 which is a term used on the GPU
				// load h0 manually - this doesn't need to happen if using the normal planet integrator
				integrator.load_h0(hd.planets);

				// make sure that we haven't overrun the lookup interval
				ASSERT(interpolator.cur_ts <= interpolator.n_ts, "sanity check fialed - interpolator cur_ts is in an illegal position")
			}
			else
			{
				cur_dt = config.dt;
				cur_tbsize = config.tbsize;

				integrator.integrate_planets_timeblock(hd.planets, cur_tbsize, t, cur_dt);
			}

			// swap new and old logs:
			// the interpolator and integrator both make sure to write to the old logs
			// so swapping the logs here brings the logs into the correct position
			// of course, the GPU integration next chunk is done using the new logs
			swap_logs();

			t_planet = dt_ms(clock, hrclock::now());

			// we only upload the planet log if any particles are going to use the planet log on the GPU
			// i.e. there are particles that could be alive
			if (dd.particle_phase_space().n_alive > 0 || hd.particles.n_encounter() > 0)
			{
				clock = hrclock::now();
				upload_planet_log();
				t_planetup = dt_ms(clock, hrclock::now());
			}
		}

		void Executor::upload_data(size_t begin, size_t length)
		{
			auto &particles = dd.particle_phase_space();
			integrator.upload_data_cuda(htd_stream, begin, length);

			memcpy_htd(particles.r, hd.particles.r(), htd_stream, begin, begin, length);
			cudaStreamSynchronize(htd_stream);
			memcpy_htd(particles.v, hd.particles.v(), htd_stream, begin, begin, length);
			cudaStreamSynchronize(htd_stream);
			memcpy_htd(particles.deathflags, hd.particles.deathflags(), htd_stream, begin, begin, length);
			cudaStreamSynchronize(htd_stream);
			memcpy_htd(particles.deathtime_index, hd.particles.deathtime_index(), htd_stream, begin, begin, length);
			cudaStreamSynchronize(htd_stream);
			memcpy_htd(particles.id, hd.particles.id(), htd_stream, begin, begin, length);
			cudaStreamSynchronize(htd_stream);
		}

		void Executor::add_job(const std::function<void()> &job)
		{
			work.push_back(std::move(job));
		}

		void Executor::alloc_packed()
		{
			size_t n_particle = hd.particles.n();
			size_t size_particle =
				VEC_EL_SIZE(hd.particles.r()) +
				VEC_EL_SIZE(hd.particles.v()) +
				VEC_EL_SIZE(hd.particles.id()) +
				VEC_EL_SIZE(hd.particles.deathflags()) +
				VEC_EL_SIZE(hd.particles.deathtime_index());

			packed_mem_size = size_particle * n_particle;
			cudaMallocHost((uint8_t **)&cpu_packed_mem, packed_mem_size);
			cudaMalloc((uint8_t **)&gpu_packed_mem, packed_mem_size);
		}

		void Executor::download_data(size_t begin, size_t length)
		{
			auto &particles = dd.particle_phase_space();

			// using a batched download instead of broken up downloads
			if (true)
			{
				nvtxRangeId_t id1 = nvtxRangeStartA("pack gpu");

				// PACK ALL GPU DATA
				size_t index_offset = 0;
				size_t length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.r());
				cudaMemcpyAsync(gpu_packed_mem + index_offset, particles.r.data().get() + begin, length_bytes, cudaMemcpyDeviceToDevice, dth_stream);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.v());
				cudaMemcpyAsync(gpu_packed_mem + index_offset, particles.v.data().get() + begin, length_bytes, cudaMemcpyDeviceToDevice, dth_stream);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.id());
				cudaMemcpyAsync(gpu_packed_mem + index_offset, particles.id.data().get() + begin, length_bytes, cudaMemcpyDeviceToDevice, dth_stream);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.deathflags());
				cudaMemcpyAsync(gpu_packed_mem + index_offset, particles.deathflags.data().get() + begin, length_bytes, cudaMemcpyDeviceToDevice, dth_stream);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.deathtime_index());
				cudaMemcpyAsync(gpu_packed_mem + index_offset, particles.deathtime_index.data().get() + begin, length_bytes, cudaMemcpyDeviceToDevice, dth_stream);
				index_offset += length_bytes;

				ASSERT(index_offset <= packed_mem_size, "")
				cudaStreamSynchronize(dth_stream);

				nvtxRangeEnd(id1);

				// DTH
				cudaMemcpyAsync(cpu_packed_mem, gpu_packed_mem, index_offset, cudaMemcpyDeviceToHost, dth_stream);
				cudaStreamSynchronize(dth_stream);

				nvtxRangeId_t id2 = nvtxRangeStartA("unpack cpu");

				// UNPACK CPU DATA

				index_offset = 0;

				length_bytes = length * VEC_EL_SIZE(hd.particles.r());
				memcpy(hd.particles.r().data() + begin, cpu_packed_mem + index_offset, length_bytes);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.v());
				memcpy(hd.particles.v().data() + begin, cpu_packed_mem + index_offset, length_bytes);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.id());
				memcpy(hd.particles.id().data() + begin, cpu_packed_mem + index_offset, length_bytes);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.deathflags());
				memcpy(hd.particles.deathflags().data() + begin, cpu_packed_mem + index_offset, length_bytes);
				index_offset += length_bytes;

				length_bytes = length * VEC_EL_SIZE(hd.particles.deathtime_index());
				memcpy(hd.particles.deathtime_index().data() + begin, cpu_packed_mem + index_offset, length_bytes);
				index_offset += length_bytes;

				nvtxRangeEnd(id2);
			}
			else
			{
				// download only the alive particle data - dead particles are handled in resync()
				// since they're dead, they don't get updated any more so no need to download again

				// note: dead particles DO need to be downloaded when using resync2 so we might as well just download everything
				memcpy_dth(hd.particles.r(), particles.r, dth_stream, begin, begin, length);
				cudaStreamSynchronize(dth_stream);
				memcpy_dth(hd.particles.v(), particles.v, dth_stream, begin, begin, length);
				cudaStreamSynchronize(dth_stream);
				memcpy_dth(hd.particles.id(), particles.id, dth_stream, begin, begin, length);
				cudaStreamSynchronize(dth_stream);
				memcpy_dth(hd.particles.deathflags(), particles.deathflags, dth_stream, begin, begin, length);
				cudaStreamSynchronize(dth_stream);
			}

			// host n_alive includes encounter particles, but not the device n_alive
			hd.particles.n_alive() = dd.particle_phase_space().n_alive + hd.particles.n_encounter();
		}

		void Executor::upload_planet_log()
		{
			// alternate the planet data id, this is to make sure we don't copy into data currently being used on GPU
			dd.planet_data_id++;

			// planet_phase_space uses planet_data_id to figure out which one to get
			auto &planets = dd.planet_phase_space();

			// copy in everything
			// for (auto r_temp : hd.planets.r_log().log){
			// 	temp_log << r_temp << std::endl;
			// }

			memcpy_htd(planets.r_log, hd.planets.r_log().log, htd_stream);
			memcpy_htd(planets.m, hd.planets.m(), htd_stream);
			memcpy_htd(planets.id, hd.planets.id(), htd_stream);

			cudaStreamSynchronize(htd_stream);

			planets.n_alive = hd.planets.n_alive();
			planets.log_len = hd.planets.r_log().len;

			integrator.upload_planet_log_cuda(htd_stream, dd.planet_data_id);
		}

		double Executor::time() const
		{
			auto now = hrclock::now();
			std::chrono::duration<double, std::milli> millis = now - starttime;
			return millis.count() / 60000;
		}

		bool Executor::loop(double *cputimeout, double *gputimeout, double *totaltimeout, size_t *nencounter)
		{
			// At the beginning of the loop the following things should be true:
			// No GPU or CPU processes are running
			// The next timeblock should be loaded on the GPU
			// cur_tbsize refers to the size of that timeblock
			// prev_tbsize refers to the size of the previous timeblock
			// planet logs (old and new are filled with prev_tbsize and cur_tbsize entries, respectively)
			// particles on the GPU have an acceleration loaded (integerator.particle_a)
			// -> this can come from either the previous GPU kernel run, or, the WHIntegrator constructor
			// the time refers to the time before the next GPU kernel, which is also equal to the time at the start of the current planet log
			// it's possible that not all particles have reached the current t, since they might be about to be stepped forward on SWIFT all the way until
			// t + cur_tbsize * dt
			// t = the time at the start of the block that is about to be calculated

			// t_initswift = t_delayswift = t_enc = t_swift = t_writeswift = t_readswift = t_io = t_planet = t_planetup = t_encup = t_sort = t_dl = t_enc2 = t_resync = 0;

			// size_t n_swift = 0;
			auto t_start = hrclock::now();

			size_t n_enc = hd.particles.n_encounter();
			*nencounter = 0;

			if (dd.particle_phase_space().n_alive > 0)
			{
				// if resolving encounters, we need the particle states at the beginning of the chunk
				// so that encounter particles can be rolled back to their initial state

				// manual diagram [1] Backup particle states (Copy particle states to duplicate array)
				if (config.resolve_encounters)
				{
					memcpy_dtd(rollback_state.r, dd.particle_phase_space().r, main_stream);
					memcpy_dtd(rollback_state.v, dd.particle_phase_space().v, main_stream);
					memcpy_dtd(rollback_state.deathflags, dd.particle_phase_space().deathflags, main_stream);
					memcpy_dtd(rollback_state.deathtime_index, dd.particle_phase_space().deathtime_index, main_stream);
					memcpy_dtd(rollback_state.id, dd.particle_phase_space().id, main_stream);

					rollback_state.n_alive = dd.particle_phase_space().n_alive;
					rollback_state.n_total = dd.particle_phase_space().n_total;
					cudaStreamSynchronize(main_stream);
				}

				cudaEventRecord(start_event, main_stream);

				// manual diagram [2] Commence GPU integration
				integrator.integrate_particles_timeblock_cuda(
					main_stream,
					dd.planet_data_id,
					dd.planet_phase_space(),
					dd.particle_phase_space(),
					cur_dt);

				cudaEventRecord(gpu_finish_event, main_stream);
			}

			auto c_start = hrclock::now();

			// do work after GPU starts
			// this is typically all file I/O
			// when encounters are enabled, handle_encounters handles the work vector
			if (config.resolve_encounters && hd.particles.n_encounter() > 0)
			{
				auto cl = hrclock::now();
				*nencounter += handle_encounters(false);
				t_enc = dt_ms(cl, hrclock::now());
			}
			else
			{
				auto cl = hrclock::now();
				nvtxRangeId_t id2 = nvtxRangeStartA("io");

				for (auto &i : work)
					i();
				work.clear();
				t_io = dt_ms(cl, hrclock::now());
				nvtxRangeEnd(id2);

				nvtxRangeId_t id3 = nvtxRangeStartA("planets");

				// The snapshot contains the planet states at the end of the current timechunk (= beginning of next timechunk)
				// this is necessary since update_planets brings all the planets one timechunk forward
				// e.g. if the integration finishes at t=1, update_planets will still bring the planets forward to t=1 + dt
				// so in order to get the correct planetary positions, we need to record the planet positions before they get updated
				hd.planets_snapshot = hd.planets.base;

				// step time forward
				t += cur_dt * static_cast<double>(cur_tbsize);

				// calculate planetary positions for the next chunk - this is REALLY important and also very subtle
				// usually, this would be called at the very end of loop() but
				// we can save some time by doing this here
				// HOWEVER, doing this here means that the planets are now one entire time chunk ahead of the particles
				// for the remainder of loop()
				// it's VERY important to make sure that we're using the correct data

				// IF A FUNCTION IS CALLED AFTER UPDATE_PLANETS: USE OLD DATA (e.g. prev_dt, prev_tbsize, pl.r_log.get<old=true>, etc...)
				update_planets();
				nvtxRangeEnd(id3);
			}
			if (cputimeout)
				*cputimeout = dt_ms(c_start, hrclock::now());

			float gputime = 0;

			// there's nothing to resync if the GPU didn't integrate any particles, i.e. dd.particles.n_alive = 0
			if (dd.particle_phase_space().n_alive > 0)
			{
				// this stuff is queued up after the integration finishes
				auto &particles = dd.particle_phase_space();

				// kill particles in encounters - we do this by setting
				// the dead bit if the encounter bit is also set
				if (!config.resolve_encounters)
				{
					auto it = particles.begin();
					thrust::async::for_each(thrust::cuda::par.on(main_stream), it, it + particles.n_alive, KillEncounterKernel());
				}

				// manual diagram [3] GPU presort
				auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), rollback_state.begin(), integrator.device_begin()));
				presort_index = thrust::partition(thrust::cuda::par.on(main_stream),
												  partition_it, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) -
								partition_it;

				// ------------------------

				// manual diagram [17] WAIT FOR GPU
				cudaStreamSynchronize(main_stream);
				cudaStreamSynchronize(htd_stream);
				cudaEventSynchronize(gpu_finish_event);

				cudaEventElapsedTime(&gputime, start_event, gpu_finish_event);
				if (gputimeout)
					*gputimeout = gputime;

				// resync_every is guaranteed to be 1 if encounters are enabled
				resync_counter++;
				if (resync_counter % config.resync_every == 0)
				{
					auto clock = hrclock::now();
					*nencounter += resync2();
					t_resync = dt_ms(clock, hrclock::now());
					if (cputimeout)
						*cputimeout += t_resync;
				}
			}

			// out_timing << std::left << std::fixed << std::setprecision(3);
			// out_timing << std::setw(8) << t;
			// out_timing << std::setw(8) << n_enc;
			// out_timing << std::setw(8) << t_enc;
			// out_timing << std::setw(8) << t_writeswift;
			// out_timing << std::setw(8) << t_initswift;
			// out_timing << std::setw(8) << t_swift;
			// out_timing << std::setw(8) << t_readswift;
			// out_timing << std::setw(8) << t_delayswift;
			// out_timing << std::setw(8) << t_io;
			// out_timing << std::setw(8) << t_planet;
			// out_timing << std::setw(8) << t_planetup;
			// out_timing << std::setw(8) << t_encup;
			// out_timing << std::setw(8) << t_sort;
			// out_timing << std::setw(8) << t_dl;
			// out_timing << std::setw(8) << t_enc2;
			// out_timing << std::setw(8) << t_resync;
			// out_timing << std::setw(8) << *cputimeout;
			// out_timing << std::setw(8) << *gputimeout;

			// out_timing << std::endl;

			// if not resolving encounters, every time is safe to end on
			// if resolving encounters, only timechunks that end the lookup interval are safe

			if (totaltimeout)
				*totaltimeout = dt_ms(t_start, hrclock::now());
			return starting_lookup_interval || !config.resolve_encounters;
		}

		size_t Executor::resync2()
		{
			// this is a simplified version of the resync algorithm which reads the entire particle arrays, this
			// means that the algorithm doesn't need to match particle ids when reading

			auto &particles = dd.particle_phase_space();

			// previous alive particles
			size_t prev_alive = particles.n_alive;

			auto clock = hrclock::now();

			cudaEventRecord(start_event, sort_stream);

			// partition twice
			// manual diagram [18] Postsort
			if (config.resolve_encounters)
			{
				auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), rollback_state.begin(), integrator.device_begin()));

				// Here particles.n_alive = all particles with non-zero flags, which is actually n_flagged.
				// Notice alive encounter particles coming from SWIFT have their deathflags reset to 0,
				particles.n_alive = thrust::partition(thrust::cuda::par.on(sort_stream),
													  partition_it + presort_index, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) -
									partition_it;

				// the second partition for encounter particles only needs to run between n_alive and prev_alive, since all the alive particles
				// will be pushed to the beginning anyways
				hd.particles.n_encounter() = (thrust::partition(thrust::cuda::par.on(sort_stream),
																partition_it + particles.n_alive, partition_it + prev_alive, DeviceParticleAlivePredicate()) -
											  partition_it) -
											 particles.n_alive;

				cudaStreamSynchronize(sort_stream);
			}
			else
			{
				particles.n_alive = presort_index;
			}

			size_t n_unflagged = particles.n_alive;

			cudaEventRecord(gpu_finish_event, sort_stream);

			cudaEventSynchronize(gpu_finish_event);
			t_sort = dt_ms(clock, hrclock::now());

			float elapsed;
			cudaEventElapsedTime(&elapsed, start_event, gpu_finish_event);
			t_sort = elapsed;

			clock = hrclock::now();
			// copy everything back - n_alive is also copied from device to host
			// manual diagram [19] Download particle data
			download_data(0, particles.n_total);

			t_dl = dt_ms(clock, hrclock::now());

			// output << std::setprecision(8) << t << " " << n_unflagged;

			// for encounter particles, use the rollback data
			// manual diagram [20] Roll-back encounte particles
			if (hd.particles.n_encounter() > 0)
			{
				memcpy_dth(hd.particles.r(), rollback_state.r, dth_stream, particles.n_alive, particles.n_alive, hd.particles.n_encounter());
				cudaStreamSynchronize(dth_stream);
				memcpy_dth(hd.particles.v(), rollback_state.v, dth_stream, particles.n_alive, particles.n_alive, hd.particles.n_encounter());
				cudaStreamSynchronize(dth_stream);
			}

			bool sync_to_end = starting_lookup_interval && hd.particles.n_encounter() > 0;
			size_t nencounter = 0;

			// handle particles that just entered encounter, and partition again
			if (sync_to_end)
			{
				clock = hrclock::now();

				size_t n_encounters = hd.particles.n_encounter();

				// manual diagram [21] Perform end-chunk encounters to sync to interval boundary
				nencounter = handle_encounters(true);

				auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), rollback_state.begin(), integrator.device_begin()));

				// since we only updated encounter particles, we do a partial partition

				particles.n_alive = thrust::partition(thrust::cuda::par.on(sort_stream),
													  partition_it + particles.n_alive - n_encounters, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) -
									partition_it;

				hd.particles.n_encounter() = (thrust::partition(thrust::cuda::par.on(sort_stream),
																partition_it + particles.n_alive, partition_it + prev_alive, DeviceParticleAlivePredicate()) -
											  partition_it) -
											 particles.n_alive;

				// and also only a partial download
				// download_data(n_unflagged, n_encounters);

				cudaStreamSynchronize(main_stream);
				cudaStreamSynchronize(htd_stream);
				cudaStreamSynchronize(dth_stream);
				cudaStreamSynchronize(sort_stream);

				download_data(0, particles.n_total);
				t_enc2 = dt_ms(clock, hrclock::now());
				// output << " " << n_encounters;
			}
			write_encounter(n_unflagged, prev_alive, t - cur_dt * static_cast<double>(cur_tbsize));
			return nencounter;
		}

		size_t Executor::handle_encounters(bool called_from_resync)
		{
			size_t encounter_start = hd.particles.n_alive() - hd.particles.n_encounter();

			// if at the beginning of a lookup interval, don't integrate the previous chunk
			// IGNORE IF CALLED_FROM_RESYNC
			// prev_len = 0 means the previous chunk isn't integrated
			size_t prev_len = prev_tbsize;
			if (starting_lookup_interval && !called_from_resync)
			{
				prev_len = 0;
			}

			// if called from resync, don't integrate the future chunk - this is because
			// the future chunk is in a different lookup interval
			size_t cur_len = called_from_resync ? 0 : cur_tbsize;

			// if called from resync, update_planets will have updated everything to the next chunk
			// already, so we need to use prev_dt
			double which_dt = called_from_resync ? prev_dt : cur_dt;

			if (called_from_resync)
			{
				ASSERT(std::abs(interpolator.t0 - t) < 1e-4, "sanity check failed: end-of-chunk encounter time")
			}

			// if called from resync, t is at the end of the timechunk, otherwise use interpolator relative t MINUS a timeblock because rel_t is the planet time, not the particle time
			double rel_t = called_from_resync ? interpolator.t0 - interpolator.t_m1 : interpolator.rel_t - interpolator.eff_dt * static_cast<double>(cur_tbsize);

			double which_t = called_from_resync ? interpolator.t0 : interpolator.t0 + rel_t;

			nvtxRangeId_t id = nvtxRangeStartA("swift");
			auto clocks = hrclock::now();

			// manual diagram [6]
			auto times1 = swift.begin_integrate(hd.planets, hd.particles, interpolator, called_from_resync, which_t, rel_t, which_dt, prev_len, cur_len);

			if (!called_from_resync)
			{
				t_writeswift = times1.first;
				t_initswift = times1.second;
			}

			// if this was called in the middle of loop, we do the work here while other processes are happening
			if (!called_from_resync)
			{
				auto clock = hrclock::now();
				nvtxRangeId_t id2 = nvtxRangeStartA("io");

				// manual diagram [8]
				for (auto &i : work)
					i();
				work.clear();

				nvtxRangeEnd(id2);
				t_io = dt_ms(clock, hrclock::now());

				// do the planet stuff
				nvtxRangeId_t id3 = nvtxRangeStartA("planets");

				// save the planet locations
				hd.planets_snapshot = hd.planets.base;
				t += cur_dt * static_cast<double>(cur_tbsize);

				// manual diagram [9] [10]
				update_planets();

				nvtxRangeEnd(id3);
			}

			// update encounter particles
			// manual diagram [11]
			auto times = swift.end_integrate(hd.particles);
			nvtxRangeEnd(id);

			if (!called_from_resync)
			{
				t_readswift = times.first;
				t_delayswift = times.second;
				t_swift = dt_ms(clocks, hrclock::now());
			}

			ASSERT(std::isnormal(hd.particles.r()[0].x), "nan")

			// need to calculate particle accelerations for the next timeblock -
			// this is because these particles did not come out of a regular GPU timechunk,
			// so accelerations are outdated

			// load accelerations (the planets already have h0 loaded, so no problem here)

			// manual diagram [12]
			// Done on CPU
			// helio_acc_particles detect all encounter particles and particles leave the inner or outer bounds, and reset their deathflags.
			integrator.helio_acc_particles(
				hd.planets,
				hd.particles,
				encounter_start,
				hd.particles.n_encounter(),
				t,
				hd.planets.r_log().len_old - 1,
				true // use the old log since we just updated the planets
			);

			// Since helio_acc_particles sets deathflags, unset them IFF in encounter since we want the GPU to detect an encounter, delayed
			// for (size_t i = encounter_start; i < encounter_start + hd.particles.n_encounter(); i++)
			// {
			// 	if ((hd.particles.deathflags()[i] & 0xFF) == 0x01)
			// 	{
			// 		hd.particles.deathflags()[i] = 0;
			// 	}
			// }

			clocks = hrclock::now();
			nvtxRangeId_t id4 = nvtxRangeStartA("io");

			// upload the changes to the GPU
			// manual diagram [13]
			upload_data(encounter_start, hd.particles.n_encounter());

			nvtxRangeEnd(id4);

			if (!called_from_resync)
			{
				t_encup = dt_ms(clocks, hrclock::now());
			}

			// set n_alive so that the resync function knows to deal with the particles that we just added back
			dd.particles.n_alive = hd.particles.n_alive();

			return hd.particles.n_encounter();
		}

		void Executor::write_encounter(size_t begin, size_t end, double prev_t)
		{
			for (size_t i = begin; i < end; i++)
			{
				auto &pa = hd.particles;
				uint32_t d_index = pa.deathtime_index()[i];
				double deathtime = static_cast<double>(d_index) * cur_dt + prev_t;
				// t = time at the end of this chunk
				// hd.particles.deathtime_map()[hd.particles.id()[i]] = static_cast<float>(t);
				// deathtime = t;

				discard_output << std::setprecision(17);
				uint16_t deathflag = pa.deathflags()[i];
				// uint16_t tempflag = 65535 - deathflag + 1;
				// std::cout << deathflag << " " << !(deathflag & 0xff) << " " << (deathflag >> 8) << std::endl;
				if (!(deathflag & 0xff) && (deathflag >> 8) !=0)
				{
					if (config.write_encounter_log)
					{
						double closest_approach = pa.deathtime_map()[pa.id()[i]];
						if (closest_approach > 0) *encounter_output << "[E] " << pa.id()[i] << " w/ " << (deathflag >> 8) << " at " << deathtime << " at_c_a_of " << std::setprecision(7) << closest_approach << std::endl;
					}
					pa.deathflags()[i] = 0;
				}
				else if (deathflag == 0x01)
				{
					double helio_dist = sqrt(pa.r()[i].lensq());
					discard_output << "[D] " << pa.id()[i] << " enter_inner_bound_at " << deathtime << " at_helio_dist_of " << std::setprecision(7) << helio_dist << std::endl;
					output << "[D] " << pa.id()[i] << " enter inner bound at " << deathtime << " at helio dist of " << helio_dist << std::endl;
				}
				else if (deathflag == 0x02)
				{
					double helio_dist = sqrt(pa.r()[i].lensq());
					discard_output << "[D] " << pa.id()[i] << " exceed_outer_bound_at " << deathtime << " at_helio_dist_of " << std::setprecision(7) << helio_dist << std::endl;
					output << "[D] " << pa.id()[i] << " exceed outer bound at " << deathtime << " at helio dist of " << helio_dist << std::endl;
				}
				else if (deathflag == 0x04)
				{
					double a, e;
					sr::convert::to_elements(hd.planets.m()[0], pa.r()[i], pa.v()[i], nullptr, &a, &e);
					discard_output << "[D] " << pa.id()[i] << " did_not_converge_at " << deathtime << " at_e_of " << std::setprecision(7) << e << std::endl;
					output << "[D] " << pa.id()[i] << " did not converge at " << deathtime << " at e of " << e << std::endl;
				}
				else if (deathflag == 0x08)
				{
					double helio_dist = sqrt(pa.r()[i].lensq());
					discard_output << "[D] " << pa.id()[i] << " orbit_unbound_at " << deathtime << " at_helio_dist_of " << std::setprecision(7) << helio_dist << std::endl;
					output << "[D] " << pa.id()[i] << " orbit unbound at " << deathtime << " at helio dist of " << helio_dist << std::endl;
				}
				else if (deathflag & 0x80)
				{
					double closest_approach = pa.deathtime_map()[pa.id()[i]];
					discard_output << "[D] " << pa.id()[i] << " absorbed_by " << (deathflag >> 8) << " at " << deathtime << " at_closest_approach_of " << std::setprecision(7) << closest_approach << std::endl;
					output << "[D] " << pa.id()[i] << " absorbed by " << (deathflag >> 8) << " at " << deathtime << " at closest approach of " << closest_approach << std::endl;
				}

				pa.deathtime_map()[pa.id()[i]] = deathtime;
			}
		}

		void Executor::finish()
		{
			cudaStreamSynchronize(main_stream);
			// swift.write_stat(sr::util::joinpath(config.outfolder, "stat.out"));

			for (auto &i : work)
				i();
			work.clear();

			output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive() << std::endl;
		}
	}
}
