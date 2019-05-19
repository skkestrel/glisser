#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <unistd.h>
#include <ctime>

#include "swift.h"
#include "util.cuh"
#include "util.h"
#include "types.h"
#include "executor.cuh"
#include "wh.cuh"
#include "convert.h"

namespace sr
{
namespace exec
{
	using namespace sr::wh;
	using namespace sr::util;
	using namespace sr::convert;
	using namespace sr::data;

	// ExecutorData stores partial data that is copied from the GPU - it's used in resync
	ExecutorData::ExecutorData() { }

	ExecutorData::ExecutorData(size_t n)
	{
		r = v = std::vector<f64_3>(n);
		deathflags = std::vector<uint16_t>(n);
		id = std::vector<uint32_t>(n);
		deathtime_index = std::vector<uint32_t>(n);
	}

	struct DeviceParticleUnflaggedPredicate
	{
		template<typename Tuple>
		__host__ __device__
		bool operator()(const Tuple& args)
		{
			uint16_t flag = thrust::get<2>(thrust::get<0>(args));
			return flag == 0;
		}
	};

	struct DeviceParticleAlivePredicate
	{
		template<typename Tuple>
		__host__ __device__
		bool operator()(const Tuple& args)
		{
			uint16_t flag = thrust::get<2>(thrust::get<0>(args));
			return (flag & 0xFE) == 0;
		}
	};

	struct KillEncounterKernel
	{
		KillEncounterKernel() { }

		template<typename Tuple>
		__host__ __device__
		void operator()(Tuple args) const
		{
			uint16_t flags = thrust::get<2>(args);
			if ((flags & 0x01) == 0x01)
			{
				flags |= 0x80;
			}
			thrust::get<2>(args) = flags;
		}
	};


	Executor::Executor(HostData& _hd, DeviceData& _dd, const Configuration& _config, std::ostream& out)
		: hd(_hd), dd(_dd), output(out), config(_config), resync_counter(0) { }

	void Executor::assert_true(bool cond, std::string msg)
	{
		if (!cond)
		{
			throw std::runtime_error(msg);
		}
	}

	void Executor::init()
	{
		// glisse only supports helio
		sr::convert::to_helio(hd);

		// setup integrator - this gives all particles an acceleration and also detects initial encounters
		integrator = sr::wh::WHCudaIntegrator(hd.planets, hd.particles, config);

		if (config.interp_planets)
		{
			// setup interpolator
			interpolator = sr::interp::Interpolator(config, hd.planets, config.planet_history_file);

			// step interp forward until we're in the correct interval
			while (interpolator.t1 <= t)
			{
				interpolator.next(hd.planets);
			}

			// make sure we call next at least once, otherwise it means that the lookup file started past
			// the current time

			// interpolator.t0 is initialized to +inf, but calling next once makes it finite
			// so if t0 is not infinity, next was called at least once
			assert_true(!std::isinf(interpolator.t0), "lookup file starts past the given Initial-Time");

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
		}

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
		cudaStreamCreate(&htd_stream);
		cudaStreamCreate(&par_stream);
		cudaStreamCreate(&dth_stream);

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


		// partition alive
		hd.particles.stable_partition_unflagged(0, hd.particles.n());

		// now, n_alive = n_unflagged
		size_t n_alive = hd.particles.n_alive();

		// bring the encounter particles to the beginning
		hd.particles.stable_partition_alive(n_alive, hd.particles.n() - n_alive);
		// now, n_alive = total alive (enc + nonenc)

		// don't need to upload if no particles
		if (hd.particles.n() > 0)
		{
			upload_data(0, hd.particles.n());
			dd.particle_phase_space().n_alive = hd.particles.n_alive();
		}

		// download data right after uploading -
		// i forgot why we need this, but I think it's just so the GPU data can be debugged
		download_data();

		starttime = std::chrono::high_resolution_clock::now();

		output << "       Starting simulation.       " << std::endl << std::endl;

		if (encounter_output)
		{
			*encounter_output << std::setprecision(17);
		}

		// upload planet data before the first timechunk
		update_planets();
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

		if (config.interp_planets)
		{
			// cur_ts = n_ts means that we've reached the end of the lookup interval
			// need to step to the next one
			if (interpolator.cur_ts == interpolator.n_ts)
			{
				// sanity check - make sure that rel_t = t1 - t0
				double diff = interpolator.rel_t - (interpolator.t1 - interpolator.t0);
				assert_true(std::abs(diff) < 1e-8, "sanity check failed - adjusting t by too much: " + std::to_string(diff));

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
			assert_true(interpolator.cur_ts <= interpolator.n_ts, "sanity check fialed - interpolator cur_ts is in an illegal position");

			std::cout << "tbsize = " << cur_tbsize << ", curstep = " << interpolator.cur_ts << "/" << interpolator.n_ts << ", prevstep = " << prev_tbsize << std::endl;
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

		// we only upload the planet log if any particles are going to use the planet log on the GPU
		// i.e. there are particles that could be alive
		if (dd.particle_phase_space().n_alive > 0 || hd.particles.n_encounter() > 0)
		{
			upload_planet_log();
		}
	}

	void Executor::upload_data(size_t begin, size_t length)
	{
		auto& particles = dd.particle_phase_space();
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

	void Executor::add_job(const std::function<void()>& job)
	{
		work.push_back(std::move(job));
	}

	void Executor::download_data()
	{
		auto& particles = dd.particle_phase_space();

		// download only the alive particle data - dead particles are handled in resync()
		// since they're dead, they don't get updated any more so no need to download again

		// note: dead particles DO need to be downloaded when using resync2 so we might as well just download everything
		memcpy_dth(hd.particles.r(), particles.r, dth_stream, 0, 0, particles.n_total);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.v(), particles.v, dth_stream, 0, 0, particles.n_total);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.id(), particles.id, dth_stream, 0, 0, particles.n_total);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.deathflags(), particles.deathflags, dth_stream, 0, 0, particles.n_total);
		cudaStreamSynchronize(dth_stream);

		// host n_alive includes encounter particles, but not the device n_alive
		hd.particles.n_alive() = dd.particle_phase_space().n_alive + hd.particles.n_encounter();
	}

	void Executor::upload_planet_log()
	{
		// alternate the planet data id, this is to make sure we don't copy into data currently being used on GPU
		dd.planet_data_id++;

		// planet_phase_space uses planet_data_id to figure out which one to get
		auto& planets = dd.planet_phase_space();

		// copy in everything
		memcpy_htd(planets.r_log, hd.planets.r_log().log, htd_stream);
		memcpy_htd(planets.m, hd.planets.m(), htd_stream);
		memcpy_htd(planets.id, hd.planets.id(), htd_stream);

		cudaStreamSynchronize(htd_stream);

		planets.log_len = hd.planets.r_log().len;

		integrator.upload_planet_log_cuda(htd_stream, dd.planet_data_id);
	}

	double Executor::time() const
	{
		auto now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> millis = now - starttime;
		return millis.count() / 60000;
	}

	void Executor::handle_encounters(bool called_from_resync)
	{
		std::cout << "encounter start, hd n = " << hd.particles.n_alive() << ", dd = " << dd.particle_phase_space().n_alive << std::endl;

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
		if (called_from_resync)
		{
			cur_len = 0;
		}

		// if called from resync, update_planets will have updated everything to the next chunk
		// already, so we need to use prev_dt
		double which_dt = called_from_resync ? prev_dt : cur_dt;

		if (called_from_resync)
		{
			assert_true(std::abs(interpolator.t0 - t) < 1e-2, "sanity check failed: end-of-chunk encounter time");
		}

		// do the thing
		sr::swift::SwiftEncounterIntegrator enc(config, t, which_dt, prev_len, cur_len);

		enc.begin_integrate(hd.planets, hd.particles, interpolator, called_from_resync);
		
		// if this was called in the middle of loop, we do the work here while other processes are happening
		if (!called_from_resync)
		{
			for (auto& i : work) i();
			work.clear();
		}
		
		// update encounter particles
		enc.end_integrate(hd.particles);

		// whether to use the old log or not depends on whether we called from resync
		size_t which_timestep_index = called_from_resync ? hd.planets.r_log().len_old : hd.planets.r_log().len;

		// need to calculate particle accelerations for the next timeblock -
		// this is because these particles did not come out of a regular GPU timechunk,
		// so accelerations are outdated

		// load accelerations (the planets already have h0 loaded, so no problem here)
		integrator.helio_acc_particles(
			hd.planets,
			hd.particles,
			encounter_start,
			hd.particles.n_encounter(),
			t + static_cast<double>(cur_len) * which_dt,
			which_timestep_index,
			called_from_resync // use the old log if called from resync, otherwise use the new log
		);
			
		// upload the changes to the GPU
		// no need to sort the particles here, resync will do all the sorting
		upload_data(encounter_start, hd.particles.n_encounter());

		// set n_alive so that the resync function knows to deal with the particles that we just added back
		// TODO this messes up the particle count
		dd.particles.n_alive = hd.particles.n_alive();

		std::cout << "encounter end, hd n = " << hd.particles.n_alive() << ", dd = " << dd.particle_phase_space().n_alive << std::endl;
	}

	void Executor::loop(double* cputimeout, double* gputimeout)
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

		std::clock_t c_start = std::clock();

		if (dd.particle_phase_space().n_alive > 0)
		{
			// if resolving encounters, we need the particle states at the beginning of the chunk
			// so that encounter particles can be rolled back to their initial state
			if (config.resolve_encounters)
			{
				memcpy_dtd(rollback_state.r, dd.particle_phase_space().r, main_stream);
				memcpy_dtd(rollback_state.v, dd.particle_phase_space().v, main_stream);
				memcpy_dtd(rollback_state.deathflags, dd.particle_phase_space().deathflags, main_stream);
				memcpy_dtd(rollback_state.deathtime_index, dd.particle_phase_space().deathtime_index, main_stream);
				memcpy_dtd(rollback_state.id, dd.particle_phase_space().id, main_stream);

				rollback_state.n_alive = dd.particle_phase_space().n_alive;
				rollback_state.n_total = dd.particle_phase_space().n_total;
			}

			cudaEventRecord(start_event, main_stream);

			// in order to integrate the particles on GPU, the particle accelerations must be set.
			// typically the accelerations are set by the previous timeblock
			// but in the case of the first timeblock, or when recovering from a close encounter, it needs to be set manually...
			integrator.integrate_particles_timeblock_cuda(
				main_stream,
				dd.planet_data_id,
				dd.planet_phase_space(),
				dd.particle_phase_space(),
				cur_dt
			);

			cudaEventRecord(gpu_finish_event, main_stream);
		}

		// do work after GPU starts
		// this is typically all file I/O
		// when encounters are enabled, handle_encounters handles the work vector
		if (config.resolve_encounters && hd.particles.n_encounter() > 0)
		{
			handle_encounters(false);
		}
		else
		{
			for (auto& i : work) i();
			work.clear();
		}

		// The snapshot contains the planet states at the end of the previous timechunk (= beginning of current timechunk)
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

		float cputime = static_cast<float>(std::clock() - c_start) / CLOCKS_PER_SEC * 1000;
		if (cputimeout) *cputimeout = cputime;

		// there's nothing to resync if the GPU didn't integrate any particles, i.e. dd.particles.n_alive = 0
		if (dd.particle_phase_space().n_alive > 0)
		{
			cudaStreamSynchronize(htd_stream);
			cudaEventSynchronize(gpu_finish_event);

			float gputime;
			cudaEventElapsedTime(&gputime, start_event, gpu_finish_event);
			if (gputimeout) *gputimeout = gputime;

			// resync_every is guaranteed to be 1 if encounters are enabled
			resync_counter++;
			if (resync_counter % config.resync_every == 0)
			{
				resync2();
			}
		}
	}

	void Executor::resync2()
	{
		std::cout << "resync begin, hd n = " << hd.particles.n_alive() << ", dd = " << dd.particle_phase_space().n_alive << std::endl;

		// this is a simplified version of the resync algorithm which reads the entire particle arrays, this
		// means that the algorithm doesn't need to match particle ids when reading

		auto& particles = dd.particle_phase_space();
		size_t prev_alive = particles.n_alive;

		// kill particles in encounters
		if (!config.resolve_encounters)
		{
			auto it = particles.begin();
			thrust::for_each(thrust::cuda::par.on(main_stream), it, it + particles.n_alive, KillEncounterKernel());
		}

		// partition twice
		if (config.resolve_encounters)
		{
			auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), rollback_state.begin(), integrator.device_begin()));

			particles.n_alive = thrust::stable_partition(thrust::cuda::par.on(main_stream),
					partition_it, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) - partition_it;
			
			// the second partition for encounter particles only needs to run between n_alive and prev_alive, since all the alive particles
			// will be pushed to the beginning anyways
			hd.particles.n_encounter() = (thrust::stable_partition(thrust::cuda::par.on(main_stream),
					partition_it + particles.n_alive, partition_it + prev_alive, DeviceParticleAlivePredicate()) - partition_it) - particles.n_alive;
		}
		else
		{
			auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), integrator.device_begin()));
			particles.n_alive = thrust::stable_partition(thrust::cuda::par.on(main_stream),
					partition_it, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) - partition_it;
		}

		// copy everything back - n_alive is also copied from device to host
		download_data();

		// set the deathtime for dead particles - let's set the encounter particles deathtimes too, just to show when they entered encounter
		// here t refers to the ending time of the timechunk
		for (size_t i = particles.n_alive; i < prev_alive; i++)
		{
			hd.particles.deathtime()[i] = static_cast<float>(t);
		}

		// for encounter particles, use the rollback data
		if (hd.particles.n_encounter() > 0)
		{
			memcpy_dth(hd.particles.r(), particles.r, dth_stream, 0, particles.n_alive, hd.particles.n_encounter());
			cudaStreamSynchronize(dth_stream);
			memcpy_dth(hd.particles.v(), particles.v, dth_stream, 0, particles.n_alive, hd.particles.n_encounter());
			cudaStreamSynchronize(dth_stream);
		}

		// handle particles that just entered encounter, and partition again
		if (starting_lookup_interval && hd.particles.n_encounter() > 0)
		{
			handle_encounters(true);

			auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), rollback_state.begin(), integrator.device_begin()));

			particles.n_alive = thrust::stable_partition(thrust::cuda::par.on(main_stream),
					partition_it, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) - partition_it;
			
			hd.particles.n_encounter() = (thrust::stable_partition(thrust::cuda::par.on(main_stream),
					partition_it + particles.n_alive, partition_it + prev_alive, DeviceParticleAlivePredicate()) - partition_it) - particles.n_alive;
		}
	}

	void Executor::resync()
	{
		// the responsibility of resync is to bring updated data from the GPU back to the CPU, EXCEPT for particles which are still alive
		// and running on the GPU

		// this function allows the CPU to know about particles in close encounters
		// using a series of partitions
		
		auto& particles = dd.particle_phase_space();

		size_t prev_alive = particles.n_alive;


		// first partition on the GPU - this moves "gpu-marked" dead particles to the end
		// gpu-marked dead particles means: the particle is really dead (e.g. out of bounds)
		// or came close to a planet (which MAY trigger an encounter if config.resolve_encounters = 1)
		// note - this first partition doesn't differentiate between close encounter and dead particles
		if (config.resolve_encounters)
		{
			// if handling encounters, partition the rollback as well
			auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), rollback_state.begin(), integrator.device_begin()));

			// particles.n_alive is equal to the number of particles that the GPU thinks is alive: n_part - n_encounter - n_dead
			particles.n_alive = thrust::stable_partition(thrust::cuda::par.on(main_stream),
					partition_it, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) - partition_it;
		}
		else
		{
			auto partition_it = thrust::make_zip_iterator(thrust::make_tuple(particles.begin(), integrator.device_begin()));
			particles.n_alive = thrust::stable_partition(thrust::cuda::par.on(main_stream),
					partition_it, partition_it + particles.n_alive, DeviceParticleUnflaggedPredicate()) - partition_it;
		}

		cudaStreamSynchronize(main_stream);

		// the difference between the previously number of alive particles and the current
		// - i.e. this is the number of new "dead" or encounter particles that have been generated
		size_t diff = prev_alive - particles.n_alive;

		// copying dead particles to CPU
		// executor data shows us particle data for diffed particles
		// - we don't need to know about every single particle on the GPU, only the changed ones
		ExecutorData exdata = ExecutorData(diff);

		memcpy_dth(exdata.r, particles.r, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(exdata.v, particles.v, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(exdata.id, particles.id, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(exdata.deathtime_index, particles.deathtime_index, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(exdata.deathflags, particles.deathflags, dth_stream, 0, particles.n_alive, diff);
		cudaStreamSynchronize(dth_stream);

		ExecutorData rollback_exdata;
	       
		// if resolving encounters, we also copy back the rollback data
		if (config.resolve_encounters)
		{
			rollback_exdata = ExecutorData(diff);
			memcpy_dth(rollback_exdata.r, rollback_state.r, dth_stream, 0, particles.n_alive, diff);
			cudaStreamSynchronize(dth_stream);
			memcpy_dth(rollback_exdata.v, rollback_state.v, dth_stream, 0, particles.n_alive, diff);
			cudaStreamSynchronize(dth_stream);
			memcpy_dth(rollback_exdata.id, rollback_state.id, dth_stream, 0, particles.n_alive, diff);
			cudaStreamSynchronize(dth_stream);
			memcpy_dth(rollback_exdata.deathtime_index, rollback_state.deathtime_index, dth_stream, 0, particles.n_alive, diff);
			cudaStreamSynchronize(dth_stream);
			memcpy_dth(rollback_exdata.deathflags, rollback_state.deathflags, dth_stream, 0, particles.n_alive, diff);
			cudaStreamSynchronize(dth_stream);
		}

		// rewrite deathflags in the case of an encounter - in any case, iterate over all of the dead particles
		for (size_t i = 0; i < diff; i++)
		{
			if ((exdata.deathflags[i] & 0x00FF) == 0x0001)
			{
				if (config.resolve_encounters)
				{
					// if encounter handling, do nothing, keep the flag though
				}
				else
				{
					// If encounters are not being dealt with, kill the particle!
					exdata.deathflags[i] = 0x0080;
				}
			}

			if (exdata.deathflags[i] & 0x0004)
			{
				output << "Warning: simulation did not converge on particle " << exdata.id[i] << std::endl;
			}

			if (exdata.deathflags[i] & 0x0002)
			{
				output << "Warning: particle " << exdata.id[i] << " OOB" << std::endl;
			}
		}

		// now we need to partition the "gpu-marked" dead particles -
		// a "dead" particle on the GPU might still be alive if it's in an encounter
		// so we partition one more time to move encounter particles to the beginning of the dead zone, but the end of the alive zone
		std::unique_ptr<std::vector<size_t>> ed_indices;
		sr::data::stable_partition_alive_indices(exdata.deathflags, 0, diff, &ed_indices);

		gather(exdata.r, *ed_indices, 0, diff);
		gather(exdata.v, *ed_indices, 0, diff);
		gather(exdata.id, *ed_indices, 0, diff);
		gather(exdata.deathflags, *ed_indices, 0, diff);
		gather(exdata.deathtime_index, *ed_indices, 0, diff);

		// gotta mirror the operation 

		if (config.resolve_encounters)
		{
			gather(rollback_exdata.r, *ed_indices, 0, diff);
			gather(rollback_exdata.v, *ed_indices, 0, diff);
			gather(rollback_exdata.id, *ed_indices, 0, diff);
			gather(rollback_exdata.deathflags, *ed_indices, 0, diff);
			gather(rollback_exdata.deathtime_index, *ed_indices, 0, diff);
		}

		// map from particle index to array index
		std::unordered_map<size_t, size_t> indices;
		for (size_t i = 0; i < prev_alive; i++)
		{
			indices[hd.particles.id()[i]] = i;
		}

		// now we update the CPU particles information with the information from the GPU
		// this is essentially copying the GPU information into CPU
		for (size_t i = 0; i < diff; i++)
		{
			size_t index = indices[exdata.id[i]];

			hd.particles.deathflags()[index] = exdata.deathflags[i];
			hd.particles.deathtime_index()[index] = exdata.deathtime_index[i];

			// decide where to copy from: either roll back to beginning of chunk, or take the values right at death on the GPU
			ExecutorData* datasource;
			if (config.resolve_encounters)
			{
				// if encounter:
				if ((exdata.deathflags[i] & 0x00FF) == 0x0001)
				{
					datasource = &rollback_exdata;
				}
				else // otherwise the particle is really dead: OOB or non converge, etc
				{
					datasource = &exdata;
				}
			}
			else
			{
				datasource = &exdata;
			}

			hd.particles.r()[index] = datasource->r[i];
			hd.particles.v()[index] = datasource->v[i];

			// if the particle is dead or in an encounter, we can record the time that it happened
			if (exdata.deathflags[i])
			{
				hd.particles.deathtime()[index] = static_cast<float>(
					t - prev_dt * static_cast<double>(config.tbsize - exdata.deathtime_index[i])
				);
			}
		}

		// now, we try to match the GPU ordering
		// first, partition by unflagged
		auto gather_indices = hd.particles.stable_partition_unflagged(0, prev_alive);
		integrator.gather_particles(*gather_indices, 0, prev_alive);

		// then, partition by alive - this will bring encounter particles, if any, between the two regions
		// in the end:
		// the array will look like
		// alive, encounter, dead
		if (config.resolve_encounters)
		{
			auto gather_indices2 = hd.particles.stable_partition_alive(0, prev_alive);
			integrator.gather_particles(*gather_indices2, 0, prev_alive);
		}

		// hd.particles.n_alive is updated after each gather
		// thus, it will point to the total number of alive + encoutner particles after gathering
		hd.particles.n_encounter() = hd.particles.n_alive() - particles.n_alive;


		// check that the CPU ids match the GPU ids
		Vu32 prev_ids(hd.particles.id().begin(), hd.particles.id().end());
		memcpy_dth(hd.particles.id(), particles.id, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);

		assert_true(prev_ids == hd.particles.id(), "particle ID mismatch occurred");

		// if the last timechunk of a lookup interval was reached, finish the encounters: run swift NOW!
		// this allows everything to be up-to-date at the end of a lookup-interval, this makes it safe to 
		// write particle histories

		// starting_lookup_interval refers to the timechunk right after this, so if the next timechunk is starting,
		// this one is ending
		if (starting_lookup_interval && hd.particles.n_encounter() > 0)
		{
			handle_encounters(true);
		}

		add_job([prev_alive, this]()
			{
				if (encounter_output)
				{
					// loop over dead particles only - 
					for (size_t i = hd.particles.n_alive(); i < prev_alive; i++)
					{
						// state
						*encounter_output << hd.particles.r()[i] << std::endl;
						*encounter_output << hd.particles.v()[i] << std::endl;

						// id, deathflags, deathtime
						*encounter_output << hd.particles.id()[i] << " "
							<< hd.particles.deathflags()[i] << " "
							// deathtime comes from 
							<< t - prev_dt * static_cast<double>(config.tbsize - hd.particles.deathtime_index()[i]) << " death"
							<< std::endl;

						// output planet states too
						*encounter_output << hd.planets.n_alive_old() << std::endl;
						*encounter_output << hd.planets.m()[0] << std::endl;
						*encounter_output << f64_3(0) << std::endl;
						*encounter_output << f64_3(0) << std::endl;
						*encounter_output << hd.planets.id()[0] << std::endl;
						for (size_t j = 1; j < hd.planets.n_alive(); j++)
						{
							*encounter_output << hd.planets.m()[j] << std::endl;
							*encounter_output << hd.planets.r_log().old[hd.particles.deathtime_index()[i] * (hd.planets.n() - 1) + j - 1] << std::endl;
							*encounter_output << hd.planets.v_log().old[hd.particles.deathtime_index()[i] * (hd.planets.n() - 1) + j - 1] << std::endl;
							*encounter_output << hd.planets.id()[j] << std::endl;
						}
					}

					// TODO output encounter particles as well!

					*encounter_output << std::flush;
				}
			});
	}


	void Executor::finish()
	{
		cudaStreamSynchronize(main_stream);

		for (auto& i : work) i();
		work.clear();

		resync();

		for (auto& i : work) i();
		work.clear();

		output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive() << std::endl;
	}
}
}
