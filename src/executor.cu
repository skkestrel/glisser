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
		assert_true(config.use_gpu, "USE_GPU is not allowed to be 0 on this version of GLISSER");

		to_helio(hd);

		integrator = sr::wh::WHCudaIntegrator(hd.planets, hd.particles, config);
		if (config.interp_planets)
		{
			interpolator = sr::interp::Interpolator(config, hd.planets, config.planet_history_file);
			interpolator.next(hd.planets);
		}

		calculate_planet_metrics(hd.planets, &e_0, nullptr);

		output << std::setprecision(7);
		output << "e_0 (planets) = " << e_0 << std::endl;
		output << "n_particle = " << hd.particles.n() << std::endl;
		output << "n_particle_alive = " << hd.particles.n_alive() << std::endl;
		output << "==================================" << std::endl;
		output << "Sending initial conditions to GPU." << std::endl;

		cudaStreamCreate(&main_stream);
		cudaStreamCreate(&htd_stream);
		cudaStreamCreate(&par_stream);
		cudaStreamCreate(&dth_stream);

		cudaEventCreate(&start_event);
		cudaEventCreate(&gpu_finish_event);

		rollback_state = DeviceParticlePhaseSpace(hd.particles.n());
		dd.particles = DeviceParticlePhaseSpace(hd.particles.n());

		dd.planets0 = DevicePlanetPhaseSpace(hd.planets.n(), config.tbsize);
		dd.planets1 = DevicePlanetPhaseSpace(hd.planets.n(), config.tbsize);
		dd.planet_data_id = 0;

		memcpy_htd(dd.planet_phase_space().m, hd.planets.m(), htd_stream);
		cudaStreamSynchronize(htd_stream);
		dd.planet_data_id++;
		memcpy_htd(dd.planet_phase_space().m, hd.planets.m(), htd_stream);
		cudaStreamSynchronize(htd_stream);

		if (hd.particles.n() > 0)
		{
			upload_data(0, hd.particles.n());
		}

		download_data();

		starttime = std::chrono::high_resolution_clock::now();

		output << "       Starting simulation.       " << std::endl << std::endl;

		if (encounter_output)
		{
			*encounter_output << std::setprecision(17);
		}

		update_planets();
	}

	void Executor::swap_logs()
	{
		hd.planets.swap_logs();
		integrator.swap_logs();
	}

	void Executor::update_planets()
	{
		prev_tbsize = cur_tbsize;

		if (config.interp_planets)
		{
			if (interpolator.cur_ts == interpolator.n_ts)
			{
				interpolator.next(hd.planets);
				assert_true(std::abs(t - interpolator.t0) < 1e-6, "adjusting t by too much");
				t = interpolator.t0;
			}

			// select the size of the next timestep
			cur_tbsize = std::min(config.tbsize, static_cast<uint32_t>(interpolator.n_ts - interpolator.cur_ts));

			interpolator.fill(hd.planets, cur_tbsize, t, current_dt());
			// need to fill h0 manually because planet acceleration was not calculated
			integrator.load_h0(hd.planets);

			interpolator.cur_ts += cur_tbsize;

			assert_true(interpolator.cur_ts <= interpolator.n_ts, "interpolator timestep is in an illegal position");
		}
		else
		{
			cur_tbsize = config.tbsize;
			integrator.integrate_planets_timeblock(hd.planets, cur_tbsize, t, current_dt());
		}

		swap_logs();

		// We only upload the planet log if any particles are going to use the planet log on the GPU
		// Cases where the planet log is not used by the particles:
		// - There are no particles alive on the GPU, AND there are no particles in close encounters on the CPU
		// since the particles that survive close encounters can make it to the GPU at the end of this timestep
		// and thus the next planet chunk will be required

		if (dd.particle_phase_space().n_alive > 0 || hd.particles.n_encounter() > 0)
		{
			upload_planet_log();
		}
	}

	void Executor::upload_data(size_t begin, size_t length)
	{
		auto& particles = dd.particle_phase_space();
		particles.n_alive = hd.particles.n_alive();
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

	void Executor::download_data(bool ignore_errors)
	{
		auto& particles = dd.particle_phase_space();

		Vu32 prev_ids(hd.particles.id().begin(), hd.particles.id().end());

		// download only the alive particle data - the dead particle data is handled in resync()
		// also, the particles that the GPU thinks are dead might be in encounter
		memcpy_dth(hd.particles.r(), particles.r, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.v(), particles.v, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.id(), particles.id, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);
		memcpy_dth(hd.particles.deathflags(), particles.deathflags, dth_stream, 0, 0, particles.n_alive);
		cudaStreamSynchronize(dth_stream);

		// This should NEVER happen. I think this is a recoverable 
		// error, by swapping particle indices on the host, but that sounds annoying...
		if (prev_ids != hd.particles.id())
		{
			/*
			for (int i = 0; i < prev_ids.size(); i++)
			{
				if (prev_ids[i] != hd.particles.id()[i])
				{
					std::cout << i << " " << prev_ids[i] << " " << hd.particles.id()[i] << std::endl;
				}
			}
			*/
			output << "WARNING! ID MISMATCH! WARNING!" << std::endl;

			assert_true(ignore_errors, "ID mismatch occurred when trying to download GPU data.");
		}

		hd.particles.n_alive() = dd.particle_phase_space().n_alive;
	}

	double Executor::current_dt()
	{
		return config.interp_planets ? interpolator.eff_dt : config.dt;
	}

	void Executor::upload_planet_log()
	{
		dd.planet_data_id++;
		auto& planets = dd.planet_phase_space();

		memcpy_htd(planets.r_log, hd.planets.r_log().slow, htd_stream);
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


		std::clock_t c_start = std::clock();
		// t = the time at the start of the block that is about to be calculated
		std::thread cpu_thread;
		
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

			// in order to integrate the particles on GPU, the particle accelerations must be set
			// typically the accelerations are set by the previous timeblock
			// but in the case of the first timeblock, or when recovering from a close encounter, it needs to be set manually...
			// the constructor of WHIntegrator does this to fulfill the roll in the first timeblock
			integrator.integrate_particles_timeblock_cuda(
				main_stream,
				dd.planet_data_id,
				dd.planet_phase_space(),
				dd.particle_phase_space(),
				current_dt()
			);
			cudaEventRecord(gpu_finish_event, main_stream);
		}

		// The queued work should begin RIGHT after the CUDA call
		if (config.resolve_encounters && hd.particles.n_encounter() > 0)
		{
			size_t encounter_start = hd.particles.n_alive() - hd.particles.n_encounter();

			if (config.enable_swift)
			{
				sr::swift::SwiftEncounterIntegrator enc(config, t, current_dt(), prev_tbsize, cur_tbsize);

				enc.begin_integrate(hd.planets, hd.particles, interpolator);
				
				// work happens here - this can happen while close encounters are happening
				for (auto& i : work) i();
				work.clear();
				

				// TODO: Particles that finished here were started from the previous timechunk's beginning,
				// this means that if there was a history output event after the previous timechunk,
				// the encounter particles will be missing. We need to fill the array
				enc.end_integrate(hd.particles);
			}
			else
			{
				for (auto& i : work) i();
				work.clear();

				for (size_t i = encounter_start; i < hd.particles.n_alive(); i++)
				{
					integrator.integrate_encounter_particle_catchup(
							hd.planets,
							hd.particles,
							i,
							0,
							t - current_dt() * static_cast<double>(prev_tbsize),
							current_dt()
					);
				}
			}

			auto gather_indices = hd.particles.stable_partition_alive(encounter_start, hd.particles.n_encounter());
			integrator.gather_particles(*gather_indices, encounter_start, hd.particles.n_encounter());
			upload_data(encounter_start, hd.particles.n_encounter());

			// Fill deathtime index with 0 so that the continuation will work
			// This makes it so that the particles which remain in an encounter at the end of this encounter will continue
			// integrating from where it left off rather than some random point
			// well, we don't use CPU close encounter handling anyways...
			thrust::fill(thrust::cuda::par.on(htd_stream), dd.particles.deathtime_index.begin() + encounter_start,
					dd.particles.deathtime_index.begin() + encounter_start + hd.particles.n_encounter(), 0);

			// need to calculate particle accelerations for the next timeblock -
			// this is because these particles did not come out of a regular GPU timechunk,
			// so accelerations are incorrect right now
#warning TODO
			// integrator.base.helio_acc_particles<false, true>(hd.planets, hd.particles, encounter_start, hd.particles.n_encounter(), 0, 0);
		}
		else
		{
			for (auto& i : work) i();
			work.clear();
		}


		// The snapshot contains the planet states at the end of the previous timestep - 
		// consider removing this? We can use hd.planets.*_log_old()[-1] to replicate this functionality
		hd.planets_snapshot = hd.planets.base;

		// step time forward since the GPU has already started
		t += current_dt() * static_cast<double>(cur_tbsize);

		update_planets();


		float cputime = static_cast<float>(std::clock() - c_start) / CLOCKS_PER_SEC * 1000;
		if (cputimeout) *cputimeout = cputime;

		if (dd.particle_phase_space().n_alive > 0)
		{
			cudaStreamSynchronize(htd_stream);
			cudaEventSynchronize(gpu_finish_event);

			float gputime;
			cudaEventElapsedTime(&gputime, start_event, gpu_finish_event);
			if (gputimeout) *gputimeout = gputime;

			// There's nothing to resync if all the particles on the device are dead!
			// Although dd.particles.n_alive can be out of sync with dd.particles.deathflags before
			// resync() is called, this is safe:
			// - The MVS kernel does not revive particles, so resync() will never INCREASE n_alive
			// - dd.particles.n_alive is adjusted by the close encounter handler just BEFORE this call


			resync_counter++;

			if (resync_counter % config.resync_every == 0)
			{
				resync();
			}
		}
	}

	void Executor::resync()
	{
		auto& particles = dd.particle_phase_space();

		size_t prev_alive = particles.n_alive;

		// first partition on the GPU - this moves "gpu-marked" dead particles to the end
		// also partition the rollback if doing close encounters

		// note - the first partition doesn't differentiate between close encounter and dead particles

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
					exdata.deathflags[i] |= 0x0080;
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
					t - current_dt() * static_cast<double>(config.tbsize - exdata.deathtime_index[i])
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


		size_t encounter_start = particles.n_alive;

		// 

		// TODO: if the last timechunk of a lookup interval was reached, finish the encounters: run swift NOW! (for half the regular time, i.e. 1/2 * 2timechunk size = 1 timechunk)

		add_job([encounter_start, diff, &exdata, this]()
			{
				if (encounter_output)
				{
					for (size_t i = hd.particles.n_encounter(); i < diff; i++)
					{
						*encounter_output << hd.particles.r()[encounter_start + i] << std::endl;
						*encounter_output << hd.particles.v()[encounter_start + i] << std::endl;
						*encounter_output << hd.particles.id()[encounter_start + i] << " "
							<< hd.particles.deathflags()[encounter_start + i] << " "
							<< t - current_dt() * static_cast<double>(config.tbsize - exdata.deathtime_index[i]) << " death"
							<< std::endl;
						*encounter_output << hd.planets.n_alive() << std::endl;

						*encounter_output << hd.planets.m()[0] << std::endl;
						*encounter_output << f64_3(0) << std::endl;
						*encounter_output << f64_3(0) << std::endl;
						*encounter_output << hd.planets.id()[0] << std::endl;
						for (size_t j = 1; j < hd.planets.n_alive(); j++)
						{
							*encounter_output << hd.planets.m()[j] << std::endl;
							*encounter_output << hd.planets.r_log().slow[exdata.deathtime_index[i] * (hd.planets.n() - 1) + j - 1] << std::endl;
							*encounter_output << hd.planets.v_log().slow[exdata.deathtime_index[i] * (hd.planets.n() - 1) + j - 1] << std::endl;
							*encounter_output << hd.planets.id()[j] << std::endl;
						}
					}

					*encounter_output << std::flush;
				}
			});
	}


	void Executor::finish()
	{
		// TODO: upon finishing the simulation there is a half step of encounters that needs to be done.
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
