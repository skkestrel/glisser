#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <algorithm>

#include "types.h"
#include "executor.cuh"
#include "wh.cuh"
#include "convert.h"
#include "wh.h"

ExecutorData::ExecutorData() { }
ExecutorData::ExecutorData(size_t n)
{
	r = v = std::vector<f64_3>(n);
	id = deathtime_index = std::vector<uint32_t>(n);
	deathflags = std::vector<uint16_t>(n);
}

template<typename T>
cudaError_t memcpy_dth(std::vector<T>& dest, const thrust::device_vector<T>& src, cudaStream_t& stream, size_t destbegin = 0, size_t srcbegin = 0, size_t len = static_cast<uint32_t>(-1))
{
	if (len == static_cast<uint32_t>(-1))
	{
		len = src.size();
	}
	if (dest.size() < destbegin + len)
	{
		throw std::exception();
	}

	return cudaMemcpyAsync(dest.data() + destbegin, src.data().get() + srcbegin, len * sizeof(T), cudaMemcpyDeviceToHost, stream);
}

template<typename T>
cudaError_t memcpy_htd(thrust::device_vector<T>& dest, const std::vector<T>& src, cudaStream_t& stream, size_t destbegin = 0, size_t srcbegin = 0, size_t len = static_cast<uint32_t>(-1))
{
	if (len == static_cast<uint32_t>(-1))
	{
		len = src.size();
	}
	if (dest.size() < destbegin + len)
	{
		throw std::exception();
	}

	return cudaMemcpyAsync(dest.data().get() + destbegin, src.data() + srcbegin, len * sizeof(T), cudaMemcpyHostToDevice, stream);
}

struct DeviceParticleUnflaggedPredicate
{
	template<typename Tuple>
	__host__ __device__
	bool operator()(const Tuple& args)
	{
		uint8_t flag = thrust::get<3>(args);
		return flag == 0;
	}
};

Executor::Executor(HostData& hd, DeviceData& dd, std::ostream& out) : hd(hd), dd(dd), tbsize(128), ce_factor(1),
	output(out), resolve_encounters(false) { }

void Executor::init()
{
	to_helio(hd);
	initialize(hd.planets, hd.particles, wh_alloc);
	calculate_planet_metrics(hd.planets, wh_alloc, &e_0, nullptr);

	output << std::setprecision(7);
	output << "e_0 (planets) = " << e_0 << std::endl;
	output << "n_particle = " << hd.particles.n << std::endl;
	output << "n_particle_alive = " << hd.particles.n_alive << std::endl;
	output << "==================================" << std::endl;
	output << "Sending initial conditions to GPU." << std::endl;

	cudaStreamCreate(&main_stream);
	cudaStreamCreate(&htd_stream);
	cudaStreamCreate(&par_stream);
	cudaStreamCreate(&dth_stream);

	upload_data();
	output << "n_particle_alive = " << dd.particle_phase_space().n_alive << std::endl;

	resync();
	download_data();

	starttime = std::chrono::high_resolution_clock::now();

	output << "       Starting simulation.       " << std::endl << std::endl;

	if (encounter_output)
	{
		*encounter_output << std::setprecision(17);
	}

	step_and_upload_planets();
	if (!resolve_encounters)
	{
		ce_factor = 1;
	}
}

void Executor::step_and_upload_planets()
{
	if (resolve_encounters)
	{

		for (size_t i = 0; i < tbsize * ce_factor; i++)
		{
			step_planets(hd.planets, wh_alloc, t, i, dt / ce_factor);
			// take the planet positions at the end of every timestep

			if (i % ce_factor == ce_factor - 1)
			{
				size_t slow_index = i / ce_factor;

				auto fast_begin = hd.planets.r_log.begin() + i * (hd.planets.n - 1);
				std::copy(fast_begin, fast_begin + (hd.planets.n - 1), hd.planets.r_log_slow.begin() + slow_index * (hd.planets.n - 1));

				fast_begin = hd.planets.v_log.begin() + i * (hd.planets.n - 1);
				std::copy(fast_begin, fast_begin + (hd.planets.n - 1), hd.planets.v_log_slow.begin() + slow_index * (hd.planets.n - 1));

				hd.planets.h0_log_slow[i / ce_factor] = hd.planets.h0_log[i];
			}
		}
	}
	else
	{
		for (size_t i = 0; i < tbsize; i++)
		{
			step_planets(hd.planets, wh_alloc, t, i, dt);
		}

		std::copy(hd.planets.h0_log.begin(), hd.planets.h0_log.end(), hd.planets.h0_log_slow.begin());
		std::copy(hd.planets.r_log.begin(), hd.planets.r_log.end(), hd.planets.r_log_slow.begin());
		std::copy(hd.planets.v_log.begin(), hd.planets.v_log.end(), hd.planets.v_log_slow.begin());
	}

	// We only upload the planet log if any particles are going to use the planet log on the GPU
	// Cases where the planet log is not used by the particles:
	// - There are no particles alive on the GPU, AND there are no particles in close encounters on the CPU
	// since the particles that survive close encounters can make it to the GPU at the end of this timestep
	// and thus the next planet chunk will be required
	if (dd.particle_phase_space().n_alive == 0 && hd.particles.n_encounter == 0)
	{
		upload_planet_log();
	}
}

void Executor::upload_data()
{
	dd.particles0 = DeviceParticlePhaseSpace(hd.particles.n);
	dd.particles1 = DeviceParticlePhaseSpace(hd.particles.n);

	dd.planets0 = DevicePlanetPhaseSpace(hd.planets.n, tbsize);
	dd.planets1 = DevicePlanetPhaseSpace(hd.planets.n, tbsize);

	dd.planet_data_id = 0;
	dd.particle_data_id = 0;

	auto& particles = dd.particle_phase_space();

	particles.n_alive = hd.particles.n_alive;

	memcpy_htd(particles.r, hd.particles.r, htd_stream);
	cudaStreamSynchronize(htd_stream);
	memcpy_htd(particles.v, hd.particles.v, htd_stream);
	cudaStreamSynchronize(htd_stream);
	memcpy_htd(particles.a, hd.particles.a, htd_stream);
	cudaStreamSynchronize(htd_stream);
	memcpy_htd(particles.deathflags, hd.particles.deathflags, htd_stream);
	cudaStreamSynchronize(htd_stream);
	memcpy_htd(particles.id, hd.particles.id, htd_stream);
	cudaStreamSynchronize(htd_stream);

	memcpy_htd(dd.planet_phase_space().m, hd.planets.m, htd_stream);
	cudaStreamSynchronize(htd_stream);

	dd.planet_data_id++;
	memcpy_htd(dd.planet_phase_space().m, hd.planets.m, htd_stream);
	cudaStreamSynchronize(htd_stream);
}

void Executor::add_job(const std::function<void()>& job)
{
	work.push_back(std::move(job));
}

void Executor::download_data()
{
	auto& particles = dd.particle_phase_space();

	Vu32 prev_ids(hd.particles.id.begin(), hd.particles.id.end());

	memcpy_dth(hd.particles.r, particles.r, dth_stream);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(hd.particles.v, particles.v, dth_stream);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(hd.particles.a, particles.a, dth_stream);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(hd.particles.id, particles.id, dth_stream);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(hd.particles.deathflags, particles.deathflags, dth_stream);
	cudaStreamSynchronize(dth_stream);

	// This should NEVER happen. I think this is a recoverable 
	// error, by swapping particle indices on the host, but that sounds annoying...
	if (prev_ids != hd.particles.id)
	{
		output << "WARNING! ID MISMATCH! WARNING!" << std::endl;
		throw std::exception();
	}

	hd.particles.n_alive = dd.particle_phase_space().n_alive;
}

void Executor::upload_planet_log()
{
	dd.planet_data_id++;
	auto& planets = dd.planet_phase_space();

	memcpy_htd(planets.h0_log, hd.planets.h0_log_slow, htd_stream);
	cudaStreamSynchronize(htd_stream);
	memcpy_htd(planets.r_log, hd.planets.r_log_slow, htd_stream);
	cudaStreamSynchronize(htd_stream);
}


double Executor::time() const
{
	auto now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> millis = now - starttime;
	return millis.count() / 60000;
}

void Executor::loop()
{
	step_particles_cuda(main_stream, dd.planet_phase_space(), dd.particle_phase_space(), tbsize, dt);

	// The queued work should begin RIGHT after the CUDA call
	for (auto& i : work)
	{
		i();
	}
	work.clear();

	// The snapshot contains the planet states at the end of the previous timestep - 
	// consider removing this? We can use hd.planets.*_log_old[-1] to replicate this functionality
	hd.planets_snapshot = HostPlanetSnapshot(hd.planets);

	// The OLD logs are required by the close encounter handler
	std::swap(hd.planets.r_log, hd.planets.r_log_old);
	std::swap(hd.planets.v_log, hd.planets.v_log_old);
	std::swap(hd.planets.r_log_slow, hd.planets.r_log_slow_old);
	std::swap(hd.planets.v_log_slow, hd.planets.v_log_slow_old);

	std::swap(hd.planets.h0_log, hd.planets.h0_log_old);
	std::swap(hd.planets.h0_log_slow, hd.planets.h0_log_slow_old);

	t += dt * tbsize;
	step_and_upload_planets();
	cudaStreamSynchronize(htd_stream);

	for (size_t i = hd.particles.n_alive; i < hd.particles.n_alive + hd.particles.n_encounter; i++)
	{
		// step_particle(..)
	}
	// hd.n_alive = hd.particles.stable_partition_alive(...)
	// memcy_htd(particles.n_alive, hd.n_encounter)
	// dd.n_alive = hd.n_alive

	cudaStreamSynchronize(main_stream);
	resync();
}

void Executor::resync()
{
	// There's nothing to resync if all the particles on the device are dead!
	// Although dd.particles.n_alive can be out of sync with dd.particles.deathflags before
	// resync() is called, this is safe:
	// - The MVS kernel does not revive particles, so resync() will never INCREASE n_alive
	// - dd.particles.n_alive is adjusted by the close encounter handler just BEFORE this call
	if (dd.particle_phase_space().n_alive == 0) return;

	auto& particles = dd.particle_phase_space();
	size_t prev_alive = particles.n_alive;

	particles.n_alive = thrust::stable_partition(thrust::cuda::par.on(par_stream), particles.begin(), particles.begin() + particles.n_alive, DeviceParticleUnflaggedPredicate())
		- particles.begin();
	cudaStreamSynchronize(par_stream);

	size_t diff = prev_alive - particles.n_alive;

	ed = ExecutorData(diff);

	memcpy_dth(ed.r, particles.r, dth_stream, 0, particles.n_alive, diff);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(ed.v, particles.v, dth_stream, 0, particles.n_alive, diff);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(ed.id, particles.id, dth_stream, 0, particles.n_alive, diff);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(ed.deathtime_index, particles.deathtime_index, dth_stream, 0, particles.n_alive, diff);
	cudaStreamSynchronize(dth_stream);
	memcpy_dth(ed.deathflags, particles.deathflags, dth_stream, 0, particles.n_alive, diff);
	cudaStreamSynchronize(dth_stream);

	std::unordered_map<size_t, size_t> indices;
	for (size_t i = 0; i < prev_alive; i++)
	{
		indices[hd.particles.id[i]] = i;
	}

	for (size_t i = 0; i < diff; i++)
	{
		size_t index = indices[ed.id[i]];
		hd.particles.r[index] = ed.r[i];
		hd.particles.v[index] = ed.v[i];
		hd.particles.deathflags[index] = ed.deathflags[i];

		if ((ed.deathflags[i] & 0x0001) && resolve_encounters)
		{
			// TODO do something here... maybe clear the death bit?
		}
		else
		{
			hd.particles.deathtime[index] = t - dt * (tbsize - ed.deathtime_index[i]);
		}
	}

	hd.particles.stable_partition_alive();
	hd.particles.n_encounter = hd.particles.n_alive - particles.n_alive;

	add_job([diff, this]()
		{
			for (size_t i = 0; i < diff; i++)
			{
				*encounter_output << ed.r[i] << std::endl;
				*encounter_output << ed.v[i] << std::endl;
				*encounter_output << ed.id[i] << " " << ed.deathflags[i] << " " << t - dt * (tbsize - ed.deathtime_index[i]) << std::endl;
				*encounter_output << hd.planets.n_alive << std::endl;

				*encounter_output << hd.planets.m[0] << std::endl;
				*encounter_output << f64_3(0) << std::endl;
				*encounter_output << f64_3(0) << std::endl;
				*encounter_output << hd.planets.id[0] << std::endl;
				for (size_t j = 1; j < hd.planets.n_alive; j++)
				{
					*encounter_output << hd.planets.m[j] << std::endl;
					*encounter_output << hd.planets.r_log_slow[ed.deathtime_index[i] * (hd.planets.n - 1) + j - 1] << std::endl;
					*encounter_output << hd.planets.v_log_slow[ed.deathtime_index[i] * (hd.planets.n - 1) + j - 1] << std::endl;
					*encounter_output << hd.planets.id[i] << std::endl;
				}
			}
		});
}


void Executor::finish()
{
	cudaStreamSynchronize(main_stream);
	resync();
	download_data();

	output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive << std::endl;
}
