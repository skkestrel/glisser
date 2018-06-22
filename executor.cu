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

Executor::Executor(HostData& hd, DeviceData& dd, std::ostream& out) : hd(hd), dd(dd), print_every(10), print_counter(0), tbsize(128), ce_factor(1),
	output(out), timing_output(nullptr), resolve_encounters(false) { }

void Executor::init()
{
	to_helio(hd);

	calculate_planet_metrics(hd.planets, &e_0, nullptr);

	output << std::setprecision(7);
	output << "e_0 (planets) = " << e_0 << std::endl;
	output << "n_particle = " << hd.particles.n << std::endl;
	output << "n_particle_alive = " << hd.particles.n_alive << std::endl;
	output << "==================================" << std::endl;
	output << "Running for half a time step.     " << std::endl;

	initialize(hd.planets, hd.particles);

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

	if (timing_output)
	{
		*timing_output << std::setprecision(17);
	}
	if (discard_output)
	{
		*discard_output << std::setprecision(17);
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
			step_planets(hd.planets, t, i, dt / ce_factor);
			// take the planet positions at the end of every timestep

			if (i % ce_factor == ce_factor - 1)
			{
				size_t slow_index = i / ce_factor;
				auto fast_begin = hd.planets.r_log.begin() + i * (hd.planets.n - 1);
				std::copy(fast_begin, fast_begin + (hd.planets.n - 1), hd.planets.r_log_slow.begin() + slow_index * (hd.planets.n - 1));

				hd.planets.h0_log_slow[i / ce_factor] = hd.planets.h0_log[i];
			}
		}
	}
	else
	{
		for (size_t i = 0; i < tbsize; i++)
		{
			step_planets(hd.planets, t, i, dt);
		}

		std::copy(hd.planets.h0_log.begin(), hd.planets.h0_log.end(), hd.planets.h0_log_slow.begin());
		std::copy(hd.planets.r_log.begin(), hd.planets.r_log.end(), hd.planets.r_log_slow.begin());
	}

	upload_planet_log();
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

	if (prev_ids != hd.particles.id)
	{
		output << "WARNING! ID MISMATCH! WARNING!" << std::endl;
		throw std::exception();
	}

	hd.particles.n_alive = dd.particle_phase_space().n_alive;

	/*
	// zip will crash the program

	auto iterator = thrust::make_zip_iterator(thrust::make_tuple(
				ps.r.begin(),
				ps.v.begin(),
				ps.deathflags.begin(), ps.deathtime.begin(), ps.id.begin()));
	thrust::copy(thrust::cuda::par.on(stream),
			iterator,
			iterator + n,
			thrust::make_zip_iterator(thrust::make_tuple(
					hd.particles.r.begin(), hd.particles.v.begin(),
					hd.particles.deathflags.begin(), hd.deathtime.begin(), hd.id.begin())));
	 */
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

	for (auto& i : work)
	{
		i();
	}
	work.clear();

	hd.planets_snapshot = HostPlanetSnapshot(hd.planets);

	t += dt * tbsize;
	step_and_upload_planets();
	cudaStreamSynchronize(htd_stream);

	if (print_counter % print_every == 0)
	{
		double e;
		calculate_planet_metrics(hd.planets, &e, nullptr);

		double elapsed = time();
		double total = elapsed * (t_f - t_0) / (t - t_0);

		output << "t=" << t << " (" << elapsed / total * 100 << "% " << elapsed << "m elapsed, " << total << "m total " << total - elapsed << "m remain)" << std::endl;
		output << "Error = " << (e - e_0) / e_0 * 100 << ", " << dd.particle_phase_space().n_alive << " particles remaining" << std::endl;
	}
	print_counter++;

	if (timing_output)
	{
		double e_;
		f64_3 l_;
		calculate_planet_metrics(hd.planets, &e_, &l_);
	
		*timing_output << "ep " << e_ << std::endl;
		*timing_output << "lp " << l_.x << " " << l_.y << " " << l_.z << std::endl;
	}

	for (size_t i = hd.particles.n_alive; i < hd.particles.n_alive + hd.particles.n_encounter; i++)
	{
		// step_particle(..)
	}
	// memcy_htd(hd.n_alive, hd.n_encounter)
	// hd.n_alive = partition(...)
	// dd.n_alive = hd.n_alive

	cudaStreamSynchronize(main_stream);
	resync();
}

void Executor::resync()
{
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

	add_job([prev_alive, diff, this]()
		{
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

				std::ostream* output_stream = nullptr;

				if ((ed.deathflags[i] & 0x0001) && resolve_encounters)
				{
					output_stream = discard_output;
				}
				else
				{
					output_stream = discard_output;
					hd.particles.deathtime[index] = t + dt * ed.deathtime_index[i];
				}

				if (output_stream)
				{
					*output_stream << ed.r[i] << std::endl;
					*output_stream << ed.v[i] << std::endl;
					*output_stream << hd.particles.deathtime[index] << " " << ed.deathflags[i] << " " << ed.id[i] << std::endl;
					*output_stream << hd.planets.n - 1 << std::endl;
					for (size_t j = 1; j < hd.planets.n; j++)
					{
						*output_stream << hd.planets.m[j] << std::endl;
						*output_stream << hd.planets.r_log_slow[ed.deathtime_index[i] * (hd.planets.n - 1) + j - 1] << std::endl;
					}
				}
			}

			hd.particles.stable_partition_alive();
		});
}


void Executor::finish()
{
	cudaStreamSynchronize(main_stream);
	resync();
	download_data();

	output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive << std::endl;
}
