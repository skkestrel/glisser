/*************************************************************

.---. .            .         .   .-,--.                        
\___  |- ,-. ,-. ,-| . . ,-. |-   `|__/ ,-. .  , ,-. ,-. . ,-. 
    \ |  ,-| |   | | | | `-. |    )| \  |-' | /  |-' |   | |-' 
`---' `' `-^ '   `-^ `-^ `-' `'   `'  ` `-' `'   `-' '   ' `-' 

*************************************************************/


// Stardust Reverie v0:
// Close encounter handling +
// timestep chunking.

// TODO:
// Figure out where to store particle death information
// How to deal with close encounters
// Profile CPU version of code


#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <thread>
#include <atomic>
#include <cstdint>
#include <limits>

#include <thrust/system/cuda/execution_policy.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>

const int SIMPL_DEGREE = 4;
const int TIMEBLOCK_SIZE = 1024;

using Hvf64 = thrust::host_vector<double>;
using Hvu32 = thrust::host_vector<uint32_t>;
using Hvu8 = thrust::host_vector<uint8_t>;
using Hvf32 = thrust::host_vector<float>;

using Dvf64 = thrust::device_vector<double>;
using Dvu32 = thrust::device_vector<uint32_t>;
using Dvu8 = thrust::device_vector<uint8_t>;
using Dvf32 = thrust::device_vector<float>;

struct f64_3
{
	double x, y, z;

	__host__ __device__
	f64_3() { }
};

using Dvf64_3 = thrust::device_vector<f64_3>;
using Hvf64_3 = thrust::host_vector<f64_3>;

void load_coefs(Hvf64& coef, double dt)
{
	coef = Hvf64(SIMPL_DEGREE * 2);

	double curt2 = std::pow(2., 1. / 3.);

	coef[0] = coef[3] = (2. + curt2 + 1. / curt2) / 6.; 
	coef[1] = coef[2] = (1. - curt2 - 1. / curt2) / 6.;
	coef[5] = coef[7] = 1. / (2. - curt2);
	coef[6] = 1. / (1. - std::pow(2., 2. / 3.));
	coef[4] = 0.;

	for (size_t i = 0; i < SIMPL_DEGREE * 2; i++)
	{
		coef[i] *= dt;
	}
}

struct SIA4Kernel
{
	const double* coefdt;
	const double* m_planet;
	const f64_3* r_planet_log;
	size_t n_planet;

	SIA4Kernel(const double* coefdt, const f64_3* r_planet_log, const double* m_planet, size_t n_planet)
	{
		this->n_planet = n_planet;
		this->coefdt = coefdt;
		this->r_planet_log = r_planet_log;
		this->m_planet = m_planet;
	}

	__host__ __device__
	f64_3 newton3d(const f64_3& r, uint8_t& flags, uint32_t& deathtime, const f64_3* planet_log_offset) const
	{
		uint32_t n_planet_local = this->n_planet;
		const double* m_planet_index = this->m_planet;

		f64_3 a;
		a.x = a.y = a.z = 0;

		for (uint32_t i = 0; i < n_planet_local; i++)
		{
			double dist;
			f64_3 d;

			{
				f64_3 planets = *(planet_log_offset + i);

				{
					d.x = r.x - planets.x;
					dist = d.x * d.x;
				}

				{
					d.y = r.y - planets.y;
					dist += d.y * d.y;
				}

				{
					d.z = r.z - planets.z;
					dist += d.z * d.z;
				}
			}

			flags = flags || (dist < 0.5 * 0.5);

#ifdef USE_RSQRT
			double mul = *m_planet_index * rsqrt(dist) / dist;
#else
			double mul = *m_planet_index / sqrt(dist) / dist;
#endif
			a.x -= d.x * mul;
			a.y -= d.y * mul;
			a.z -= d.z * mul;

			m_planet_index++;
		}
		return a;
	}

	template<typename Tuple>
	__host__ __device__
	void operator()(Tuple args) const
	{ 
		f64_3 r = thrust::get<0>(args);
		f64_3 v = thrust::get<1>(args);
		uint8_t flags = thrust::get<2>(args);
		uint32_t deathtime = thrust::get<3>(args);

		uint32_t log_offset = 0;
		uint32_t n_planet_local = n_planet;

		const f64_3* r_planet_log_local = this->r_planet_log;

		f64_3 running_r = r;
		f64_3 running_v = v;

		for (uint32_t step = 0; step < TIMEBLOCK_SIZE; step++)
		{
			double adt = coefdt[0];
			double bdt;


			r.x += adt * v.x;
			r.y += adt * v.y;
			r.z += adt * v.z;

#pragma unroll
			for (uint32_t coef = 1; coef < SIMPL_DEGREE; coef++)
			{
				adt = coefdt[coef];
				bdt = coefdt[coef + SIMPL_DEGREE];

				f64_3 a = newton3d(r, flags, deathtime, r_planet_log_local + log_offset);

				{
					v.x += bdt * a.x;
				}

				{
					v.y += bdt * a.y;
				}

				{
					v.z += bdt * a.z;
				}

				r.x += adt * v.x;
				r.y += adt * v.y;
				r.z += adt * v.z;

				log_offset += n_planet_local;
			}

			flags = flags || (r.x * r.x + r.y * r.y + r.z * r.z > 100 * 100);

			deathtime = deathtime * flags + step * !flags;
			running_r.x = running_r.x * flags + r.x * !flags;
			running_r.y = running_r.y * flags + r.y * !flags;
			running_r.z = running_r.z * flags + r.z * !flags;
			running_v.x = running_v.x * flags + v.x * !flags;
			running_v.y = running_v.y * flags + v.y * !flags;
			running_v.z = running_v.z * flags + v.z * !flags;
		}
		
		thrust::get<0>(args) = running_r;
		thrust::get<1>(args) = running_v;
		thrust::get<2>(args) = flags;
		thrust::get<3>(args) = deathtime;
	}
};
 

void sia4_planet(const Hvf64& coefdt, Hvf64_3& r, Hvf64_3& v, Hvf64& m, Hvf64_3& r_log);
void newton3d_planet(Hvf64_3& f, const Hvf64_3& r, const Hvf64& m);

double energy_planet(const Hvf64_3& r, const Hvf64_3& v, const Hvf64& m);

struct DevicePhaseSpace
{
	Dvf64_3 r, v;
	Dvu8 flags;
	Dvu32 deathtime;

	// TODO shove this into flags later
	Dvu32 id;
};

struct DeviceData
{
	// double buffer for HtD transfer of planet locations
	Dvf64_3 r_planet_log0, r_planet_log1;

	Dvf64 m_planet;

	// double buffer for killing dead particles
	DevicePhaseSpace ps0, ps1;

	Dvu32 gather_indices;

	size_t log_buffer_id, phase_space_id;
	size_t n_part_alive;

	Dvf64 coefdt;

	DeviceData() { }
};

struct HostData
{
	Hvf64_3 r_planet, v_planet;
	Hvf64_3 r_planet_log;

	Hvf64 m_planet;

	Hvf64_3 r_part, v_part;

	Hvf64 coefdt;

	Hvu8 part_flags;
	Hvu32 deathtime;
	Hvu32 id;

	double t, dt, t_f;
	size_t n_part, n_planet;

	HostData() { }
};

bool load_data(HostData& hd, std::string plin, std::string icsin, size_t max_particle = 0, bool readmomenta = true)
{
	std::ifstream plinfile(plin), icsinfile(icsin);

	plinfile >> hd.n_planet;

	hd.r_planet = Hvf64_3(hd.n_planet);
	hd.v_planet = Hvf64_3(hd.n_planet);
	hd.m_planet = Hvf64(hd.n_planet);

	size_t lsize = hd.n_planet * (SIMPL_DEGREE - 1) * TIMEBLOCK_SIZE;
	hd.r_planet_log = Hvf64_3(lsize);

	for (size_t i = 0; i < hd.n_planet; i++)
	{
		plinfile >> hd.m_planet[i];
		plinfile >> hd.r_planet[i].x >> hd.r_planet[i].y >> hd.r_planet[i].z;
		plinfile >> hd.v_planet[i].x >> hd.v_planet[i].y >> hd.v_planet[i].z;

		if (readmomenta)
		{
			hd.v_planet[i].x /= hd.m_planet[i];
			hd.v_planet[i].y /= hd.m_planet[i];
			hd.v_planet[i].z /= hd.m_planet[i];
		}
	}

	icsinfile >> hd.n_part;
	if (max_particle > 0) hd.n_part = std::min(hd.n_part, max_particle);

	hd.r_part = Hvf64_3(hd.n_part);
	hd.v_part = Hvf64_3(hd.n_part);
	hd.part_flags = Hvu8(hd.n_part);
	hd.deathtime = Hvu32(hd.n_part);
	hd.id = Hvu32(hd.n_part);

	for (size_t i = 0; i < hd.n_part; i++)
	{
		icsinfile >> hd.r_part[i].x >> hd.r_part[i].y >> hd.r_part[i].z;
		icsinfile >> hd.v_part[i].x >> hd.v_part[i].y >> hd.v_part[i].z;

		std::string s;
		icsinfile >> s;
		if (!isdigit(s[0]))
		{
			icsinfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			hd.deathtime[i] = 0;
			hd.id[i] = i;
			hd.part_flags[i] = 0;
		}
		else
		{
			hd.deathtime[i] = std::stod(s);
			icsinfile >> hd.part_flags[i] >> hd.id[i];
		}
	}

	load_coefs(hd.coefdt, hd.dt);
	return false;
}

void save_data(const HostData& hd, const DeviceData& dd, std::string plout, std::string icsout, std::string infoout)
{
	std::ofstream ploutfile(plout), icsoutfile(icsout), infooutfile(infoout);

	ploutfile << hd.n_planet << std::endl;
	ploutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.n_planet; i++)
	{
		ploutfile << hd.m_planet[i] << std::endl;
		ploutfile << hd.r_planet[i].x << " " << hd.r_planet[i].y << " " << hd.r_planet[i].z << std::endl;
		ploutfile << hd.v_planet[i].x << " " << hd.v_planet[i].y << " " << hd.v_planet[i].z << std::endl;
	}

	icsoutfile << hd.n_part << std::endl;
	icsoutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.n_part; i++)
	{
		icsoutfile << hd.r_part[i].x << " " << hd.r_part[i].y << " " << hd.r_part[i].z << std::endl;
		icsoutfile << hd.v_part[i].x << " " << hd.v_part[i].y << " " << hd.v_part[i].z << std::endl;
		icsoutfile << hd.deathtime[i] << " " << hd.part_flags[i] << " " << hd.id[i] << std::endl;
	}

	infooutfile << std::setprecision(17);
	infooutfile << hd.t << " " << dd.n_part_alive << std::endl;
}

void transfer_data(const HostData& hd, DeviceData& dd)
{
	dd.log_buffer_id = 0;
	dd.phase_space_id = 0;

	dd.r_planet_log0 = Dvf64_3(hd.r_planet_log.size());
	dd.r_planet_log1 = Dvf64_3(hd.r_planet_log.size());
	dd.m_planet = Dvf64(hd.n_planet);

	dd.ps0.r = dd.ps0.v = Dvf64_3(hd.n_part);
	dd.ps0.flags = Dvu8(hd.n_part);
	dd.ps0.deathtime = Dvu32(hd.n_part);
	dd.ps0.id = Dvu32(hd.n_part);

	dd.ps1.r = dd.ps1.v = Dvf64_3(hd.n_part);
	dd.ps1.flags = Dvu8(hd.n_part);
	dd.ps1.deathtime = Dvu32(hd.n_part);
	dd.ps1.id = Dvu32(hd.n_part);

	dd.gather_indices = Dvu32(hd.n_part);

	dd.n_part_alive = 0;
	for (size_t i = 0; i < hd.n_part; i++)
	{
		if (~hd.part_flags[i] & 0x0001)
		{
			dd.n_part_alive++;
		}
	}

	dd.coefdt = Dvf64(hd.coefdt.size());

	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;
	thrust::copy(hd.r_part.begin(), hd.r_part.end(), ps.r.begin());
	thrust::copy(hd.v_part.begin(), hd.v_part.end(), ps.v.begin());

	thrust::copy(hd.m_planet.begin(), hd.m_planet.end(), dd.m_planet.begin());
	thrust::copy(hd.coefdt.begin(), hd.coefdt.end(), dd.coefdt.begin());

	thrust::copy(hd.id.begin(), hd.id.end(), ps.id.begin());
	thrust::copy(hd.deathtime.begin(), hd.deathtime.end(), ps.deathtime.begin());
	thrust::copy(hd.part_flags.begin(), hd.part_flags.end(), ps.flags.begin());
}

void recover_data(HostData& hd, const DeviceData& dd, cudaStream_t& stream)
{
	const DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;

	size_t n = hd.n_part;

	/*
	auto iterator = thrust::make_zip_iterator(thrust::make_tuple(
				ps.r.begin(),
				ps.v.begin(),
				ps.flags.begin(), ps.deathtime.begin(), ps.id.begin()));
	thrust::copy(thrust::cuda::par.on(stream),
			iterator,
			iterator + n,
			thrust::make_zip_iterator(thrust::make_tuple(
					hd.r_part.begin(), hd.v_part.begin(),
					hd.part_flags.begin(), hd.deathtime.begin(), hd.id.begin())));
	 */

	thrust::copy(thrust::cuda::par.on(stream), ps.r.begin(), ps.r.begin() + n, hd.r_part.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.v.begin(), ps.v.begin() + n, hd.v_part.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.flags.begin(), ps.flags.begin() + n, hd.part_flags.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.deathtime.begin(), ps.deathtime.begin() + n, hd.deathtime.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.id.begin(), ps.id.begin() + n, hd.id.begin());
}

void convert_to_barycentric(HostData& hd)
{
	f64_3 r;
	f64_3 v;
	r.x = r.y = r.z = 0;
	v.x = v.y = v.z = 0;
	double totalm = 0;

	for (size_t i = 0; i < hd.n_planet; i++)
	{
		r.x += hd.r_planet[i].x * hd.m_planet[i];
		r.y += hd.r_planet[i].y * hd.m_planet[i];
		r.z += hd.r_planet[i].z * hd.m_planet[i];

		v.x += hd.v_planet[i].x * hd.m_planet[i];
		v.y += hd.v_planet[i].y * hd.m_planet[i];
		v.z += hd.v_planet[i].z * hd.m_planet[i];

		totalm += hd.m_planet[i];
	}

	r.x /= totalm;
	r.y /= totalm;
	r.z /= totalm;
	v.x /= totalm;
	v.y /= totalm;
	v.z /= totalm;

	for (size_t i = 0; i < hd.n_planet; i++)
	{
		hd.r_planet[i].x -= r.x;
		hd.r_planet[i].y -= r.y;
		hd.r_planet[i].z -= r.z;

		hd.v_planet[i].x -= v.x;
		hd.v_planet[i].y -= v.y;
		hd.v_planet[i].z -= v.z;
	}

	for (size_t i = 0; i < hd.n_part; i++)
	{
		hd.r_part[i].x -= r.x;
		hd.r_part[i].y -= r.y;
		hd.r_part[i].z -= r.z;

		hd.v_part[i].x -= v.x;
		hd.v_part[i].y -= v.y;
		hd.v_part[i].z -= v.z;
	}
}

size_t prune(cudaStream_t& main_stream, cudaStream_t& work_stream, HostData& hd, DeviceData& dd)
{
	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;
	DevicePhaseSpace& other_ps = dd.phase_space_id % 2 ? dd.ps1 : dd.ps0;

	thrust::copy(thrust::cuda::par.on(work_stream), ps.flags.begin(), ps.flags.begin() + dd.n_part_alive, hd.part_flags.begin());
	Hvu32 indices = Hvu32(hd.n_part);

	// start placing alive particles at the front and dead particles at the back
	uint32_t front_counter = 0;
	uint32_t back_counter = dd.n_part_alive - 1;

	for (size_t i = 0; i < dd.n_part_alive; i++)
	{
		if (hd.part_flags[i] & 0x0001) // dead particle
		{
			indices[back_counter--] = i;
		}
		else
		{
			indices[front_counter++] = i;
		}
	}
	for (size_t i = dd.n_part_alive; i < hd.n_part; i++)
	{
		indices[i] = i;
	}

	cudaStreamSynchronize(work_stream);
	thrust::copy(thrust::cuda::par.on(work_stream), indices.begin(), indices.end(), dd.gather_indices.begin());

	size_t pruned = dd.n_part_alive - front_counter;
	dd.n_part_alive = front_counter;

	cudaStreamSynchronize(main_stream);
	cudaStreamSynchronize(work_stream);

	thrust::gather(thrust::cuda::par.on(work_stream),
			dd.gather_indices.begin(),
			dd.gather_indices.end(),
			thrust::make_zip_iterator(thrust::make_tuple(
					ps.r.begin(), ps.v.begin(),
					ps.flags.begin(), ps.deathtime.begin(), ps.id.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(
					other_ps.r.begin(), other_ps.v.begin(),
					other_ps.flags.begin(), other_ps.deathtime.begin(), other_ps.id.begin())));

	dd.phase_space_id++;

	return pruned;
}

int main(int argv, char** argc)
{
	if (argv < 4)
	{
		std::cerr << "Please specify time step and final time " << argc[0] << " <CURTIME> <TIMESTEP> <FINALTIME> [<MAXPARTICLES>]" << std::endl;
		return -1;
	}

	HostData hd;
	DeviceData dd;

	hd.t = std::stod(argc[1]);
	float t0 = hd.t;
	hd.dt = std::stod(argc[2]);
	hd.t_f = std::stod(argc[3]);

	size_t max_particle = 0;
	if (argv >= 5) max_particle = static_cast<size_t>(std::stoi(argc[4]));

	if (load_data(hd, "pl.in", "ics.in", max_particle, false)) return -1;

	convert_to_barycentric(hd);

	double e0 = energy_planet(hd.r_planet, hd.v_planet, hd.m_planet);
	std::cout << std::setprecision(17);
	std::cout << "e0 (planets) = " << e0 << std::endl;
	std::cout << "t = " << hd.t << std::endl;
	std::cout << "dt = " << hd.dt << std::endl;
	std::cout << "t_f = " << hd.t_f << std::endl;
	std::cout << "n_particle = " << hd.n_part << std::endl;
	std::cout << "==================================" << std::endl;
	std::cout << "Sending initial conditions to GPU." << std::endl;

	transfer_data(hd, dd);
	std::cout << "       Starting simulation.       " << std::endl << std::endl;

	size_t print_every = 10;
	size_t print_counter = 0;

	size_t prune_every = 10;
	size_t prune_counter = 1;

	cudaStream_t work_stream, htd_stream, dth_stream, gather_stream;
	cudaStreamCreate(&work_stream);
	cudaStreamCreate(&htd_stream);
	cudaStreamCreate(&gather_stream);
	cudaStreamCreate(&dth_stream);

	std::atomic<bool> is_pulling_data(false);

	size_t pull_every = 500;
	size_t pull_counter = 0;

	recover_data(hd, dd, dth_stream);
	cudaStreamSynchronize(dth_stream);

	std::cout << "Saving to disk." << std::endl;
	save_data(hd, dd, "pl.part.out", "ics.part.out", "info.part.out");

	auto starttime = std::chrono::high_resolution_clock::now();


	std::ofstream timelog("time.out");
	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);
	timelog << "start " << std::put_time(&tm, "%c %Z") << std::endl;

	while (hd.t < hd.t_f)
	{
		if (print_counter % print_every == 0)
		{
			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> millis = now - starttime;
			double elapsed = millis.count() / 60000;
			double total = millis.count() / 60000 * (hd.t_f - t0) / (hd.t - t0);
			std::cout << "t=" << hd.t << " (" << (hd.t - t0) / (hd.t_f - t0) * 100 << "% " << elapsed << "m elapsed, " << total << "m total " << total - elapsed << "m remain)" << std::endl;
			std::cout << "Error = " << (energy_planet(hd.r_planet, hd.v_planet, hd.m_planet) - e0) / e0 * 100 << ", " << dd.n_part_alive << " particles remaining" << std::endl;

			timelog << std::setprecision(17) << "timing " << hd.t << " " << millis.count() / 60000 << " " << dd.n_part_alive << std::endl;

			// TODO report on barycenter drift
		}

		if (pull_counter >= pull_every && !is_pulling_data)
		{
			pull_counter = 0;
		}

		sia4_planet(hd.coefdt, hd.r_planet, hd.v_planet, hd.m_planet, hd.r_planet_log);

		// Switching between two planet log buffers allows the overhead of copying planet locations to the GPU be removed.
		Dvf64_3& r_log_buffer = dd.log_buffer_id % 2 ? dd.r_planet_log0 : dd.r_planet_log1;

		thrust::copy(thrust::cuda::par.on(htd_stream), hd.r_planet_log.begin(), hd.r_planet_log.end(), r_log_buffer.begin());
		dd.log_buffer_id++;

		if (prune_counter % prune_every == 0)
		{
			cudaStreamSynchronize(dth_stream);
			prune(work_stream, dth_stream, hd, dd);
		}

		DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;

		size_t n = dd.n_part_alive;

		cudaStreamSynchronize(htd_stream);
		cudaStreamSynchronize(gather_stream);
		cudaStreamSynchronize(work_stream);

		auto iterator = thrust::make_zip_iterator(thrust::make_tuple(ps.r.begin(), ps.v.begin(),
					ps.flags.begin(), ps.deathtime.begin(), ps.id.begin()));

		thrust::for_each(thrust::cuda::par.on(work_stream),
				iterator,
				iterator + n,
				SIA4Kernel(dd.coefdt.data().get(), r_log_buffer.data().get(), dd.m_planet.data().get(), hd.n_planet));
		recover_data(hd, dd, dth_stream);

		print_counter++;
		prune_counter++;
		pull_counter++;
		hd.t += hd.dt * TIMEBLOCK_SIZE;
	}


	cudaStreamSynchronize(work_stream);
	prune(work_stream, gather_stream, hd, dd);
	cudaStreamSynchronize(gather_stream);

	std::cout << "Simulation finished. t = " << hd.t << ". n_particle = " << dd.n_part_alive << std::endl;

	cudaStreamSynchronize(dth_stream);
	recover_data(hd, dd, dth_stream);


	auto now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> millis = now - starttime;
	double elapsed = millis.count() / 60000;
	double total = millis.count() / 60000 * (hd.t_f - t0) / (hd.t - t0);
	std::cout << "t=" << hd.t << " (" << (hd.t - t0) / (hd.t_f - t0) * 100 << "% " << elapsed << "m elapsed, " << total << "m total " << total - elapsed << "m remain)" << std::endl;
	std::cout << "Error = " << (energy_planet(hd.r_planet, hd.v_planet, hd.m_planet) - e0) / e0 * 100 << ", " << dd.n_part_alive << " particles remaining" << std::endl;

	timelog << std::setprecision(17) << "timing " << hd.t << " " << millis.count() / 60000 << " " << dd.n_part_alive << std::endl;

	save_data(hd, dd, "pl.out", "ics.out", "info.out");

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return 0;
}

void sia4_planet(const Hvf64& coefdt, Hvf64_3& r, Hvf64_3& v, Hvf64& m, Hvf64_3& r_log)
{
	size_t n = m.size();

	Hvf64_3 a(n);

	for (size_t step = 0; step < TIMEBLOCK_SIZE; step++)
	{
		double adt = coefdt[0];
		double bdt;

		for (size_t i = 0; i < n; i++)
		{
			r[i].x += adt * v[i].x;
			r[i].y += adt * v[i].y;
			r[i].z += adt * v[i].z;
		}

		for (size_t coef = 1; coef < SIMPL_DEGREE; coef++)
		{
			// Write planet locations to log
			for (size_t planet = 0; planet < n; planet++)
			{
				int log_index = step * (SIMPL_DEGREE - 1) * n + (coef - 1) * n + planet;

				r_log[log_index].x = r[planet].x;
				r_log[log_index].y = r[planet].y;
				r_log[log_index].z = r[planet].z;
			}

			adt = coefdt[coef];
			bdt = coefdt[coef + SIMPL_DEGREE];

			newton3d_planet(a, r, m);

			for (size_t i = 0; i < n; i++)
			{
				v[i].x += bdt * a[i].x;
				v[i].y += bdt * a[i].y;
				v[i].z += bdt * a[i].z;
			}
			for (size_t i = 0; i < n; i++)
			{
				r[i].x += adt * v[i].x;
				r[i].y += adt * v[i].y;
				r[i].z += adt * v[i].z;
			}
		}
	}
}

void newton3d_planet(Hvf64_3& a, const Hvf64_3& r, const Hvf64& m)
{
	size_t n = m.size();

	for (size_t i = 0; i < n; i++)
	{
		a[i].x = a[i].y = a[i].z = 0;
	}
	
	for (size_t i = 0; i < n - 1; i++)
	{
		for (size_t j = i + 1; j < n; j++)
		{
			double dx = r[i].x - r[j].x;
			double dy = r[i].y - r[j].y;
			double dz = r[i].z - r[j].z;
			double dist = dx * dx + dy * dy + dz * dz;

			/* rij holds the cube of the distance from i to j */
			double rij = dist * std::sqrt(dist);

			a[i].x -= m[j] * dx / rij;
			a[j].x += m[i] * dx / rij;

			a[i].y -= m[j] * dy / rij;
			a[j].y += m[i] * dy / rij;

			a[i].z -= m[j] * dz / rij;
			a[j].z += m[i] * dz / rij;
		}
	}
}

double energy_planet(const Hvf64_3& r, const Hvf64_3& v, const Hvf64& m)
{
	size_t n = m.size();

	double ke = 0.0;
	double pe = 0.0;

	for (size_t i = 0; i < n; i++)
	{
		ke += 0.5 * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z) * m[i];
	}

	for (size_t i = 0; i < n - 1; i++)
	{
		for (size_t j = i + 1; j < n; j++)
		{
			double dx = r[i].x - r[j].x;
			double dy = r[i].y - r[j].y;
			double dz = r[i].z - r[j].z;

			pe -= m[i] * m[j] / std::sqrt(dx * dx + dy * dy + dz * dz);
		}
	}

	return ke + pe;
}
