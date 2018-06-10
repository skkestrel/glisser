#include <iomanip>
#include "executor.h"
#include "convert.h"
#include "wh.h"
#include <fstream>

Executor::Executor(HostData& hd, DeviceData& dd, std::ostream& out) : hd(hd), dd(dd), print_every(10), print_counter(0), tbsize(128), output(out) { }

void Executor::init()
{
	to_helio(hd);

	e_0 = 0; // TODO

	output << std::setprecision(17);
	output << "e_0 (planets) = " << e_0 << std::endl;
	output << "t_0 = " << t << std::endl;
	output << "dt = " << dt << std::endl;
	output << "t_f = " << t_f << std::endl;
	output << "n_particle = " << hd.particles.n << std::endl;
	output << "==================================" << std::endl;
	output << "Running for half a time step.     " << std::endl;

	initialize(hd.planets, hd.particles);
	first_step(hd.planets, hd.particles, dt);

	output << "==================================" << std::endl;
	output << "Sending initial conditions to GPU." << std::endl;

	cudaStreamCreate(&work_stream);
	cudaStreamCreate(&htd_stream);
	cudaStreamCreate(&gather_stream);
	cudaStreamCreate(&dth_stream);

	upload_data();
	download_data();

	cudaStreamSynchronize(dth_stream);

	starttime = std::chrono::high_resolution_clock::now();

	output << "       Starting simulation.       " << std::endl << std::endl;

}

void Executor::upload_data()
{
	dd.particles0 = DeviceParticlePhaseSpace(hd.particles.n);
	dd.particles1 = DeviceParticlePhaseSpace(hd.particles.n);

	dd.planets0 = DevicePlanetPhaseSpace(hd.planets.n, hd.tbsize);
	dd.planets1 = DevicePlanetPhaseSpace(hd.planets.n, hd.tbsize);

	/*
	dd.log_buffer_id = 0;
	dd.phase_space_id = 0;

	dd.planets.r_log0 = Dvf64_3(hd.planets.r_log.size());
	dd.planets.r_log1 = Dvf64_3(hd.planets.r_log.size());
	dd.planets.m = Dvf64(hd.n_planet);

	dd.ps0.r = dd.ps0.v = Dvf64_3(hd.n_part);
	dd.ps0.flags = Dvu8(hd.n_part);
	dd.ps0.deathtime = Dvu32(hd.n_part);
	dd.ps0.id = Dvu32(hd.n_part);

	dd.ps1.r = dd.ps1.v = Dvf64_3(hd.n_part);
	dd.ps1.flags = Dvu8(hd.n_part);
	dd.ps1.deathtime = Dvu32(hd.n_part);
	dd.ps1.id = Dvu32(hd.n_part);

	dd.gather_indices = Dvu32(hd.n_part);

	dd.particles.n_alive = 0;
	for (size_t i = 0; i < hd.n_part; i++)
	{
		if (~hd.particles.flags[i] & 0x0001)
		{
			dd.particles.n_alive++;
		}
	}

	dd.coefdt = Dvf64(hd.coefdt.size());

	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;
	thrust::copy(hd.particles.r.begin(), hd.particles.r.end(), ps.r.begin());
	thrust::copy(hd.particles.v.begin(), hd.particles.v.end(), ps.v.begin());

	thrust::copy(hd.planets.m.begin(), hd.planets.m.end(), dd.planets.m.begin());
	thrust::copy(hd.coefdt.begin(), hd.coefdt.end(), dd.coefdt.begin());

	thrust::copy(hd.id.begin(), hd.id.end(), ps.id.begin());
	thrust::copy(hd.deathtime.begin(), hd.deathtime.end(), ps.deathtime.begin());
	thrust::copy(hd.particles.flags.begin(), hd.particles.flags.end(), ps.flags.begin());
	*/
}

void Executor::download_data()
{
	/*
	const DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;

	size_t n = hd.n_part;


	thrust::copy(thrust::cuda::par.on(stream), ps.r.begin(), ps.r.begin() + n, hd.particles.r.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.v.begin(), ps.v.begin() + n, hd.particles.v.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.flags.begin(), ps.flags.begin() + n, hd.particles.flags.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.deathtime.begin(), ps.deathtime.begin() + n, hd.particles.deathtime.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.id.begin(), ps.id.begin() + n, hd.particles.id.begin());
	*/

	/*
	// zip will crash the program

	auto iterator = thrust::make_zip_iterator(thrust::make_tuple(
				ps.r.begin(),
				ps.v.begin(),
				ps.flags.begin(), ps.deathtime.begin(), ps.id.begin()));
	thrust::copy(thrust::cuda::par.on(stream),
			iterator,
			iterator + n,
			thrust::make_zip_iterator(thrust::make_tuple(
					hd.particles.r.begin(), hd.particles.v.begin(),
					hd.particles.flags.begin(), hd.deathtime.begin(), hd.id.begin())));
	 */
}


double Executor::time() const
{
	auto now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> millis = now - starttime;
	return millis.count() / 60000;
}

void Executor::loop()
{
	if (print_counter % print_every == 0)
	{
		double e = 0;
		double elapsed = time();
		double total = elapsed * (t_f - t_0) / (t - t_0);

		output << "t=" << t << " (" << elapsed / total * 100 << "% " << elapsed << "m elapsed, " << total << "m total " << total - elapsed << "m remain)" << std::endl;
		output << "Error = " << (e - e_0) / e_0 * 100 << ", " << hd.particles.n_alive << " particles remaining" << std::endl;
	}
	print_counter++;

	// advance planets
	for (size_t i = 0; i < tbsize; i++)
	{
		step_planets(hd.planets, t, i, dt);
	}

	/*
	// copy planet log to the buffer currently NOT in use
	Dvf64_3& r_log_buffer = dd.log_buffer_id % 2 ? dd.r_planet_log0 : dd.r_planet_log1;
	thrust::copy(thrust::cuda::par.on(htd_stream), hd.r_planet_log.begin(), hd.r_planet_log.end(), r_log_buffer.begin());
	dd.log_buffer_id++;
	*/

	for (size_t i = 0; i < tbsize; i++)
	{
		step_particles(hd.planets, hd.particles, t, i, dt);
	}

	resync();


	/*
	// pruning might have moved things to the other buffer
	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;


	size_t n = hd.particles.n_alive;
	cudaStreamSynchronize(htd_stream);
	cudaStreamSynchronize(gather_stream);
	cudaStreamSynchronize(work_stream);
	*/


	/*
	std::ofstream anilog("animation.out", std::ios_base::app);

	to_bary(hd);
	for (size_t j = 0; j < hd.planets.n; j++)
	{
		anilog << hd.planets.r[j].x << " " << hd.planets.r[j].y << " " << hd.planets.r[j].z << std::endl << std::flush;
	}
	for (size_t j = 0; j < hd.particles.n; j++)
	{
		anilog << hd.particles.r[j].x << " " << hd.particles.r[j].y << " " << hd.particles.r[j].z << std::endl << std::flush;
	}
	to_helio(hd);
	*/


	t += dt * tbsize;
}

void Executor::resync()
{
	hd.particles.n_alive = hd.particles.n;

	for (size_t i = 0; i < hd.particles.n; i++)
	{
		if (hd.particles.flags[i] > 0)
		{
			hd.particles.n_alive--;
		}
	}
	
	/*
	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;
	DevicePhaseSpace& other_ps = dd.phase_space_id % 2 ? dd.ps1 : dd.ps0;

	thrust::copy(thrust::cuda::par.on(work_stream), ps.flags.begin(), ps.flags.begin() + hd.particles.n_alive, hd.part_flags.begin());
	Hvu32 indices = Hvu32(hd.n_part);

	// start placing alive particles at the front and dead particles at the back
	uint32_t front_counter = 0;
	uint32_t back_counter = hd.particles.n_alive - 1;

	for (size_t i = 0; i < hd.particles.n_alive; i++)
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
	for (size_t i = hd.particles.n_alive; i < hd.n_part; i++)
	{
		indices[i] = i;
	}

	cudaStreamSynchronize(work_stream);
	thrust::copy(thrust::cuda::par.on(work_stream), indices.begin(), indices.end(), dd.gather_indices.begin());

	size_t resyncd = hd.particles.n_alive - front_counter;
	hd.particles.n_alive = front_counter;

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
*/
}


void Executor::finish()
{
	cudaStreamSynchronize(work_stream);
	resync();

	cudaStreamSynchronize(gather_stream);

	output << "Simulation finished. t = " << t << ". n_particle = " << hd.particles.n_alive << std::endl;

	cudaStreamSynchronize(dth_stream);
	download_data();
}
