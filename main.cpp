/*************************************************************

.---. .            .         .   .-,--.                        
\___  |- ,-. ,-. ,-| . . ,-. |-   `|__/ ,-. .  , ,-. ,-. . ,-. 
    \ |  ,-| |   | | | | `-. |    )  \  |-' | /  |-' |   | |-' 
`---' `' `-^ '   `-^ `-^ `-' `'   `-  ` `-' `'   `-' '   ' `-' 

*************************************************************/

#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <ctime>
#include <thread>
#include <atomic>
#include <iomanip>

#include "types.h"
#include "data.h"
#include "wh.h"
#include "convert.h"

#ifdef GCC
#define cudaStreamCreate(...)
#define cudaStreamSynchronize(...)
#else
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>
#endif

const int SIMPL_DEGREE = 4;
const int TIMEBLOCK_SIZE = 1;

size_t prune(cudaStream_t& main_stream, cudaStream_t& work_stream, HostData& hd, DeviceData& dd)
{
	(void) main_stream;
	(void) work_stream;
	
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

	size_t pruned = hd.particles.n_alive - front_counter;
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

	return pruned;
*/
	return 0;
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

	if (load_data(hd, "pl.in", "ics.in", TIMEBLOCK_SIZE, max_particle, false)) return -1;

	double e0 = 0;
	std::cout << std::setprecision(17);
	std::cout << "e0 (planets) = " << e0 << std::endl;
	std::cout << "t = " << hd.t << std::endl;
	std::cout << "dt = " << hd.dt << std::endl;
	std::cout << "t_f = " << hd.t_f << std::endl;
	std::cout << "n_particle = " << hd.particles.n << std::endl;
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

	to_helio(hd);
	initialize(hd.planets, hd.particles);
	first_step(hd.planets, hd.particles, hd.dt);
	
	while (hd.t < hd.t_f)
	{
		if (print_counter % print_every == 0)
		{
			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> millis = now - starttime;
			double elapsed = millis.count() / 60000;
			double total = millis.count() / 60000 * (hd.t_f - t0) / (hd.t - t0);
			std::cout << "t=" << hd.t << " (" << (hd.t - t0) / (hd.t_f - t0) * 100 << "% " << elapsed << "m elapsed, " << total << "m total " << total - elapsed << "m remain)" << std::endl;
			std::cout << "Error = " << (0 - e0) / e0 * 100 << ", " << hd.particles.n_alive << " particles remaining" << std::endl;

			timelog << std::setprecision(17) << "timing " << hd.t << " " << millis.count() / 60000 << " " << hd.particles.n_alive << std::endl;
		}

		if (pull_counter >= pull_every && !is_pulling_data)
		{
			pull_counter = 0;
		}

		// advance planets
		for (size_t i = 0; i < TIMEBLOCK_SIZE; i++)
		{
			step_planets(hd.planets, hd.t, i, hd.dt);
		}

		/*
		// copy planet log to the buffer currently NOT in use
		Dvf64_3& r_log_buffer = dd.log_buffer_id % 2 ? dd.r_planet_log0 : dd.r_planet_log1;
		thrust::copy(thrust::cuda::par.on(htd_stream), hd.r_planet_log.begin(), hd.r_planet_log.end(), r_log_buffer.begin());
		dd.log_buffer_id++;
		*/

		if (prune_counter % prune_every == 0)
		{
/*
			cudaStreamSynchronize(dth_stream);
*/
			prune(work_stream, dth_stream, hd, dd);
		}

/*
		// pruning might have moved things to the other buffer
		DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;


		size_t n = hd.particles.n_alive;
		cudaStreamSynchronize(htd_stream);
		cudaStreamSynchronize(gather_stream);
		cudaStreamSynchronize(work_stream);

		// recover data from gpu (optional)
		// recover_data(hd, dd, dth_stream);
*/

		for (size_t i = 0; i < TIMEBLOCK_SIZE; i++)
		{
			step_particles(hd.planets, hd.particles, hd.t, i, hd.dt);
		}

		print_counter++;
		prune_counter++;
		pull_counter++;
		hd.t += hd.dt * TIMEBLOCK_SIZE;
	}

	cudaStreamSynchronize(work_stream);
	prune(work_stream, gather_stream, hd, dd);
	cudaStreamSynchronize(gather_stream);

	std::cout << "Simulation finished. t = " << hd.t << ". n_particle = " << hd.particles.n_alive << std::endl;

	cudaStreamSynchronize(dth_stream);
	recover_data(hd, dd, dth_stream);


	auto now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> millis = now - starttime;
	double elapsed = millis.count() / 60000;
	double total = millis.count() / 60000 * (hd.t_f - t0) / (hd.t - t0);
	std::cout << "t=" << hd.t << " (" << (hd.t - t0) / (hd.t_f - t0) * 100 << "% " << elapsed << "m elapsed, " << total << "m total " << total - elapsed << "m remain)" << std::endl;
	std::cout << "Error = " << (0 - e0) / e0 * 100 << ", " << hd.particles.n_alive << " particles remaining" << std::endl;

	timelog << std::setprecision(17) << "timing " << hd.t << " " << millis.count() / 60000 << " " << hd.particles.n_alive << std::endl;

	save_data(hd, dd, "pl.out", "ics.out", "info.out");

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return 0;
}
