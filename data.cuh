#pragma once
#include "types.cuh"

struct DeviceParticlePhaseSpace
{
	Dvf64_3 r, v, a;

	Dvu32 flags;
	Dvu32 deathtime_index;
	Dvu32 id;

	size_t n_total;
	size_t n_alive;
	size_t n_encounter;

	inline DeviceParticlePhaseSpace() { }
	inline DeviceParticlePhaseSpace(size_t n) : r(n), v(n), a(n), flags(n), deathtime_index(n), id(n), n_total(n), n_alive(n) { }

	using iterator_tuple = decltype(thrust::make_tuple(r.begin(), v.begin(), a.begin(), flags.begin(), deathtime_index.begin(), id.begin()));
	using iterator = thrust::zip_iterator<iterator_tuple>;

	using mvs_iterator_tuple = decltype(thrust::make_tuple(r.begin(), v.begin(), a.begin(), flags.begin(), deathtime_index.begin()));
	using mvs_iterator = thrust::zip_iterator<mvs_iterator_tuple>;

	inline iterator begin()
	{
		return thrust::make_zip_iterator(thrust::make_tuple(
			r.begin(), v.begin(), a.begin(), flags.begin(), deathtime_index.begin(), id.begin()));
	}

	inline mvs_iterator mvs_begin()
	{
		return thrust::make_zip_iterator(thrust::make_tuple(
			r.begin(), v.begin(), a.begin(), flags.begin(), deathtime_index.begin()));
	}
};

struct DevicePlanetPhaseSpace
{
	Dvf64 m;
	Dvf64_3 r_log;
	Dvf64_3 h0_log;

	size_t n_total;
	size_t n_alive;

	inline DevicePlanetPhaseSpace() { }
	inline DevicePlanetPhaseSpace(size_t n, size_t tbsize) : m(n), r_log(n * tbsize), h0_log(tbsize),
		n_total(n), n_alive(n) { }
};

struct DeviceData
{
	// double buffer for HtD transfer of planet locations
	DeviceParticlePhaseSpace particles0, particles1;
	DevicePlanetPhaseSpace planets0, planets1;

	size_t particle_data_id, planet_data_id;

	inline DeviceData() { }
	inline DeviceParticlePhaseSpace& particle_phase_space() { return particle_data_id % 2 ? particles1 : particles0; }
	inline DevicePlanetPhaseSpace& planet_phase_space() { return planet_data_id % 2 ? planets1 : planets0; }
};

