#pragma once
#include "types.h"
#include <string>
#include <utility>
#include <cstdlib>
using size_t = std::size_t;

struct HostParticlePhaseSpace
{
	size_t n, n_alive, n_encounter;

	Vf64_3 r, v;
	Vf64_3 a;

	Vu16 deathflags;
	Vf32 deathtime;
	Vu32 id;

	inline HostParticlePhaseSpace() { }
	inline HostParticlePhaseSpace(size_t n) : n(n), n_alive(n), n_encounter(0), r(n), v(n), a(n), deathflags(n), deathtime(n), id(n) { }

	void stable_partition_alive();
};

struct HostPlanetPhaseSpace
{
	size_t n, n_alive;
	Vf64 m, eta;
	Vf64_3 r, v, rj, vj;
	Vf64_3 a;

	Vf64_3 r_log;
	Vf64_3 h0_log;

	Vf64_3 r_log_prev;
	Vf64_3 h0_log_prev;

	f64_3 bary_r, bary_v;

	Vu16 deathflags;
	Vf32 deathtime;
	Vu32 id;

	inline HostPlanetPhaseSpace() { }
	inline HostPlanetPhaseSpace(size_t n, size_t tb_size) :
		n(n), n_alive(n), m(n), eta(n), r(n), v(n), rj(n), vj(n),
		a(n),
		r_log((n - 1) * tb_size), h0_log(tb_size),
		r_log_prev((n - 1) * tb_size), h0_log_prev(tb_size),
       		deathflags(n), deathtime(n), id(n) { }
};

struct HostData
{
	HostParticlePhaseSpace particles, encounter_particles;
	HostPlanetPhaseSpace planets;

	size_t tbsize;

	inline HostData() { }
};


bool load_data(HostData& hd, std::string plin, std::string icsin, size_t tbsize, size_t max_particle = 0, bool readmomenta = true);
void save_data(const HostData& hd, std::string plout, std::string icsout);
