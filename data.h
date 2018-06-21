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
	Vf32 deathtime; Vu32 id; 
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
	Vf64_3 r_log_slow;
	Vf64_3 h0_log_slow;

	f64_3 bary_r, bary_v;

	Vu16 deathflags;
	Vf32 deathtime;
	Vu32 id;

	inline HostPlanetPhaseSpace() { }
	inline HostPlanetPhaseSpace(size_t n, size_t tb_size, size_t ce_factor):
		n(n), n_alive(n), m(n), eta(n), r(n), v(n), rj(n), vj(n),
		a(n),
		r_log((n - 1) * tb_size * ce_factor), h0_log(tb_size * ce_factor),
		r_log_slow((n - 1) * tb_size), h0_log_slow(tb_size),
       		deathflags(n), deathtime(n), id(n) { }
};

struct HostPlanetSnapshot
{
	size_t n, n_alive;
	Vf64 m;
	Vf64_3 r, v;
	Vu32 id;

	inline HostPlanetSnapshot() { }
	inline HostPlanetSnapshot(const HostPlanetPhaseSpace& o) :
		n(o.n), n_alive(o.n_alive), m(o.m), r(o.r), v(o.v), id(o.id) { }
};

struct HostData
{
	HostParticlePhaseSpace particles;
	HostPlanetPhaseSpace planets;
	HostPlanetSnapshot planets_snapshot;

	inline HostData() { }
};

struct Configuration
{
	double t, t_0, t_f, dt;
	size_t tbsize, ce_factor, print_every, dump_every, periodic_every, max_particle;
	bool resolve_encounters, readmomenta, enable_ascii_track, enable_binary_track;
	std::string icsin, plin;
	std::string outfolder;

	Configuration();
};

bool load_data(HostData& hd, std::string plin, std::string icsin, size_t tbsize, size_t ce_factor, size_t max_particle = 0, bool readmomenta = true);
void save_data(const HostData& hd, std::string plout, std::string icsout);
bool read_configuration(std::istream& in, Configuration* out);
void write_configuration(std::ostream& in, const Configuration& config);
