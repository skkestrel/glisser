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

	Vu32 id;

	inline HostPlanetPhaseSpace() { }
	inline HostPlanetPhaseSpace(size_t n, size_t tb_size, size_t ce_factor):
		n(n), n_alive(n), m(n), eta(n), r(n), v(n), rj(n), vj(n),
		a(n),
		r_log((n - 1) * tb_size * ce_factor), h0_log(tb_size * ce_factor),
		r_log_slow((n - 1) * tb_size), h0_log_slow(tb_size),
       		id(n) { }
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
	size_t tbsize, ce_factor, print_every, dump_every, track_every, max_particle;
	bool resolve_encounters, readmomenta, enable_ascii_track, enable_binary_track, readhybrid, writehybrid, dumpbinary, writehybridbinary, readhybridbinary;

	std::string icsin, plin, hybridin;
	std::string outfolder;

	Configuration();
};

template<typename T>
void write_binary(std::ostream& o, const T& t);
template<typename T>
void read_binary(std::istream& i, T& t);

bool load_data(HostData& hd, const Configuration& config);
void save_data(const HostData& hd, const Configuration& config, bool dump=false, size_t dumpnum=0);
bool read_configuration(std::istream& in, Configuration* out);
void write_configuration(std::ostream& in, const Configuration& config);
std::string joinpath(const std::string& base, const std::string& file);
