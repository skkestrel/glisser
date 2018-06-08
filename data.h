#pragma once
#include "types.h"

struct DevicePhaseSpace
{
};

struct DeviceData
{
	// double buffer for HtD transfer of planet locations
	Dvf64_3 r_planet_log0, r_planet_log1;

	Dvf64 m_planet;

	DevicePhaseSpace ps0, ps1;

	Dvu32 gather_indices;

	size_t log_buffer_id, phase_space_id;
	size_t n_part_alive;

	Dvf64 coefdt;

	inline DeviceData() { }
};

struct HostParticlePhaseSpace
{
	size_t n;
	Hvf64_3 r, v, rj, vj;
	Hvf64_3 a;

	Hvu8 flags;
	Hvu32 deathtime;
	Hvu32 id;

	inline HostParticlePhaseSpace() { }
	inline HostParticlePhaseSpace(size_t n) : n(n), r(n), v(n), a(n), flags(n), deathtime(n), id(n) { }
};

struct HostPlanetPhaseSpace
{
	size_t n;
	Hvf64 m, eta;
	Hvf64_3 r, v, rj, vj;
	Hvf64_3 a;

	Hvf64_3 r_log;
	Hvf64_3 h0_log;

	f64_3 bary_r, bary_v;
	f64_3 h0;

	Hvu8 flags;
	Hvu32 deathtime;
	Hvu32 id;

	inline HostPlanetPhaseSpace() { }
	inline HostPlanetPhaseSpace(size_t n, size_t tb_size) :
		n(n), m(n), eta(n), r(n), v(n), rj(n), vj(n),
		a(n),
		r_log(n * tb_size), h0_log(tb_size),
       		flags(n), deathtime(n), id(n) { }
};

struct HostData
{
	HostParticlePhaseSpace particles;
	HostPlanetPhaseSpace planets;

	Hvf64 coefdt;

	double t, dt, t_f;

	inline HostData() { }
};


bool load_data(HostData& hd, std::string plin, std::string icsin, size_t tbsize, size_t max_particle = 0, bool readmomenta = true);
void save_data(const HostData& hd, const DeviceData& dd, std::string plout, std::string icsout, std::string infoout);
void transfer_data(const HostData& hd, DeviceData& dd);
void recover_data(HostData& hd, const DeviceData& dd, cudaStream_t& stream);
