#include "types.h"

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

	inline DeviceData() { }
};

struct HostData
{
	HostParticlePhaseSpace particles;
	HostPlanetPhaseSpace planets;

	Hvf64 coefdt;

	Hvu8 part_flags;
	Hvu32 deathtime;
	Hvu32 id;

	double t, dt, t_f;

	inline HostData() { }
};

struct HostParticlePhaseSpace
{
	size_t n;
	Hvf64_3 r, v, rj, vj;
	Hvf64_3 a;
}

struct HostPlanetPhaseSpace
{
	size_t n;
	Hvf64 m, eta;
	Hvf64_3 r, v, rj, vj;
	Hvf64_3 a;
	f64_3 bary_r, bary_v;
}

bool load_data(HostData& hd, std::string plin, std::string icsin, size_t max_particle = 0, bool readmomenta = true);
void save_data(const HostData& hd, const DeviceData& dd, std::string plout, std::string icsout, std::string infoout);
void transfer_data(const HostData& hd, DeviceData& dd);
void recover_data(HostData& hd, const DeviceData& dd, cudaStream_t& stream);
