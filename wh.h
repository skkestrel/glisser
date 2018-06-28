#pragma once
#include "data.h"
#include "integrator.h"

class WHIntegrator : public Integrator
{
public:
	Vf64 inverse_helio_cubed, inverse_jacobi_cubed;
	Vf64 dist, energy, vdotr;
	Vu8 mask;
	Vf64 mu, eta;

	Vf64_3 planet_rj, planet_vj;
	Vf64_3 planet_a, particle_a;

	Vf64_3 planet_h0_log, planet_h0_log_old;
	Vf64_3 planet_h0_log_slow, planet_h0_log_slow_old;

	Vf64 planet_rh;

	size_t tbsize;
	size_t encounter_n1, encounter_n2;
	double encounter_r1, encounter_r2;

	bool resolve_encounters;

	double dt;

	WHIntegrator();
	WHIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config);

	void integrate_planets_timeblock(HostPlanetPhaseSpace& pl, float64_t t) override;
	void integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t) override;
	void gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length) override;

	void integrate_encounter_particle_catchup(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t particle_deathtime_index, size_t planet_index) override;

	template<bool old>
	void integrate_encounter_particle_step(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t planet_index, size_t timestep_index, double t);

	void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t timestep_index);
	void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, size_t timestep_index);

	void drift_single(float64_t t, float64_t mu, f64_3* r, f64_3* v) const;
	void drift(float64_t t, const Vu8& mask, const Vf64& mu, Vf64_3& r, Vf64_3& v, size_t start, size_t n);

	template<bool old>
	void nonhelio_acc_encounter_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, size_t particle_index, float64_t time, size_t timestep_index, size_t central_planet_index);

	template<bool encounter, bool old>
	void helio_acc_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, float64_t time, size_t timestep_index);

	void helio_acc_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, size_t begin, size_t length, float64_t time, size_t timestep_index);
	void helio_acc_planets(HostPlanetPhaseSpace& p, size_t index);

	template<bool slow, bool old>
	inline f64_3& h0_log_at(size_t timestep);
	template<bool slow, bool old>
	inline const f64_3& h0_log_at(size_t timestep) const;
};

template<>
inline const f64_3& WHIntegrator::h0_log_at<false, false>(size_t timestep) const
{
	return planet_h0_log[timestep];
}
template<>
inline const f64_3& WHIntegrator::h0_log_at<true, false>(size_t timestep) const
{
	return planet_h0_log_slow[timestep];
}
template<>
inline const f64_3& WHIntegrator::h0_log_at<false, true>(size_t timestep) const
{
	return planet_h0_log_old[timestep];
}
template<>
inline const f64_3& WHIntegrator::h0_log_at<true, true>(size_t timestep) const
{
	return planet_h0_log_slow_old[timestep];
}
template<>
inline f64_3& WHIntegrator::h0_log_at<false, false>(size_t timestep)
{
	return planet_h0_log[timestep];
}
template<>
inline f64_3& WHIntegrator::h0_log_at<true, false>(size_t timestep)
{
	return planet_h0_log_slow[timestep];
}
template<>
inline f64_3& WHIntegrator::h0_log_at<false, true>(size_t timestep)
{
	return planet_h0_log_old[timestep];
}
template<>
inline f64_3& WHIntegrator::h0_log_at<true, true>(size_t timestep)
{
	return planet_h0_log_slow_old[timestep];
}

void calculate_planet_metrics(const HostPlanetPhaseSpace& p, double* energy, f64_3* l);
