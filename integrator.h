#pragma once
#include "data.h"

class Integrator
{
public:
	virtual void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t index, float64_t dt) = 0;
	virtual void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, size_t timestep_index, float64_t dt) = 0;
	virtual void integrate_encounter_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t n_timesteps, float64_t dt) = 0;
};
