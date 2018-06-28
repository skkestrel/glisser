#pragma once
#include "data.h"

class Integrator
{
public:
	virtual void integrate_planets_timeblock(HostPlanetPhaseSpace& pl, float64_t t) = 0;
	virtual void integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t) = 0;
	virtual void integrate_encounter_particle_catchup(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t particle_deathtime_index, size_t planet_index) = 0;
	virtual void gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length) = 0;
};
