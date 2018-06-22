#pragma once
#include "data.h"

struct WHIntermediate
{
	Vf64 inverse_helio_cubed, inverse_jacobi_cubed;
	Vf64 dist, energy, vdotr;
	Vu8 mask;
	Vf64 mu;

	Vf64_3 r, v;

	WHIntermediate() { }
	WHIntermediate(size_t n) : inverse_helio_cubed(n), inverse_jacobi_cubed(n),
       		dist(n), energy(n), vdotr(n),
		mask(n), mu(n),
		r(n), v(n) { }
};

void step_planets(HostPlanetPhaseSpace& pl, WHIntermediate& alloc, float64_t t, size_t index, float64_t dt);
void step_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, float64_t t, size_t timestep_index, float64_t dt);
void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, WHIntermediate& alloc, size_t begin, size_t length, float64_t t, size_t timestep_index, float64_t dt);
void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, WHIntermediate& alloc);
void calculate_planet_metrics(const HostPlanetPhaseSpace& p, WHIntermediate& alloc, double* energy, f64_3* l);
