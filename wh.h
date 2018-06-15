#pragma once
#include "data.h"

void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t index, float64_t dt);
void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, float64_t t, size_t index, float64_t dt);
void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa);
void drift(float64_t t, Vu8& mask, Vf64& mu, Vf64_3& r, Vf64_3& v, size_t start, size_t n);

void calculate_planet_metrics(const HostPlanetPhaseSpace& p, double* energy, f64_3* l);
