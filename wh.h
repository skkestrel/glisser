#pragma once
#include "data.h"

void step_planets(HostPlanetPhaseSpace& pl, size_t index, float64_t dt);
void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t index, float64_t dt);
void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa);
void drift(float64_t t, Hvf64& mu, Hvf64_3& r, Hvf64_3& v, size_t start, size_t n);
void first_step(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, float64_t dt);
