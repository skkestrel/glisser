#pragma once
#include "data.h"

void step_planets(HostPlanetPhaseSpace& pl, size_t index, float64_t dt);
void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t index, float64_t dt);
void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa);
