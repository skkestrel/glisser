#pragma once
#include "types.h"
#include "data.h"

void jacobi_to_helio_planets(HostPlanetPhaseSpace& pl);
void helio_to_jacobi_r_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p);
void helio_to_jacobi_v_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p);
void helio_to_jacobi_r_planets(HostPlanetPhaseSpace& p);
void helio_to_jacobi_v_planets(HostPlanetPhaseSpace& p);
void to_bary(HostData& hd);
void to_helio(HostData& hd);
