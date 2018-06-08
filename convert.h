#pragma once
#include "types.h"
#include "data.h"

void helio_to_jacobi_r_part(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p);
void helio_to_jacobi_v_part(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p);
void helio_to_jacobi_r(HostPlanetPhaseSpace& p);
void helio_to_jacobi_v(HostPlanetPhaseSpace& p);
void to_barycentric(HostData& hd);
void to_heliocentric(HostData& hd);
