#pragma once
#include "data.h"

void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t index, float64_t dt);
void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, float64_t t, size_t index, float64_t dt);
void step_particles_cuda(cudaStream_t& stream, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa, size_t tbsize, float64_t dt);
void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa);
void drift(float64_t t, Hvu8& mask, Hvf64& mu, Hvf64_3& r, Hvf64_3& v, size_t start, size_t n);

void calculate_planet_metrics(const HostPlanetPhaseSpace& p, double* energy, f64_3* l);
