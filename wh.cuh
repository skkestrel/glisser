#pragma once
#include "types.cuh"
#include "data.cuh"
#include "wh.h"

void step_particles_cuda(cudaStream_t& stream, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa, size_t tbsize, float64_t dt);
