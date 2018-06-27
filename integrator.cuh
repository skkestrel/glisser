#pragma once
#include "integrator.h"

class CudaIntegrator : public Integrator
{
public:
	virtual void step_particles_cuda(cudaStream_t& stream, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa, size_t tbsize, float64_t dt) = 0;
	virtual void upload_data(cudaStream_t& stream) = 0;
};
