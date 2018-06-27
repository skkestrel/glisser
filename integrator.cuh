#pragma once
#include "integrator.h"

class CudaIntegrator : public Integrator
{
public:
	virtual void integrate_particles_timeblock_cuda(cudaStream_t& stream, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa) = 0;
	virtual void upload_data_cuda(cudaStream_t& stream) = 0;
};
