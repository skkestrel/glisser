#pragma once
#include "types.cuh"
#include "data.cuh"
#include "integrator.cuh"
#include "wh.h"

class WHCudaIntegrator : public CudaIntegrator
{
public:
	WHIntegrator base;
	Dvf64_3 device_particle_a;

	void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t index, float64_t dt) override;
	void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, size_t timestep_index, float64_t dt) override;
	void integrate_encounter_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t n_timesteps, float64_t dt) override;

	void step_particles_cuda(cudaStream_t& stream, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa, size_t tbsize, float64_t dt) override;
	void upload_data(cudaStream_t& stream) override;

	WHCudaIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa);
};
