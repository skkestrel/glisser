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
	Dvf64 device_planet_rh;

	using device_iterator = decltype(device_particle_a.begin());

	void integrate_planets_timeblock(HostPlanetPhaseSpace& pl, float64_t t) override;
	void integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t) override;
	void integrate_encounter_particle_catchup(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t particle_deathtime_index) override;
	void gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length) override;

	void integrate_particles_timeblock_cuda(cudaStream_t& stream, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa) override;
	void upload_data_cuda(cudaStream_t& stream) override;

	WHCudaIntegrator();
	WHCudaIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config);

	inline device_iterator device_begin()
	{
		return device_particle_a.begin();
	}
};
