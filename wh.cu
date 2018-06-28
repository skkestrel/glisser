#include "types.cuh"
#include "wh.cuh"
#include "convert.h"
#include "util.cuh"

const size_t MAXKEP = 5;
const float64_t TOLKEP = 1E-14;

struct MVSKernel
{
	const float64_t* planet_m;
	const float64_t mu;
	const f64_3* planet_h0_log;
	const f64_3* planet_r_log;
	const size_t planet_n;
	const size_t tbsize;
	const float64_t dt;

	const float64_t r2;
	const float64_t* planet_rh;

	MVSKernel(const DevicePlanetPhaseSpace& planets, const Dvf64_3& h0_log, const Dvf64& planet_rh, double r2, size_t tbsize, float64_t dt)
		: planet_m(planets.m.data().get()), mu(planets.m[0]), planet_h0_log(h0_log.data().get()),
		planet_r_log(planets.r_log.data().get()), planet_n(planets.n_alive), tbsize(tbsize), dt(dt),
		r2(r2), planet_rh(planet_rh.data().get())
	{ }

	__host__ __device__
	void kepeq(double dM, double ecosEo, double esinEo, double* dE, double* sindE, double* cosdE, uint16_t& flags) const
	{
		double f, fp, delta;

		*sindE = sin(*dE);
		*cosdE = cos(*dE);

		for (size_t i = 0; i < MAXKEP; i++)
		{
			f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
			fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
			delta = -f / fp;

			*dE += delta;
			*sindE = sin(*dE);
			*cosdE = cos(*dE);
		}

		flags = static_cast<uint16_t>(flags | ((fabs(delta) > TOLKEP) << 3));
	}

	__host__ __device__
	void drift(f64_3& r, f64_3& v, uint16_t& flags) const
	{
		float64_t dist = sqrt(r.lensq());
		float64_t vdotr = v.x * r.x + v.y * r.y + v.z * r.z;

		float64_t energy = v.lensq() * 0.5 - mu / dist;

		flags = static_cast<uint16_t>(flags | ((energy >= 0) << 2));

		float64_t a = -0.5 * mu / energy;
		float64_t n_ = sqrt(mu / (a * a * a));
		float64_t ecosEo = 1.0 - dist / a;
		float64_t esinEo = vdotr / (n_ * a * a);
		float64_t e = sqrt(ecosEo * ecosEo + esinEo * esinEo);

		// subtract off an integer multiple of complete orbits
		float64_t dM = this->dt * n_ - M_2PI * (int) (dt * n_ / M_2PI);

		// remaining time to advance
		float64_t dt = dM / n_;

		// call kepler equation solver with initial guess in dE already
		float64_t dE = dM - esinEo + esinEo * cos(dM) + ecosEo * sin(dM);
		float64_t sindE, cosdE;
		kepeq(dM, ecosEo, esinEo, &dE, &sindE, &cosdE, flags);

		float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
		float64_t f = 1.0 + a * (cosdE - 1.0) / dist;
		float64_t g = dt + (sindE - dE) / n_;
		float64_t fdot = -n_ * sindE * a / (dist * fp);
		float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

		f64_3 r0 = r;
		r = r0 * f + v * g;
		v = r0 * fdot + v * gdot;
	}

	template<typename Tuple>
	__host__ __device__
	void operator()(Tuple args) const
	{
		f64_3 r = thrust::get<0>(thrust::get<0>(args));
		f64_3 v = thrust::get<1>(thrust::get<0>(args));
		uint16_t flags = thrust::get<2>(thrust::get<0>(args));

		uint32_t deathtime_index = 0;
		f64_3 a = thrust::get<1>(args);

		size_t tbsize = this->tbsize;
		const f64_3* h0_log = this->planet_h0_log;
		const f64_3* r_log = this->planet_r_log;
		const float64_t* m = this->planet_m;
		const float64_t* rh = this->planet_rh;
		float64_t r2 = this->r2;
		float64_t dt = this->dt;

		for (uint32_t step = 0; step < static_cast<uint32_t>(tbsize); step++)
		{
			if (flags == 0)
			{
				// kick
				v = v + a * (dt / 2);

				drift(r, v, flags);

				a = h0_log[step];

				// planet 0 is not counted
				for (uint32_t i = 1; i < static_cast<uint32_t>(planet_n); i++)
				{
					f64_3 dr = r - r_log[step * (planet_n - 1) + i - 1];

					float64_t rad = dr.lensq();

					if (rad < rh[i] * rh[i] * r2 * r2 && flags == 0)
					{
						flags = static_cast<uint16_t>(flags | (i << 8) | 0x0001);
					}

					float64_t inv3 = 1. / (rad * sqrt(rad));
					float64_t fac = m[i] * inv3;

					a -= dr * fac;
				}

				float64_t rad = r.lensq();
				if (rad < rh[0] * rh[0] * r2 * r2)
				{
					flags = flags | 0x0001;
				}
				if (rad > 200 * 200)
				{
					flags = flags | 0x0002;
				}


				v = v + a * (dt / 2);

				deathtime_index = step + 1;
			}
		}

		thrust::get<0>(thrust::get<0>(args)) = r;
		thrust::get<1>(thrust::get<0>(args)) = v;
		thrust::get<2>(thrust::get<0>(args)) = flags;
		thrust::get<3>(thrust::get<0>(args)) = deathtime_index;

		thrust::get<1>(args) = a;
	}
};

WHCudaIntegrator::WHCudaIntegrator() { }

WHCudaIntegrator::WHCudaIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config)
	: base(pl, pa, config)
{
	device_h0_log_0 = Dvf64_3(config.tbsize);
	device_h0_log_1 = Dvf64_3(config.tbsize);
	device_particle_a = Dvf64_3(pa.n);

	device_planet_rh = Dvf64(pl.n);

	memcpy_htd(device_planet_rh, base.planet_rh, 0);
	cudaStreamSynchronize(0);
}

Dvf64_3& WHCudaIntegrator::device_h0_log(size_t planet_data_id)
{
	return planet_data_id % 2 ? device_h0_log_1 : device_h0_log_0;
}

void WHCudaIntegrator::upload_planet_log_cuda(cudaStream_t stream, size_t planet_data_id)
{
	memcpy_htd(device_h0_log(planet_data_id), base.planet_h0_log.get<true, false>(), stream);
	cudaStreamSynchronize(stream);
}

void WHCudaIntegrator::gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length)
{
	base.gather_particles(indices, begin, length);
}

void WHCudaIntegrator::integrate_planets_timeblock(HostPlanetPhaseSpace& pl, float64_t t)
{
	base.integrate_planets_timeblock(pl, t);
}

void WHCudaIntegrator::integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t)
{
	base.integrate_particles_timeblock(pl, pa, begin, length, t);
}

void WHCudaIntegrator::integrate_encounter_particle_catchup(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t particle_deathtime_index, size_t planet_index)
{
	base.integrate_encounter_particle_catchup(pl, pa, particle_index, particle_deathtime_index, planet_index);
}

void WHCudaIntegrator::upload_data_cuda(cudaStream_t stream, size_t begin, size_t length)
{
	memcpy_htd(device_particle_a, base.particle_a, stream, begin, begin, length);
	cudaStreamSynchronize(stream);
}

void WHCudaIntegrator::integrate_particles_timeblock_cuda(cudaStream_t stream, size_t planet_data_id, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa)
{
	if (pa.n_alive > 0)
	{
		auto it = thrust::make_zip_iterator(thrust::make_tuple(pa.begin(), device_begin()));
		thrust::for_each(thrust::cuda::par.on(stream), it, it + pa.n_alive, MVSKernel(pl, device_h0_log(planet_data_id), device_planet_rh, base.encounter_r2, base.tbsize, base.dt));
	}
}
