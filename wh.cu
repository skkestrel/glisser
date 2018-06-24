#include "types.cuh"
#include "wh.cuh"
#include "wh.h"
#include "convert.h"

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

	MVSKernel(const DevicePlanetPhaseSpace& planets, size_t tbsize, float64_t dt)
		: planet_m(planets.m.data().get()), mu(planets.m[0]), planet_h0_log(planets.h0_log.data().get()),
		planet_r_log(planets.r_log.data().get()), planet_n(planets.n_alive), tbsize(tbsize), dt(dt)
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
		f64_3 r = thrust::get<0>(args);
		f64_3 v = thrust::get<1>(args);
		f64_3 a = thrust::get<2>(args);
		uint32_t deathtime_index = 0;
		uint16_t flags = thrust::get<3>(args);

		size_t tbsize = this->tbsize;
		const f64_3* h0_log = this->planet_h0_log;
		const f64_3* r_log = this->planet_r_log;
		const float64_t* m = this->planet_m;
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

					if (rad < 0.5 * 0.5)
					{
						flags = static_cast<uint16_t>(flags | (i << 8) | 0x0001);
					}

					float64_t inv3 = 1. / (rad * sqrt(rad));
					float64_t fac = m[i] * inv3;

					a -= dr * fac;
				}

				float64_t rad = r.lensq();
				if (rad < 0.5 * 0.5)
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

		thrust::get<0>(args) = r;
		thrust::get<1>(args) = v;
		thrust::get<2>(args) = a;
		thrust::get<3>(args) = flags;
		thrust::get<4>(args) = deathtime_index;
	}
};

void step_particles_cuda(cudaStream_t& stream, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa, size_t tbsize, float64_t dt)
{
	if (pa.n_alive > 0)
	{
		thrust::for_each(thrust::cuda::par.on(stream), pa.mvs_begin(), pa.mvs_begin() + pa.n_alive, MVSKernel(pl, tbsize, dt));
	}
}
