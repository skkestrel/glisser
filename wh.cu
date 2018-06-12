#include "wh.h"
#include "convert.h"

#include <iostream>
#include <cmath>

const size_t MAXKEP = 10;
const float64_t TOLKEP = 1E-13;

struct MVSKernel
{
	const float64_t* planet_m;
	float64_t mu;
	const f64_3* planet_h0_log;
	const f64_3* planet_r_log;
	size_t planet_n;
	size_t tbsize;
	float64_t dt;

	MVSKernel(const DevicePlanetPhaseSpace& planets, size_t tbsize, float64_t dt)
		: planet_m(planets.m.data().get()), mu(planets.m[0]), planet_h0_log(planets.h0_log.data().get()),
		planet_r_log(planets.r_log.data().get()), planet_n(planets.n_alive), tbsize(tbsize), dt(dt)
	{ }

	__host__ __device__
	void kepeq(double dM, double ecosEo, double esinEo, double* dE, double* sindE, double* cosdE, uint8_t& flags) const
	{
		double f, fp, delta;

		*sindE = sin(*dE);
		*cosdE = cos(*dE);

		// TODO maxkep?
		for (size_t i = 0; i < 10; i++)
		{
			f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
			fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
			delta = -f / fp;

			*dE += delta;
			*sindE = sin(*dE);
			*cosdE = cos(*dE);
		}

		flags = flags | ((fabs(delta) < TOLKEP) << 2);
	}

	__host__ __device__
	void drift(f64_3& r, f64_3& v, uint8_t& flags) const
	{
		float64_t dist = sqrt(r.lensq());
		float64_t vdotr = v.x * r.x + v.y * r.y + v.z * r.z;

		float64_t energy = v.lensq() * 0.5 - mu / dist;

		flags = flags | ((energy >= 0) << 1);

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
		uint32_t deathtime_index = static_cast<uint8_t>(-1);
		uint8_t flags = 0;

		size_t tbsize = this->tbsize;
		const f64_3* h0_log = this->planet_h0_log;
		const f64_3* r_log = this->planet_r_log;
		const float64_t* m = this->planet_m;
		float64_t dt = this->dt;

		for (size_t step = 0; step < tbsize; step++)
		{
			f64_3 r_temp = r;
			f64_3 v_temp = v;
			drift(r_temp, v_temp, flags);


			f64_3 a = *h0_log;
			h0_log++;

			// planet 0 is not counted
			r_log++;
			for (size_t i = 1; i < planet_n; i++)
			{
				f64_3 dr = r - *r_log;
				r_log++;

				float64_t rad = dr.lensq();

				flags = flags | (rad < 0.5 * 0.5);

				float64_t inv3 = 1. / (rad * sqrt(rad));
				float64_t fac = m[i] * inv3;

				a -= dr * fac;
			}

			v_temp = v_temp + a * dt;

			if (flags == 0)
			{
				r = r_temp;
				v = v_temp;
				deathtime_index = step;
			}
		}

		thrust::get<0>(args) = r;
		thrust::get<1>(args) = v;
		thrust::get<2>(args) = flags;
		thrust::get<3>(args) = deathtime_index;
	}
};

void helio_acc_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, float64_t time, size_t index)
{
	for (size_t i = 0; i < p.n; i++)
	{
		p.a[i] = pl.h0_log[index];

		for (size_t j = 1; j < pl.n; j++)
		{
			f64_3 dr = p.r[i] - pl.r[j];
			float64_t rji2 = dr.lensq();
			float64_t irij3 = 1. / (rji2 * std::sqrt(rji2));
			float64_t fac = pl.m[j] * irij3;

			p.a[i] -= dr * fac;

			if (rji2 < 0.5 * 0.5)
			{
				p.deathtime[i] = time;
				p.flags[i] = j;
			}
		}
	}
}

void helio_acc_planets(HostPlanetPhaseSpace& p, size_t index)
{
	Hvf64 inverse_helio_cubed(p.n), inverse_jacobi_cubed(p.n);

	for (size_t i = 1; i < p.n; i++)
	{
		float64_t r2 = p.r[i].lensq();
		inverse_helio_cubed[i] = 1. / (std::sqrt(r2) * r2);
		r2 = p.rj[i].lensq();
		inverse_jacobi_cubed[i] = 1. / (std::sqrt(r2) * r2);
        }
	
        // compute common heliocentric acceleration
	f64_3 a_common(0);
	for (size_t i = 2; i < p.n; i++)    
	{
		float64_t mfac = p.m[i] * inverse_helio_cubed[i];
		a_common -= p.r[i] * mfac;
        }

        // Load this into all the arrays
	for (size_t i = 1; i < p.n; i++)    
	{
		p.a[i] = a_common;
        }

	p.h0_log[index] = a_common - p.r[1] * p.m[1] * inverse_helio_cubed[1];
	
	// Now do indirect acceleration ; note that planet 1 does not receive a contribution 
	for (size_t i = 2; i < p.n; i++)    
	{
		p.a[i] += (p.rj[i] * inverse_jacobi_cubed[i] - p.r[i] * inverse_helio_cubed[i]) * p.m[0];
        }
	
	/* next term ; again, first planet does not participate */
	f64_3 a_accum(0);
	for (size_t i = 2; i < p.n; i++)    
	{
		float64_t mfac = p.m[i] * p.m[0] * inverse_jacobi_cubed[i] / p.eta[i-1];
		a_accum += p.r[i] * mfac;
		p.a[i] += a_accum;
        }

	/* Finally, incorporate the direct accelerations */
	for (size_t i = 1; i < p.n - 1; i++)    
	{
		for (size_t j = i + 1; j < p.n; j++)    
		{
			f64_3 dr = p.r[j] - p.r[i];
			float64_t r2 = dr.lensq();
			float64_t irij3 = 1. / (r2 * std::sqrt(r2));

			float64_t mfac = p.m[i] * irij3;
			p.a[j] -= dr * mfac;

			// acc. on i is just negative, with m[j] instead
			mfac = p.m[j] * irij3;
			p.a[i] += dr * mfac;
		}
	}
}

size_t kepeq(double dM, double ecosEo, double esinEo, double* dE, double* sindE, double* cosdE)
{
	double f,fp, delta;

	*sindE = std::sin( *dE);
	*cosdE = std::cos( *dE);

	size_t i;
	for (i = 0; i < MAXKEP; i++)
	{
		f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
		fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
		delta = -f / fp;
		if (std::fabs(delta) < TOLKEP)
		{
			goto done;
		}

		*dE += delta;
		*sindE = std::sin(*dE);
		*cosdE = std::cos(*dE);
	}

	throw std::exception();

done:
	return i;
}

void drift(float64_t t, Hvu8& mask, Hvf64& mu, Hvf64_3& r, Hvf64_3& v, size_t start, size_t n)
{
	// save initial values
	Hvf64_3 r0 = r;
	Hvf64_3 v0 = v;

	Hvf64 dist(n);
	Hvf64 vsq(n);
	Hvf64 vdotr(n);
	for (size_t i = start; i < start + n; i++)
	{
		dist[i] = std::sqrt(r[i].lensq());
		vsq[i] = v[i].lensq();
		vdotr[i] = v0[i].x * r0[i].x + v0[i].y * r0[i].y + v0[i].z * r0[i].z;
	}

	Hvf64 energy = std::move(vsq);
	// vsq dies!
	for (size_t i = start; i < start + n; i++)
	{
		energy[i] *= 0.5;
		energy[i] -= mu[i] / dist[i];
	}

	for (size_t i = start; i < start + n; i++)
	{
		if (mask[i]) continue;
		if (energy[i] >= 0)
		{
			std::cerr << "unbound orbit" << std::endl;
			throw std::exception();
		}
		else
		{
			// maybe parallelize this
			float64_t a = -0.5 * mu[i] / energy[i];
			float64_t n_ = std::sqrt(mu[i] / (a * a * a));
			float64_t ecosEo = 1.0 - dist[i] / a;
			float64_t esinEo = vdotr[i] / (n_ * a * a);
			float64_t e = std::sqrt(ecosEo * ecosEo + esinEo * esinEo);

			// subtract off an integer multiple of complete orbits
			float64_t dM = t * n_ - M_2PI * (int) (t * n_ / M_2PI);

			// remaining time to advance
			float64_t dt = dM / n_;

			// call kepler equation solver with initial guess in dE already
			float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
			float64_t sindE, cosdE;
			kepeq(dM, ecosEo, esinEo, &dE, &sindE, &cosdE);

			float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
			float64_t f = 1.0 + a * (cosdE - 1.0) / dist[i];
			float64_t g = dt + (sindE - dE) / n_;
			float64_t fdot = -n_ * sindE * a / (dist[i] * fp);
			float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

			r[i] = r0[i] * f + v0[i] * g;
			v[i] = r0[i] * fdot + v0[i] * gdot;
		}
	}
}

void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa)
{
	(void) pa;

	helio_to_jacobi_r_planets(pl);
	helio_to_jacobi_v_planets(pl);
}

void first_step(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, float64_t dt)
{
	helio_acc_planets(pl, 0);
	helio_acc_particles(pl, pa, 0, 0);

	for (size_t i = 1; i < pl.n; i++)
	{
		pl.v[i] += pl.a[i] * (dt / 2);
	}

	for (size_t i = 0; i < pa.n; i++)
	{
		pa.v[i] += pa.a[i] * (dt / 2);
	}
}

void step_particles_cuda(const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa, size_t tbsize, float64_t dt)
{
	thrust::for_each(pa.mvs_begin(), pa.mvs_begin() + pa.n_alive, MVSKernel(pl, tbsize, dt));
}

void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, float64_t t, size_t index, float64_t dt)
{
	Hvf64 mu(pa.n_alive);
	Hvu8 mask(pa.n_alive);
	for (size_t i = 0; i < pa.n_alive; i++)
	{
		mu[i] = pl.m[0];
		mask[i] = pa.flags[i] > 0;
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(dt, mask, mu, pa.r, pa.v, 0, pa.n_alive);

	// find the accelerations of the heliocentric velocities
	helio_acc_particles(pl, pa, t, index);

	for (size_t i = 0; i < pa.n_alive; i++)
	{
		pa.v[i] += pa.a[i] * dt;
	}
}

void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t index, float64_t dt)
{
	(void) t;

	// Convert the heliocentric velocities to Jacobi velocities 
	helio_to_jacobi_v_planets(pl);

	Hvf64 mu(pl.n);
	Hvu8 mask(pl.n);
	for (size_t i = 1; i < pl.n; i++)
	{
		// Each Jacobi Kepler problem has a different central mass
		mu[i] = pl.m[0] * pl.eta[i] / pl.eta[i - 1];
		mask[i] = 0;
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(dt, mask, mu, pl.rj, pl.vj, 1, pl.n - 1);

	// convert Jacobi vectors to helio. ones for acceleration calc 
	jacobi_to_helio_planets(pl);

	// find the accelerations of the heliocentric velocities
	helio_acc_planets(pl, index);
	std::copy(pl.r.begin(), pl.r.end(), pl.r_log.begin() + (pl.n - 1) * index);

	for (size_t i = 1; i < pl.n; i++)
	{
		pl.v[i] += pl.a[i] * dt;
	}
}

double energy_planets(const HostPlanetPhaseSpace& p)
{
	double ke = 0.0;
	double pe = 0.0;

	for (size_t i = 0; i < p.n; i++)
	{
		ke += 0.5 * (p.v[i].x * p.v[i].x + p.v[i].y * p.v[i].y + p.v[i].z * p.v[i].z) * p.m[i];
	}

	for (size_t i = 0; i < p.n - 1; i++)
	{
		for (size_t j = i + 1; j < p.n; j++)
		{
			double dx = p.r[i].x - p.r[j].x;
			double dy = p.r[i].y - p.r[j].y;
			double dz = p.r[i].z - p.r[j].z;

			pe -= p.m[i] * p.m[j] / std::sqrt(dx * dx + dy * dy + dz * dz);
		}
	}

	return ke + pe;
}

f64_3 l_planets(const HostPlanetPhaseSpace& p)
{
	f64_3 l(0);

	for (size_t i = 0; i < p.n; i++)
	{
		l += p.r[i].cross(p.v[i]) * p.m[i];
	}

	return l;
}
