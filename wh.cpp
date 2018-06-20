#include "wh.h"
#include "convert.h"

#include <iostream>
#include <iomanip>
#include <cmath>

const size_t MAXKEP = 10;
const float64_t TOLKEP = 1E-14;

void helio_acc_particle_ce(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, size_t i, float64_t time, size_t timestep_index)
{
	p.a[i] = pl.h0_log[timestep_index];

	for (size_t j = 1; j < pl.n_alive; j++)
	{
		f64_3 dr = p.r[i] - pl.r[j];
		float64_t rji2 = dr.lensq();
		float64_t irij3 = 1. / (rji2 * std::sqrt(rji2));
		float64_t fac = pl.m[j] * irij3;

		p.a[i] -= dr * fac;
	}

	float64_t rji2 = p.r[i].lensq();
	if (rji2 > 200 * 200)
	{
		p.deathtime[i] = time;
		p.deathflags[i] = p.deathflags[i] | 0x0002;
	}
}

void helio_acc_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, size_t begin, size_t length, float64_t time, size_t timestep_index)
{
	for (size_t i = begin; i < begin + length; i++)
	{
		p.a[i] = pl.h0_log[timestep_index];

		for (size_t j = 1; j < pl.n_alive; j++)
		{
			f64_3 dr = p.r[i] - pl.r[j];
			float64_t rji2 = dr.lensq();
			float64_t irij3 = 1. / (rji2 * std::sqrt(rji2));
			float64_t fac = pl.m[j] * irij3;

			p.a[i] -= dr * fac;

			if (rji2 < 0.5 * 0.5)
			{
				p.deathtime[i] = time;
				p.deathflags[i] = p.deathflags[i] | (j << 8) | 0x0001;
			}
		}

		float64_t rji2 = p.r[i].lensq();
		if (rji2 < 0.5 * 0.5)
		{
			p.deathtime[i] = time;
			p.deathflags[i] = p.deathflags[i] | 0x0001;
		}
		if (rji2 > 200 * 200)
		{
			p.deathtime[i] = time;
			p.deathflags[i] = p.deathflags[i] | 0x0002;
		}
	}
}

void helio_acc_planets(HostPlanetPhaseSpace& p, size_t index)
{
	Vf64 inverse_helio_cubed(p.n_alive), inverse_jacobi_cubed(p.n_alive);

	for (size_t i = 1; i < p.n_alive; i++)
	{
		float64_t r2 = p.r[i].lensq();
		inverse_helio_cubed[i] = 1. / (std::sqrt(r2) * r2);
		r2 = p.rj[i].lensq();
		inverse_jacobi_cubed[i] = 1. / (std::sqrt(r2) * r2);
        }
	
        // compute common heliocentric acceleration
	f64_3 a_common(0);
	for (size_t i = 2; i < p.n_alive; i++)    
	{
		float64_t mfac = p.m[i] * inverse_helio_cubed[i];
		a_common -= p.r[i] * mfac;
        }

        // Load this into all the arrays
	for (size_t i = 1; i < p.n_alive; i++)    
	{
		p.a[i] = a_common;
        }

	p.h0_log[index] = a_common - p.r[1] * p.m[1] * inverse_helio_cubed[1];
	
	// Now do indirect acceleration ; note that planet 1 does not receive a contribution 
	for (size_t i = 2; i < p.n_alive; i++)    
	{
		p.a[i] += (p.rj[i] * inverse_jacobi_cubed[i] - p.r[i] * inverse_helio_cubed[i]) * p.m[0];
        }
	
	/* next term ; again, first planet does not participate */
	f64_3 a_accum(0);
	for (size_t i = 2; i < p.n_alive; i++)    
	{
		float64_t mfac = p.m[i] * p.m[0] * inverse_jacobi_cubed[i] / p.eta[i-1];
		a_accum += p.rj[i] * mfac;
		p.a[i] += a_accum;
        }

	/* Finally, incorporate the direct accelerations */
	for (size_t i = 1; i < p.n_alive - 1; i++)    
	{
		for (size_t j = i + 1; j < p.n_alive; j++)    
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

void drift_single(float64_t t, float64_t mu, f64_3* r, f64_3* v)
{
	float64_t dist, vsq, vdotr;
	dist = std::sqrt(r->lensq());
	vsq = v->lensq();
	vdotr = v->x * r->x + v->y * r->y + v->z * r->z;

	float64_t energy = vsq;
	energy *= 0.5;
	energy -= mu / dist;

	if (energy >= 0)
	{
		// TODO
		std::cerr << "unbound orbit of particle, energy = " << energy << std::endl;
		throw std::exception();
	}
	else
	{
		f64_3 r0 = *r;
		f64_3 v0 = *v;

		// maybe parallelize this
		float64_t a = -0.5 * mu / energy;
		float64_t n_ = std::sqrt(mu / (a * a * a));
		float64_t ecosEo = 1.0 - dist / a;
		float64_t esinEo = vdotr / (n_ * a * a);
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
		float64_t f = 1.0 + a * (cosdE - 1.0) / dist;
		float64_t g = dt + (sindE - dE) / n_;
		float64_t fdot = -n_ * sindE * a / (dist * fp);
		float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

		*r = r0 * f + v0 * g;
		*v = r0 * fdot + v0 * gdot;
	}
}

void drift(float64_t t, Vu8& mask, Vf64& mu, Vf64_3& r, Vf64_3& v, size_t start, size_t n)
{
	Vf64 dist(n);
	Vf64 vsq(n);
	Vf64 vdotr(n);
	for (size_t i = start; i < start + n; i++)
	{
		dist[i] = std::sqrt(r[i].lensq());
		vsq[i] = v[i].lensq();
		vdotr[i] = v[i].x * r[i].x + v[i].y * r[i].y + v[i].z * r[i].z;
	}

	Vf64 energy = std::move(vsq);
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
			std::cerr << "unbound orbit of planet " << i << " energy = " << energy[i] << std::endl;

			for (size_t j = start; j < start + n; j++)
			{
				std::cerr << "p " << r[j].x << " " << r[j].y << " " << r[j].z << std::endl;
				std::cerr << "v " << v[j].x << " " << v[j].y << " " << v[j].z << std::endl;
			}
			
			throw std::exception();
		}
		else
		{
			f64_3 r0 = r[i];
			f64_3 v0 = v[i];

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

			r[i] = r0 * f + v0 * g;
			v[i] = r0 * fdot + v0 * gdot;
		}
	}
}

void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa)
{
	(void) pa;

	helio_to_jacobi_r_planets(pl);
	helio_to_jacobi_v_planets(pl);

	helio_acc_planets(pl, 0);
	helio_acc_particles(pl, pa, 0, pa.n_alive, 0, 0);
}

void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, size_t timestep_index, float64_t dt)
{
	Vu8 mask(length);
	for (size_t i = begin; i < begin + length; i++)
	{
		mask[i] = !((pa.deathflags[i] & 0x0001) || (pa.deathflags[i] == 0));
		if (!mask[i]) pa.v[i] += pa.a[i] * (dt / 2);
	}

	Vf64 mu(length);
	for (size_t i = begin; i < begin + length; i++)
	{
		mu[i] = pl.m[0];
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(dt, mask, mu, pa.r, pa.v, begin, length);

	// find the accelerations of the heliocentric velocities
	helio_acc_particles(pl, pa, begin, length, t, timestep_index);

	for (size_t i = begin; i < begin + length; i++)
	{
		if (!mask[i]) pa.v[i] += pa.a[i] * (dt / 2);
	}
}

void step_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, float64_t t, size_t timestep_index, float64_t dt)
{
	size_t i = particle_index;
	if ((pa.deathflags[i] & 0x0001) || (pa.deathflags[i] == 0))
	{
		return;
	}

	pa.v[i] += pa.a[i] * (dt / 2);

	// Drift all the particles along their Jacobi Kepler ellipses
	drift_single(dt, pl.m[0], &pa.r[i], &pa.v[i]);

	// find the accelerations of the heliocentric velocities
	helio_acc_particle_ce(pl, pa, i, t, timestep_index);

	pa.v[i] += pa.a[i] * (dt / 2);
}

void step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t index, float64_t dt)
{
	(void) t;

	for (size_t i = 1; i < pl.n_alive; i++)
	{
		pl.v[i] += pl.a[i] * (dt / 2);
	}

	// Convert the heliocentric velocities to Jacobi velocities 
	helio_to_jacobi_v_planets(pl);

	Vf64 mu(pl.n_alive);
	Vu8 mask(pl.n_alive);
	for (size_t i = 1; i < pl.n_alive; i++)
	{
		// Each Jacobi Kepler problem has a different central mass
		mu[i] = pl.m[0] * pl.eta[i] / pl.eta[i - 1];
		mask[i] = 0;
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(dt, mask, mu, pl.rj, pl.vj, 1, pl.n_alive - 1);

	// convert Jacobi vectors to helio. ones for acceleration calc 
	jacobi_to_helio_planets(pl);

	// find the accelerations of the heliocentric velocities
	helio_acc_planets(pl, index);
	std::copy(pl.r.begin() + 1, pl.r.end(), pl.r_log.begin() + (pl.n_alive - 1) * index);

	for (size_t i = 1; i < pl.n_alive; i++)
	{
		pl.v[i] += pl.a[i] * (dt / 2);
	}
}

void calculate_planet_metrics(const HostPlanetPhaseSpace& p, double* energy, f64_3* l)
{
	Vf64_3 r(p.n_alive), v(p.n_alive);
	f64_3 bary_r, bary_v;

	find_barycenter(p.r, p.v, p.m, p.n_alive, bary_r, bary_v);

	for (size_t i = 0; i < p.n_alive; i++)
	{
		r[i] = p.r[i] - bary_r;
		v[i] = p.v[i] - bary_v;
	}

	if (energy)
	{
		double ke = 0.0;
		double pe = 0.0;

		for (size_t i = 0; i < p.n_alive; i++)
		{
			ke += 0.5 * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z) * p.m[i];
		}

		for (size_t i = 0; i < p.n_alive - 1; i++)
		{
			for (size_t j = i + 1; j < p.n_alive; j++)
			{
				double dx = r[i].x - r[j].x;
				double dy = r[i].y - r[j].y;
				double dz = r[i].z - r[j].z;

				pe -= p.m[i] * p.m[j] / std::sqrt(dx * dx + dy * dy + dz * dz);
			}
		}

		*energy = ke + pe;
	}

	if (l)
	{
		*l = f64_3(0);

		for (size_t i = 0; i < p.n_alive; i++)
		{
			*l += r[i].cross(v[i]) * p.m[i];
		}
	}
}
