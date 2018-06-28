#include "wh.h"
#include "convert.h"

#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <iostream>

const size_t MAXKEP = 10;
const float64_t TOLKEP = 1E-14;

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

void calculate_planet_metrics(const HostPlanetPhaseSpace& p, double* energy, f64_3* l)
{
	f64_3 bary_r, bary_v;

	find_barycenter(p.r, p.v, p.m, p.n_alive, bary_r, bary_v);

	Vf64_3 r(p.n_alive), v(p.n_alive);
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

WHIntegrator::WHIntegrator() { }
WHIntegrator::WHIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config)
{
	size_t max = pl.n > pa.n ? pl.n : pa.n;

	inverse_helio_cubed = inverse_jacobi_cubed = Vf64(pl.n);
	dist = energy = vdotr = Vf64(max);
	mu = Vf64(max);
	mask = Vu8(max);
	eta = Vf64(pl.n);
	planet_rj = planet_vj = Vf64_3(pl.n);
	planet_a = Vf64_3(pl.n);
	particle_a = Vf64_3(pa.n);
	planet_rh = Vf64(pl.n);

	encounter_n1 = config.wh_ce_n1;
	encounter_n2 = config.wh_ce_n2;
	encounter_r1 = config.wh_ce_r1;
	encounter_r2 = config.wh_ce_r2;

	tbsize = config.tbsize;
	resolve_encounters = config.resolve_encounters;
	dt = config.dt;

	size_t ce_factor = resolve_encounters ? encounter_n1 * encounter_n2 : 1;
	planet_h0_log = planet_h0_log_old = Vf64_3(tbsize * ce_factor);

	if (resolve_encounters)
	{
		planet_h0_log_slow = planet_h0_log_slow_old = Vf64_3(tbsize);
	}


	eta[0] = pl.m[0];
	for (size_t i = 1; i < pl.n; i++)
	{
		eta[i] = eta[i - 1] + pl.m[i];
	}

	if (resolve_encounters)
	{
		planet_rh[0] = 0.5;
		for (size_t i = 1; i < pl.n; i++)
		{
			double a, e;
			to_elements(pl.m[0] + pl.m[i], pl.r[i], pl.v[i], nullptr, &a, &e);

			planet_rh[i] = a * (1 - e) * std::pow(pl.m[i] / (3 * pl.m[0]), 1. / 3);
			std::cout << planet_rh[i] << std::endl;
		}
	}
	else
	{
		encounter_r2 = 1;
		for (size_t i = 0; i < pl.n; i++)
		{
			planet_rh[i] = config.cull_radius;
		}
	}

	helio_to_jacobi_r_planets(pl, eta, planet_rj);
	helio_to_jacobi_v_planets(pl, eta, planet_vj);

	helio_acc_planets(pl, 0);
	helio_acc_particles(pl, pa, 0, pa.n_alive, 0, 0);
}

void WHIntegrator::integrate_planets_timeblock(HostPlanetPhaseSpace& pl, float64_t t)
{
	std::swap(planet_h0_log, planet_h0_log_old);
	std::swap(planet_h0_log_slow, planet_h0_log_slow_old);

	pl.n_alive = pl.n_alive_old;

	if (resolve_encounters)
	{
		size_t fast_factor = encounter_n1 * encounter_n2;

		for (size_t i = 0; i < tbsize * fast_factor; i++)
		{
			step_planets(pl, t, i);
			// take the planet positions at the end of every timestep

			if ((i + 1) % (fast_factor) == 0)
			{
				size_t slow_index = i / fast_factor;

				auto fast_begin = pl.r_log.begin() + i * (pl.n - 1);
				std::copy(fast_begin, fast_begin + (pl.n - 1), pl.r_log_slow.begin() + slow_index * (pl.n - 1));

				fast_begin = pl.v_log.begin() + i * (pl.n - 1);
				std::copy(fast_begin, fast_begin + (pl.n - 1), pl.v_log_slow.begin() + slow_index * (pl.n - 1));

				planet_h0_log_slow[i / fast_factor] = planet_h0_log[i];
			}

			t += dt / static_cast<double>(fast_factor);
		}
	}
	else
	{
		for (size_t i = 0; i < tbsize; i++)
		{
			step_planets(pl, t, i);
			t += dt;
		}
	}
}

void WHIntegrator::integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t)
{
	for (size_t i = 0; i < tbsize; i++)
	{
		step_particles(pl, pa, begin, length, t, i);
		t += dt;
	}
}

template<bool old>
void WHIntegrator::nonhelio_acc_encounter_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, float64_t time, size_t timestep_index, size_t central_planet_index)
{
	particle_a[particle_index] = f64_3(0);

	for (size_t i = 0; i < old ? pl.n_alive_old : pl.n_alive; i++)    
	{
		if (i == central_planet_index) continue;

		f64_3 dr = pl.r_log_at<false, old>(timestep_index, i) - pl.r_log_at<false, old>(timestep_index, central_planet_index);
		float64_t rsq = dr.lensq();
		double inverse_helio_cubed = 1. / (std::sqrt(rsq) * rsq);

#pragma GCC warning "TODO"
		particle_a[i] -= dr * pl.m[i] * inverse_helio_cubed;
        }

	for (size_t j = 0; j < old ? pl.n_alive_old : pl.n_alive; j++)
	{
		if (j == central_planet_index) continue;

		f64_3 dr = pa.r[particle_index] - pl.r_log_at<false, old>(timestep_index, j);
		float64_t planet_rji2 = dr.lensq();
		float64_t irij3 = 1. / (planet_rji2 * std::sqrt(planet_rji2));
		float64_t fac = pl.m[j] * irij3;

		particle_a[particle_index] -= dr * fac;
	}

	float64_t planet_rji2 = pa.r[particle_index].lensq();
	if (planet_rji2 > 2000 * 2000)
	{
		pa.deathtime[particle_index] = static_cast<float>(time);
		pa.deathflags[particle_index] = pa.deathflags[particle_index] | 0x0002;
	}
}

template<bool encounter, bool old>
void WHIntegrator::helio_acc_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, float64_t time, size_t timestep_index)
{
	particle_a[particle_index] = h0_log_at<!encounter, old>(timestep_index);

	for (size_t j = 1; j < pl.n_alive; j++)
	{
		f64_3 dr = pa.r[particle_index] - pl.r_log_at<!encounter, old>(timestep_index, j);
		float64_t planet_rji2 = dr.lensq();
		float64_t irij3 = 1. / (planet_rji2 * std::sqrt(planet_rji2));
		float64_t fac = pl.m[j] * irij3;

		particle_a[particle_index] -= dr * fac;

		if (planet_rji2 < planet_rh[j] * planet_rh[j] * encounter_r2 * encounter_r2)
		{
			pa.deathtime[particle_index] = static_cast<float>(time);
			pa.deathflags[particle_index] = static_cast<uint16_t>(pa.deathflags[particle_index] | (j << 8) | 0x0001);
		}
	}

	float64_t planet_rji2 = pa.r[particle_index].lensq();
	if (planet_rji2 < planet_rh[0] * planet_rh[0] * encounter_r2 * encounter_r2)
	{
		pa.deathtime[particle_index] = static_cast<float>(time);
		pa.deathflags[particle_index] = pa.deathflags[particle_index] | 0x0001;
	}
	if (planet_rji2 > 2000 * 2000)
	{
		pa.deathtime[particle_index] = static_cast<float>(time);
		pa.deathflags[particle_index] = pa.deathflags[particle_index] | 0x0002;
	}
}

void WHIntegrator::helio_acc_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t time, size_t timestep_index)
{
	for (size_t i = begin; i < begin + length; i++)
	{
		helio_acc_particle<false, false>(pl, pa, i, time, timestep_index);
	}
}

void WHIntegrator::helio_acc_planets(HostPlanetPhaseSpace& p, size_t index)
{
	for (size_t i = 1; i < p.n_alive; i++)
	{
		float64_t r2 = p.r[i].lensq();
		inverse_helio_cubed[i] = 1. / (std::sqrt(r2) * r2);
		r2 = this->planet_rj[i].lensq();
		inverse_jacobi_cubed[i] = 1. / (std::sqrt(r2) * r2);
        }
	
        // compute common heliocentric acceleration
	f64_3 a_common(0);
	for (size_t i = 2; i < p.n_alive; i++)    
	{
		float64_t mfac = p.m[i] * this->inverse_helio_cubed[i];
		a_common -= p.r[i] * mfac;
        }

        // Load this into all the arrays
	for (size_t i = 1; i < p.n_alive; i++)    
	{
		planet_a[i] = a_common;
        }

	planet_h0_log[index] = a_common - p.r[1] * p.m[1] * this->inverse_helio_cubed[1];
	
	// Now do indirect acceleration ; note that planet 1 does not receive a contribution 
	for (size_t i = 2; i < p.n_alive; i++)    
	{
		planet_a[i] += (this->planet_rj[i] * this->inverse_jacobi_cubed[i] - p.r[i] * this->inverse_helio_cubed[i]) * p.m[0];
        }
	
	/* next term ; again, first planet does not participate */
	f64_3 a_accum(0);
	for (size_t i = 2; i < p.n_alive; i++)    
	{
		float64_t mfac = p.m[i] * p.m[0] * this->inverse_jacobi_cubed[i] / eta[i-1];
		a_accum += this->planet_rj[i] * mfac;
		planet_a[i] += a_accum;
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
			planet_a[j] -= dr * mfac;

			// acc. on i is just negative, with m[j] instead
			mfac = p.m[j] * irij3;
			planet_a[i] += dr * mfac;
		}
	}
}

void WHIntegrator::drift_single(float64_t t, float64_t mu, f64_3* r, f64_3* v) const
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
		std::ostringstream ss;
		ss << "unbound orbit of particle, energy = " << energy << std::endl;
		throw std::runtime_error(ss.str());
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
		// float64_t e = std::sqrt(ecosEo * ecosEo + esinEo * esinEo);

		// subtract off an integer multiple of complete orbits
		float64_t dM = t * n_ - M_2PI * (int) (t * n_ / M_2PI);

		// remaining time to advance
		float64_t adv_dt = dM / n_;

		// call kepler equation solver with initial guess in dE already
		float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
		float64_t sindE, cosdE;
		kepeq(dM, ecosEo, esinEo, &dE, &sindE, &cosdE);

		float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
		float64_t f = 1.0 + a * (cosdE - 1.0) / dist;
		float64_t g = adv_dt + (sindE - dE) / n_;
		float64_t fdot = -n_ * sindE * a / (dist * fp);
		float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

		*r = r0 * f + v0 * g;
		*v = r0 * fdot + v0 * gdot;
	}
}

void WHIntegrator::drift(float64_t t, const Vu8& mask, const Vf64& mu, Vf64_3& r, Vf64_3& v, size_t start, size_t n)
{
	for (size_t i = start; i < start + n; i++)
	{
		this->dist[i] = std::sqrt(r[i].lensq());
		this->energy[i] = v[i].lensq();
		this->vdotr[i] = v[i].x * r[i].x + v[i].y * r[i].y + v[i].z * r[i].z;
	}

	for (size_t i = start; i < start + n; i++)
	{
		this->energy[i] *= 0.5;
		this->energy[i] -= mu[i] / this->dist[i];
	}

	for (size_t i = start; i < start + n; i++)
	{
		if (mask[i]) continue;
		if (this->energy[i] >= 0)
		{
			std::ostringstream ss;
			ss << "unbound orbit of planet " << i << " energy = " << this->energy[i] << std::endl;

			for (size_t j = start; j < start + n; j++)
			{
				ss << "p " << r[j].x << " " << r[j].y << " " << r[j].z << std::endl;
				ss << "v " << v[j].x << " " << v[j].y << " " << v[j].z << std::endl;
			}
			
			throw std::runtime_error(ss.str());
		}
		else
		{
			f64_3 r0 = r[i];
			f64_3 v0 = v[i];

			// maybe parallelize this
			float64_t a = -0.5 * mu[i] / this->energy[i];
			float64_t n_ = std::sqrt(mu[i] / (a * a * a));
			float64_t ecosEo = 1.0 - this->dist[i] / a;
			float64_t esinEo = this->vdotr[i] / (n_ * a * a);
			// float64_t e = std::sqrt(ecosEo * ecosEo + esinEo * esinEo);

			// subtract off an integer multiple of complete orbits
			float64_t dM = t * n_ - M_2PI * (int) (t * n_ / M_2PI);

			// remaining time to advance
			float64_t adv_dt = dM / n_;

			// call kepler equation solver with initial guess in dE already
			float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
			float64_t sindE, cosdE;
			kepeq(dM, ecosEo, esinEo, &dE, &sindE, &cosdE);

			float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
			float64_t f = 1.0 + a * (cosdE - 1.0) / this->dist[i];
			float64_t g = adv_dt + (sindE - dE) / n_;
			float64_t fdot = -n_ * sindE * a / (this->dist[i] * fp);
			float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

			r[i] = r0 * f + v0 * g;
			v[i] = r0 * fdot + v0 * gdot;
		}
	}
}

void WHIntegrator::step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, size_t timestep_index)
{
	for (size_t i = begin; i < begin + length; i++)
	{
		this->mask[i] = !((pa.deathflags[i] & 0x0001) || (pa.deathflags[i] == 0));
		if (!this->mask[i]) pa.v[i] += particle_a[i] * (dt / 2);
	}

	for (size_t i = begin; i < begin + length; i++)
	{
		this->mu[i] = pl.m[0];
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(dt, this->mask, this->mu, pa.r, pa.v, begin, length);

	// find the accelerations of the heliocentric velocities
	helio_acc_particles(pl, pa, begin, length, t, timestep_index);

	for (size_t i = begin; i < begin + length; i++)
	{
		if (!this->mask[i]) pa.v[i] += particle_a[i] * (dt / 2);
	}
}

template<bool old>
void WHIntegrator::integrate_encounter_particle_step(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t planet_index, size_t timestep_index, double t)
{
	size_t time_factor = encounter_n1 * encounter_n2;
	f64_3 planet_rel = pa.r[particle_index] -
		pl.r_log_at<true, old>(timestep_index, planet_index);
	double rsq = planet_rel.lensq();

	double rh = planet_rh[planet_index];

	if (rsq < rh * rh * encounter_r2 * encounter_r2)
	{
		if (rsq < rh * rh * encounter_r1 * encounter_r1)
		{
			pa.deathflags[particle_index] |= 0x0008;

			if (false)
			{
			for (size_t i = 0; i < time_factor; i++)
			{
				size_t log_index = timestep_index * time_factor + i;
					
				double little_dt = dt / static_cast<double>(time_factor);
				pa.v[particle_index] += particle_a[particle_index] * (little_dt / 2);

				f64_3 r_logged = pl.r_log_at<false, old>(log_index, planet_index);
				f64_3 v_logged = pl.v_log_at<false, old>(log_index, planet_index);

				f64_3 rel_r = pa.r[particle_index] - r_logged;
				f64_3 rel_v = pa.v[particle_index] - v_logged;

				// Drift all the particles along their Jacobi Kepler ellipses
				drift_single(little_dt, pl.m[planet_index], &rel_r, &rel_v);
				pa.r[particle_index] = rel_r + r_logged;
				pa.v[particle_index] = rel_v + v_logged;

				// find the accelerations of the heliocentric velocities
				nonhelio_acc_encounter_particle<old>(pl, pa, particle_index, t, log_index, planet_index);

				pa.v[particle_index] += particle_a[particle_index] * (little_dt / 2);
			}
			}
		}
		else
		{
			for (size_t i = 0; i < encounter_n1; i++)
			{
				double little_dt = dt / static_cast<double>(encounter_n1);
				pa.v[particle_index] += particle_a[particle_index] * (little_dt / 2);

				// Drift all the particles along their Jacobi Kepler ellipses
				drift_single(little_dt, pl.m[0], &pa.r[particle_index], &pa.v[particle_index]);

				// find the accelerations of the heliocentric velocities
				helio_acc_particle<true, old>(pl, pa, particle_index, t, timestep_index * time_factor + i * encounter_n2);

				pa.v[particle_index] += particle_a[particle_index] * (little_dt / 2);
			}
		}
	}
	else
	{
		pa.v[particle_index] += particle_a[particle_index] * (dt / 2);

		// Drift all the particles along their Jacobi Kepler ellipses
		drift_single(dt, pl.m[0], &pa.r[particle_index], &pa.v[particle_index]);

		// find the accelerations of the heliocentric velocities
		helio_acc_particle<false, old>(pl, pa, particle_index, t, timestep_index);

		pa.v[particle_index] += particle_a[particle_index] * (dt / 2);
	}
}

void WHIntegrator::integrate_encounter_particle_catchup(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t particle_deathtime_index, size_t planet_index)
{
	// Umm...
	if (planet_index == 0)
	{
		pa.deathflags[particle_index] = 0x0080;
		return;
	}

#pragma GCC warning "TODO"
	double t = 0;

	for (size_t i = particle_deathtime_index; i < tbsize; i++)
	{
		integrate_encounter_particle_step<true>(pl, pa, particle_index, planet_index, i, t);
		t += dt;
	}

	for (size_t i = 0; i < tbsize; i++)
	{
		integrate_encounter_particle_step<false>(pl, pa, particle_index, planet_index, i, t);
		t += dt;
	}
}

void WHIntegrator::gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length)
{
	gather(particle_a, indices, begin, length);
}

void WHIntegrator::step_planets(HostPlanetPhaseSpace& pl, float64_t t, size_t timestep_index)
{
	(void) t;

	size_t fast_factor = encounter_n1 * encounter_n2;
	double new_dt = resolve_encounters ? dt / static_cast<double>(fast_factor) : dt;

	for (size_t i = 1; i < pl.n_alive; i++)
	{
		pl.v[i] += planet_a[i] * (new_dt / 2);
	}

	// Convert the heliocentric velocities to Jacobi velocities 
	helio_to_jacobi_v_planets(pl, eta, planet_vj);

	for (size_t i = 1; i < pl.n_alive; i++)
	{
		// Each Jacobi Kepler problem has a different central mass
		this->mu[i] = pl.m[0] * eta[i] / eta[i - 1];
		this->mask[i] = 0;
        }

	// Drift all the particles along their Jacobi Kepler ellipses
	drift(new_dt, this->mask, this->mu, this->planet_rj, this->planet_vj, 1, pl.n_alive - 1);

	// convert Jacobi vectors to helio. ones for acceleration calc 
	jacobi_to_helio_planets(eta, planet_rj, planet_vj, pl);

	// find the accelerations of the heliocentric velocities
	helio_acc_planets(pl, timestep_index);
	auto& r_log = pl.r_log_slow;
	auto& v_log = pl.v_log_slow;
	if (resolve_encounters)
	{
		r_log = pl.r_log;
		v_log = pl.v_log;
	}

	std::copy(pl.r.begin() + 1, pl.r.end(), r_log.begin() + (pl.n_alive - 1) * timestep_index);
	std::copy(pl.v.begin() + 1, pl.v.end(), v_log.begin() + (pl.n_alive - 1) * timestep_index);

	for (size_t i = 1; i < pl.n_alive; i++)
	{
		pl.v[i] += planet_a[i] * (new_dt / 2);
	}
}
