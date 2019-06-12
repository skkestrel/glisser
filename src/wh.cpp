#include "wh.h"
#include "convert.h"
#include <iomanip>
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace sr
{
namespace wh
{
	const uint32_t MAXKEP = 10;
	const float64_t TOLKEP = 1E-14;

	using namespace sr::data;

	bool kepeq(double dM, double esinEo, double ecosEo, double* dE, double* sindE, double* cosdE, uint32_t* iterations)
	{
		double f,fp, delta;

		*sindE = std::sin( *dE);
		*cosdE = std::cos( *dE);

		for (uint32_t i = 0; i < MAXKEP; i++)
		{
			f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
			fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
			delta = -f / fp;
			if (std::fabs(delta) < TOLKEP)
			{
				*iterations = i;
				return false;
			}

			*dE += delta;
			*sindE = std::sin(*dE);
			*cosdE = std::cos(*dE);
		}

		return true;
	}


	bool kepeq_fixed(double dM, double esinEo, double ecosEo, double* dE, double* sindE, double* cosdE, uint32_t its)
	{
		double f, fp, delta = 1;

		*sindE = std::sin( *dE);
		*cosdE = std::cos( *dE);

		for (uint32_t i = 0; i < its; i++)
		{
			f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
			fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
			delta = -f / fp;

			*dE += delta;
			*sindE = std::sin(*dE);
			*cosdE = std::cos(*dE);
		}

		return std::fabs(delta) >= TOLKEP;
	}

	//  Subroutine for solving kepler's equation in difference form for an
	//  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
	//  for the criteria.
	//  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
	//  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
	//  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.
	//
	//	Input:
	//	    dm		==> increment in mean anomaly M (real*8 scalar)
	//	    es,ec       ==> ecc. times sin and cos of E_0 (real*8 scalars)
	//
	//       Output:
	//            x          ==> solution to Kepler's difference eqn (real*8 scalar)
	//            s,c        ==> sin and cosine of x (real*8 scalars)
	//
	void kepmd(double dm, double es, double ec, double* x, double* s, double* c)
	{
		const double A0 = 39916800.;
		const double A1 = 6652800.;
		const double A2 = 332640.;
		const double A3 = 7920.;
		const double A4 = 110.;

		// calc initial guess for root
		double fac1 = 1. / (1. - ec);
		double q = fac1 * dm;
		double fac2 = es * es * fac1 - ec / 3.;
		*x = q * (1. - 0.5 * fac1 * q * (es - q * fac2));

		// excellent approx. to sin and cos of x for small x.
		double y = *x * *x;
		*s = *x * (A0 - y * (A1 - y * (A2 - y * (A3 - y * (A4 - y))))) / A0;
		*c = std::sqrt(1. - *s * *s);

		// Compute better value for the root using quartic Newton method
		double f = *x - ec * *s + es * (1. - *c) - dm;
		double fp = 1. - ec * *c + es * *s;
		double fpp = ec * *s + es * *c;
		double fppp = ec * *c - es * *s;
		double dx = -f / fp;
		dx = -f / (fp + 0.5 * dx * fpp);
		dx = -f / (fp + 0.5 * dx * fpp + 0.16666666666666666 * dx * dx * fppp);
		*x = *x + dx;

		y = *x * *x;
		*s = *x * (A0 - y * (A1 - y * (A2 - y * (A3 - y * (A4 - y))))) / A0;
		*c = std::sqrt(1. - *s * *s);
	}

	// Returns the real root of cubic often found in solving kepler
	// problem in universal variables.
	//
	//             Input:
	//                 dt            ==>  time step (real scalar)
	//                 r0            ==>  Distance between `Sun' and paritcle
	//                                     (real scalar)
	//                 mu            ==>  Reduced mass of system (real scalar)
	//                 alpha         ==>  Twice the binding energy (real scalar)
	//                 u             ==>  Vel. dot radial vector (real scalar)
	//             Output:
	//                 s             ==>  solution of cubic eqn for the  
	//                                    universal variable
	//                 iflg          ==>  success flag ( = 0 if O.K.) (integer)
	//
	bool kepu_p3solve(double dt, double r0, double mu, double alpha, double u, double* s)
	{
		double denom = (mu - alpha * r0) / 6.;
		double a2 = 0.5 * u / denom;
		double a1 = r0 / denom;
		double a0 = -dt / denom;

		double q = (a1 - a2 * a2 / 3.) / 3.;
		double r = (a1 * a2 - 3. * a0) / 6. - (a2 * a2 * a2) / 27.;
		double sq2 = q * q * q + r * r;

		if (sq2 >= 0)
		{
			double sq = std::sqrt(sq2);
			double p1, p2;

			if (r + sq <= 0)
			{
				p1 = -std::pow(-r - sq, 1. / 3.);
			}
			else
			{
				p1 = std::pow(r + sq, 1. / 3.);
			}

			if (r - sq <= 0)
			{
				p2 = -std::pow(-r + sq, 1. / 3.);
			}
			else
			{
				p2 = std::pow(r - sq, 1. / 3.);
			}

			*s = p1 + p2 - a2 / 3.;
			return false;
		}
		else
		{
			*s = 0;
			return true;
		}
	}


	// Initial guess for solving kepler's equation using universal variables.
	//
	//             Input:
	//                 dt            ==>  time step (real scalor)
	//                 r0            ==>  Distance between `Sun' and paritcle
	//                                     (real scalor)
	//                 mu            ==>  Reduced mass of system (real scalor)
	//                 alpha         ==>  energy (real scalor)
	//                 u             ==>  angular momentun  (real scalor)
	//             Output:
	//                 s             ==>  initial guess for the value of 
	//                                    universal variable
	void kepu_guess(double dt, double r0, double mu, double alpha, double u, double* s)
	{
		if (alpha > 0)
		{
			// find initial guess for elliptic motion
			if (dt / r0 <= 0.4)
			{
				*s = dt / r0 - (dt * dt * u) / (2. * r0 * r0 * r0);
			}
			else
			{
				double a = mu / alpha;
				double en = std::sqrt(mu / (a * a * a));
				double ec = 1.0 - r0 / a;
				double es = u / (en * a * a);
				double e = std::sqrt(ec * ec + es * es);
				double y = en * dt - es;
				double sy, cy;
				sincos(y, &sy, &cy);
				double sigma = es * cy + ec * sy > 0 ? 1 : -1;
				double x = y + sigma * .85 * e;
				*s = x / std::sqrt(alpha);
			}
		}
		else
		{
			// find initial guess for hyperbolic motion.
			if (kepu_p3solve(dt, r0, mu, alpha, u, s))
			{
				*s = dt / r0;
			}
		}
	}

	// subroutine for the calculation of stumpff functions
	// see Danby p.172  equations 6.9.15
	//
	//             Input:
	//                 x             ==>  argument
	//             Output:
	//                 c0,c1,c2,c3   ==>  c's from p171-172
	//                                       (real scalors)
	void kepu_stumpff(double x, double* c0, double* c1, double* c2, double* c3)
	{
	      uint32_t n = 0;
	      double xm = 0.1;
	      while (std::abs(x) >= xm)
	      {
		      n++;
		      x /= 4.;
	      }

	      // Huh?
	      *c2 = (1. - x * (1. - x * (1. - x * (1. - x * (1. - x * (1. - x / 182.) / 132.) / 90.) / 56.) / 30.) / 12.) / 2.;
	      *c3 = (1. - x * (1. - x * (1. - x * (1. - x * (1. - x * (1. - x / 210.) / 156.) / 110.) / 72.) / 42.) / 20.) / 6.;
	      *c1 = 1. - x * *c3;
	      *c0 = 1. - x * *c2;

	      if (n != 0)
	      {
		      for (uint32_t i = 0; i < n; i++)
		      {
			      *c3 = (*c2 + *c0 * *c3)/4.;
			      *c2 = *c1 * *c1 / 2.;
			      *c1 = *c0 * *c1;
			      *c0 = 2. * *c0 * *c0 - 1.;
			      x *= 4;
		      }
	      }
	}

	const double DANBYB = 1e-13;
	bool kepu_new(double s, double dt, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3)
	{
		for (uint8_t nc = 0; nc < 7; nc++)
		{
			double x = s * s * alpha;
			double c0;
			kepu_stumpff(x, &c0, c1, c2, c3);
			*c1 = *c1 * s;
			*c2 = *c2 * s * s;
			*c3 = *c3 * s * s * s;
			double f = r0 * *c1 + u * *c2 + mu * *c3 - dt;
			*fp = r0 * c0 + u * *c1 + mu * *c2;
			double fpp = (-r0 * alpha + mu) * *c1 + u * c0;
			double fppp = (- r0 * alpha + mu) * c0 - u * alpha * *c1;
			double ds = - f / *fp;
			ds = -f / (*fp + ds * fpp / 2.);
			ds = -f / (*fp + ds * fpp / 2. + ds * ds * fppp / 6.);
			s = s + ds;
			double fdt = f/dt;

			// quartic convergence
			// newton's method succeeded
			if (fdt * fdt < DANBYB * DANBYB)
			{
				return false;
			}
		}

		// newton's method failed
		return true;
	}

	// Returns the value of the function f of which we are trying to find the root
	// in universal variables.
	//
	//             Input:
	//                 dt            ==>  time step (real scalar)
	//                 r0            ==>  Distance between `Sun' and particle
	//                                     (real scalar)
	//                 mu            ==>  Reduced mass of system (real scalar)
	//                 alpha         ==>  Twice the binding energy (real scalar)
	//                 u             ==>  Vel. dot radial vector (real scalar)
	//                 s             ==>  Approx. root of f 
	//             Output:
	//                 f             ==>  function value
	//
	double kepu_fchk(double dt, double r0, double mu, double alpha, double u, double s)
	{
		double c0, c1, c2, c3;

		double x = s * s * alpha;
		kepu_stumpff(x, &c0, &c1, &c2, &c3);
		c1 = c1 * s;
		c2 = c2 * s * s;
		c3 = c3 * s * s * s;
		return r0 * c1 + u * c2 + mu * c3 - dt;
	}


	// subroutine for solving kepler's equation in universal variables.
	// using LAGUERRE'S METHOD
	//
	//             Input:
	//                 s             ==>  inital value of universal variable
	//                 dt            ==>  time step (real scalor)
	//                 r0            ==>  Distance between `Sun' and paritcle
	//                                     (real scalor)
	//                 mu            ==>  Reduced mass of system (real scalor)
	//                 alpha         ==>  energy (real scalor)
	//                 u             ==>  angular momentun  (real scalor)
	//             Output:
	//                 s             ==>  final value of universal variable
	//                 fp            ==>  f' from p170  
	//                                       (real scalors)
	//                 c1,c2,c3      ==>  c's from p171-172
	//                                       (real scalors)
	//                 iflgn          ==>  =0 if converged; !=0 if not
	const uint32_t NLAG1 = 50;
	const uint32_t NLAG2 = 400;
	const uint32_t NTMP = NLAG2 + 1;
	bool kepu_lag(double s, double dt, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3)
	{
		uint32_t ncmax;

		// To get close approch needed to take lots of iterations if alpha<0
		if (alpha < 0.0)
		{
			ncmax = NLAG2;
		}
		else
		{
			// Huh? Is this supposed to be NLAG1?
			// ncmax = NLAG1;
			ncmax = NLAG2;
		}

		double ln = 5.0;
		// start laguere's method
		for (uint32_t nc = 0; nc <= ncmax; nc++)
		{
			double c0;
			double x = s * s * alpha;
			kepu_stumpff(x, &c0, c1, c2, c3);
			*c1 = *c1 * s;
			*c2 = *c2 * s * s;
			*c3 = *c3 * s * s * s;
			double f = r0 * *c1 + u * *c2 + mu * *c3 - dt;
			*fp = r0 * c0 + u * *c1 + mu * *c2;
			double fpp = (-40. * alpha + mu) * *c1 + u * c0;
			double ds = -ln * f / (*fp + (*fp > 0 ? 1 : -1) * std::sqrt(std::abs((ln - 1.) * (ln - 1.) * *fp * *fp - (ln - 1.) * ln * f * fpp)));
			s = s + ds;

			double fdt = f / dt;

			// quartic convergence
			if (fdt * fdt < DANBYB * DANBYB)
			{
				return false;
			}
		}

		return true;
	}

	// subroutine for solving kepler's equation using universal variables.
	//
	//             Input:
	//                 dt            ==>  time step (real scalor)
	//                 r0            ==>  Distance between `Sun' and paritcle
	//                                     (real scalor)
	//                 mu            ==>  Reduced mass of system (real scalor)
	//                 alpha         ==>  energy (real scalor)
	//                 u             ==>  angular momentun  (real scalor)
	//             Output:
	//                 fp            ==>  f' from p170  
	//                                       (real scalors)
	//                 c1,c2,c3      ==>  c's from p171-172
	//                                       (real scalors)
	//                 iflg          ==>  =0 if converged; !=0 if not
	bool kepu(double dt, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3)
	{
		double s;
		kepu_guess(dt, r0, mu, alpha, u, &s);

		double st = s;
		// store initial guess for possible use later in
		// laguerre's method, in case newton's method fails.

		if (kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3))
		{
			double fo = kepu_fchk(dt, r0, mu, alpha, u, st);
			double fn = kepu_fchk(dt, r0, mu, alpha, u, s);
			if (std::abs(fo) < std::abs(fn))
			{
				s = st;
			}
			return kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3);
		}

		return false;
	}

	void calculate_planet_metrics(const HostPlanetPhaseSpace& pl, double* energy, f64_3* l)
	{
		f64_3 bary_r, bary_v;
		double bary_m;

		sr::convert::find_barycenter(pl.r(), pl.v(), pl.m(), pl.n_alive(), bary_r, bary_v, bary_m);

		Vf64_3 r(pl.n_alive()), v(pl.n_alive());
		for (size_t i = 0; i < pl.n_alive(); i++)
		{
			r[i] = pl.r()[i] - bary_r;
			v[i] = pl.v()[i] - bary_v;
		}

		if (energy)
		{
			double ke = 0.0;
			double pe = 0.0;

			for (size_t i = 0; i < pl.n_alive(); i++)
			{
				ke += 0.5 * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z) * pl.m()[i];
			}

			for (size_t i = 0; i < pl.n_alive() - 1; i++)
			{
				for (size_t j = i + 1; j < pl.n_alive(); j++)
				{
					double dx = r[i].x - r[j].x;
					double dy = r[i].y - r[j].y;
					double dz = r[i].z - r[j].z;

					pe -= pl.m()[i] * pl.m()[j] / std::sqrt(dx * dx + dy * dy + dz * dz);
				}
			}

			*energy = ke + pe;
		}

		if (l)
		{
			*l = f64_3(0);

			for (size_t i = 0; i < pl.n_alive(); i++)
			{
				*l += r[i].cross(v[i]) * pl.m()[i];
			}
		}
	}

	WHIntegrator::WHIntegrator() { }
	WHIntegrator::WHIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config)
	{
		planet_inverse_helio_cubed = planet_inverse_jacobi_cubed = Vf64(pl.n());
		planet_dist = planet_energy = planet_vdotr = Vf64(pl.n());
		particle_dist = particle_energy = particle_vdotr = Vf64(pa.n());

		planet_mu = Vf64(pl.n());
		particle_mu = Vf64(pa.n());

		planet_mask = Vu8(pl.n());
		particle_mask = Vu8(pa.n());

		planet_eta = Vf64(pl.n());
		planet_rj = planet_vj = Vf64_3(pl.n());
		planet_a = Vf64_3(pl.n());
		particle_a = Vf64_3(pa.n());

		planet_rh = Vf64(pl.n());

		planet_h0_log = sr::util::History<Vf64_3>(config.tbsize);

		outer_radius = config.outer_radius;

		planet_eta[0] = pl.m()[0];
		for (size_t i = 1; i < pl.n(); i++)
		{
			planet_eta[i] = planet_eta[i - 1] + pl.m()[i];
		}

		if (std::abs(config.encounter_sphere_factor) > 1e-6)
		{
			// reinterpret cull radius to be the sun's radius
			planet_rh[0] = config.cull_radius;

			// for the rest of the planets, use the factor times the hill sphere
			for (size_t i = 1; i < pl.n(); i++)
			{
				double a, e;
				sr::convert::to_elements(pl.m()[0] + pl.m()[i], pl.r()[i], pl.v()[i], nullptr, &a, &e);

				planet_rh[i] = config.encounter_sphere_factor * a * (1 - e) * std::pow(pl.m()[i] / (3 * pl.m()[0]), 1. / 3);
			}
		}
		else
		{
			for (size_t i = 0; i < pl.n(); i++)
			{
				planet_rh[i] = config.cull_radius;
			}
		}

		sr::convert::helio_to_jacobi_r_planets(pl, planet_eta, planet_rj);
		sr::convert::helio_to_jacobi_v_planets(pl, planet_eta, planet_vj);

		// copy the planet locations into first entry of the old log - this is used just so that helio_acc_particles works below
		std::copy(pl.r().begin() + 1, pl.r().end(), pl.r_log().old.begin());

		// calculate initial accelerations for particles - this must be done so that the GPU integration works
		
		// this function populates h0 for the first timestep
		helio_acc_planets(pl, 0);

		// use the old log since we copied planet locations into the old log above
		helio_acc_particles<true>(pl, pa, 0, pa.n_alive(), 0, 0);

		// detect initial encounters
		for (size_t i = 0; i < pa.n_alive(); i++)
		{
			for (size_t j = 0; j < pl.n_alive(); j++)
			{
				f64_3 dr = pa.r()[i] - pl.r()[j];

				double rji2 = dr.lensq();
				if (rji2 < planet_rh[j] * planet_rh[j])
				{
					pa.deathtime_map()[pa.id()[i]] = static_cast<float>(config.t_0);
					pa.deathflags()[i] = static_cast<uint16_t>((j << 8) | 0x0001);
				}
			}

			double rji2 = pa.r()[i].lensq();

			if (rji2 < planet_rh[0] * planet_rh[0])
			{
				pa.deathtime_map()[pa.id()[i]] = static_cast<float>(config.t_0);
				pa.deathflags()[i] = 0x0080;
			}

			if (rji2 > config.outer_radius * config.outer_radius)
			{
				pa.deathtime_map()[pa.id()[i]] = static_cast<float>(config.t_0);
				pa.deathflags()[i] = 0x0002;
			}
		}
	}

	void WHIntegrator::swap_logs()
	{
		planet_h0_log.swap_logs();
	}

	void WHIntegrator::integrate_planets_timeblock(HostPlanetPhaseSpace& pl, size_t nsteps, float64_t t, float64_t dt)
	{
		for (size_t i = 0; i < nsteps; i++)
		{
			step_planets(pl, t, dt, i);
			t += dt;
		}

		pl.r_log().len_old = nsteps;
		pl.v_log().len_old = nsteps;
		planet_h0_log.len_old = nsteps;
	}

	void WHIntegrator::integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, float64_t dt)
	{
		for (size_t i = begin; i < begin + length; i++)
		{
			this->particle_mu[i] = pl.m()[0];
		}

		for (size_t i = 0; i < pl.r_log().len; i++)
		{
			step_particles(pl, pa, begin, length, t, dt, i);
			t += dt;
		}
	}

	template<bool old>
	void WHIntegrator::helio_acc_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, float64_t time, size_t timestep_index)
	{
		f64_3& a = particle_a[particle_index];
		a = planet_h0_log.get<old>()[timestep_index];

		for (size_t j = 1; j < pl.n_alive(); j++)
		{
			f64_3 dr = pa.r()[particle_index] - pl.r_log().get<old>()[pl.log_index_at<old>(timestep_index, j)];
#ifdef USE_FMA
			float64_t planet_rji2 = std::fma(dr.x, dr.x, std::fma(dr.y, dr.y, dr.z * dr.z));
#else
			float64_t planet_rji2 = dr.lensq();
#endif

			float64_t irij3 = 1. / (planet_rji2 * std::sqrt(planet_rji2));
			float64_t fac = pl.m()[j] * irij3;

#ifdef USE_FMA
			a.x = std::fma(-dr.x, fac, a.x);
			a.y = std::fma(-dr.y, fac, a.y);
			a.z = std::fma(-dr.z, fac, a.z);
#else
			a -= dr * fac;
#endif

			if (planet_rji2 < planet_rh[j] * planet_rh[j])
			{
				pa.deathtime_map()[pa.id()[particle_index]] = static_cast<float>(time);
				pa.deathflags()[particle_index] = static_cast<uint16_t>((j << 8) | 0x0001);
			}
		}

		float64_t planet_rji2 = pa.r()[particle_index].lensq();

		if (planet_rji2 < planet_rh[0] * planet_rh[0])
		{
			pa.deathtime_map()[pa.id()[particle_index]] = static_cast<float>(time);
			pa.deathflags()[particle_index] = 0x0001;
		}

		if (planet_rji2 > outer_radius * outer_radius)
		{
			pa.deathtime_map()[pa.id()[particle_index]] = static_cast<float>(time);
			pa.deathflags()[particle_index] = 0x0002;
		}
	}

	template<bool old>
	void WHIntegrator::helio_acc_particles(
			const HostPlanetPhaseSpace& pl,
			HostParticlePhaseSpace& pa,
			size_t begin,
			size_t length,
			float64_t time,
			size_t timestep_index)
	{
		for (size_t i = begin; i < begin + length; i++)
		{
			helio_acc_particle<old>(pl, pa, i, time, timestep_index);
		}
	}

	void WHIntegrator::helio_acc_particles(
			const HostPlanetPhaseSpace& pl,
			HostParticlePhaseSpace& pa,
			size_t begin,
			size_t len,
			float64_t time,
			size_t timestep_index,
			bool old)
	{
		if (old)
		{
			helio_acc_particles<true>(pl, pa, begin, len, time, timestep_index);
		}
		else
		{
			helio_acc_particles<false>(pl, pa, begin, len, time, timestep_index);
		}

	}


	// Loads h0 into the old log at time given by index
	void WHIntegrator::helio_acc_planets(HostPlanetPhaseSpace& p, size_t index)
	{
		for (size_t i = 1; i < p.n_alive(); i++)
		{
			float64_t r2 = p.r()[i].lensq();
			this->planet_inverse_helio_cubed[i] =
				1. / (std::sqrt(r2) * r2);
			r2 = this->planet_rj[i].lensq();
			this->planet_inverse_jacobi_cubed[i] =
				1. / (std::sqrt(r2) * r2);
		}
		
		// compute common heliocentric acceleration
		f64_3 a_common(0);
		for (size_t i = 2; i < p.n_alive(); i++)    
		{
			float64_t mfac = p.m()[i] * this->planet_inverse_helio_cubed[i];
			a_common -= p.r()[i] * mfac;
		}

		// Load this into all the arrays
		for (size_t i = 1; i < p.n_alive(); i++)    
		{
			planet_a[i] = a_common;
		}

		planet_h0_log.get<true>()[index] = a_common - p.r()[1] * p.m()[1] * this->planet_inverse_helio_cubed[1];
		
		// Now do indirect acceleration ; note that planet 1 does not receive a contribution 
		for (size_t i = 2; i < p.n_alive(); i++)    
		{
			planet_a[i] += (this->planet_rj[i] * this->planet_inverse_jacobi_cubed[i] - p.r()[i] * this->planet_inverse_helio_cubed[i]) * p.m()[0];
		}
		
		/* next term ; again, first planet does not participate */
		f64_3 a_accum(0);
		for (size_t i = 2; i < p.n_alive(); i++)    
		{
			float64_t mfac = p.m()[i] * p.m()[0] * this->planet_inverse_jacobi_cubed[i] / this->planet_eta[i-1];
			a_accum += this->planet_rj[i] * mfac;
			planet_a[i] += a_accum;
		}

		/* Finally, incorporate the direct accelerations */
		for (size_t i = 1; i < p.n_alive() - 1; i++)    
		{
			for (size_t j = i + 1; j < p.n_alive(); j++)    
			{
				f64_3 dr = p.r()[j] - p.r()[i];
				float64_t r2 = dr.lensq();
				float64_t irij3 = 1. / (r2 * std::sqrt(r2));

				float64_t mfac = p.m()[i] * irij3;
				planet_a[j] -= dr * mfac;

				// acc. on i is just negative, with m[j] instead
				mfac = p.m()[j] * irij3;
				planet_a[i] += dr * mfac;
			}
		}
	}

	template<bool danby>
	bool WHIntegrator::drift_single(float64_t dt, float64_t mu, f64_3* r, f64_3* v)
	{
		float64_t dist, vsq, vdotr;
		dist = std::sqrt(r->lensq());
		vsq = v->lensq();
		vdotr = v->x * r->x + v->y * r->y + v->z * r->z;

		float64_t energy = vsq;
		energy *= 0.5;
		energy -= mu / dist;

		if (energy < 0)
		{
			float64_t a = -0.5 * mu / energy;
			float64_t n_ = std::sqrt(mu / (a * a * a));
			float64_t ecosEo = 1.0 - dist / a;
			float64_t esinEo = vdotr / (n_ * a * a);
			float64_t esq = ecosEo * ecosEo + esinEo * esinEo;

			// subtract off an integer multiple of complete orbits
			float64_t dM = dt * n_ - M_2PI * (int) (dt * n_ / M_2PI);

			// remaining time to advance
			dt = dM / n_;
			// maybe parallelize this

			if (danby)
			{
				if ((dM * dM > 0.16) || (esq > 0.36) || (esq * dM * dM > 0.0016)) goto hyperbolic;

				double dE, sindE, cosdE;
				kepmd(dM, esinEo, ecosEo, &dE, &sindE, &cosdE);

				double fchk = dE - ecosEo * sindE + esinEo * (1. - cosdE) - dM;
				if (fchk * fchk > DANBYB)
				{
					return true;
				}

				float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
				float64_t f = 1.0 + a * (cosdE - 1.0) / dist;
				float64_t g = dt + (sindE - dE) / n_;
				float64_t fdot = -n_ * sindE * a / (dist * fp);
				float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

				f64_3 r0 = *r;
				f64_3 v0 = *v;
				*r = r0 * f + v0 * g;
				*v = r0 * fdot + v0 * gdot;
				return false;
			}
			else
			{
				// call kepler equation solver with initial guess in dE already
				float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
				float64_t sindE, cosdE;

				uint32_t its;
				if (kepeq(dM, esinEo, ecosEo, &dE, &sindE, &cosdE, &its))
				{
					return true;
				}

				float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
				float64_t f = 1.0 + a * (cosdE - 1.0) / dist;
				float64_t g = dt + (sindE - dE) / n_;
				float64_t fdot = -n_ * sindE * a / (dist * fp);
				float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

				f64_3 r0 = *r;
				f64_3 v0 = *v;
				*r = r0 * f + v0 * g;
				*v = r0 * fdot + v0 * gdot;
				return false;
			}
		}
		
	hyperbolic:
		double fp, c1, c2, c3;
		if (kepu(dt, dist, mu, -2 * energy, vdotr, &fp, &c1, &c2, &c3))
		{
			return true;
		}

		double f = 1.0 - (mu / dist) * c2;
		double g = dt - mu * c3;
		double fdot = -(mu / (fp * dist)) * c1;
		double gdot = 1. - (mu / fp) * c2;

		f64_3 r0 = *r;
		f64_3 v0 = *v;
		*r = r0 * f + v0 * g;
		*v = r0 * fdot + v0 * gdot;
		return false;
	}


	const uint32_t HP_DRIFT_ITERATIONS = 10;
	bool WHIntegrator::drift_single_hp(float64_t dt, float64_t mu, f64_3* r, f64_3* v)
	{
		if (drift_single<true>(dt, mu, r, v))
		{
			dt /= static_cast<double>(HP_DRIFT_ITERATIONS);

			for (uint8_t i = 0; i < 10; i++)
			{
				if (drift_single<true>(dt, mu, r, v))
				{
					return true;
				}
			}
		}

		return false;
	}

	template<bool fixediterations>
	void WHIntegrator::drift(float64_t t, Vf64_3& r, Vf64_3& v, size_t start, size_t n, Vf64& dist, Vf64& energy, Vf64& vdotr, Vf64& mu, Vu8& mask)
	{
		for (size_t i = start; i < start + n; i++)
		{
			dist[i] = std::sqrt(r[i].lensq());
			energy[i] = v[i].lensq();
			vdotr[i] = v[i].x * r[i].x + v[i].y * r[i].y + v[i].z * r[i].z;
		}

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
				std::ostringstream ss;
				ss << "unbound orbit of planet " << i << " energy = " << energy[i] << std::endl;

				/*
				for (size_t j = start; j < start + n; j++)
				{
					ss << "p " << r[j].x << " " << r[j].y << " " << r[j].z << std::endl;
					ss << "v " << v[j].x << " " << v[j].y << " " << v[j].z << std::endl;
				}
				*/
				
				throw std::runtime_error(ss.str());
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
				// float64_t e = std::sqrt(ecosEo * ecosEo + esinEo * esinEo);

				// subtract off an integer multiple of complete orbits
				float64_t dM = t * n_ - M_2PI * (int) (t * n_ / M_2PI);

				// remaining time to advance
				float64_t adv_dt = dM / n_;

				// call kepler equation solver with initial guess in dE already
				float64_t dE = dM - esinEo + esinEo * std::cos(dM) + ecosEo * std::sin(dM);
				float64_t sindE, cosdE;

				bool error;
				if (fixediterations)
				{
					error = kepeq_fixed(dM, esinEo, ecosEo, &dE, &sindE, &cosdE, 5);
				}
				else
				{
					uint32_t its;
					error = kepeq(dM, esinEo, ecosEo, &dE, &sindE, &cosdE, &its);
				}

				if (error)
				{
					throw std::runtime_error("Unconverging kepler");
				}

				float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
				float64_t f = 1.0 + a * (cosdE - 1.0) / dist[i];
				float64_t g = adv_dt + (sindE - dE) / n_;
				float64_t fdot = -n_ * sindE * a / (dist[i] * fp);
				float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

				r[i] = r0 * f + v0 * g;
				v[i] = r0 * fdot + v0 * gdot;
			}
		}
	}

	void WHIntegrator::step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, float64_t dt, size_t timestep_index)
	{
		for (size_t i = begin; i < begin + length; i++)
		{
			this->particle_mask[i] = pa.deathflags()[i] != 0;

			if (!this->particle_mask[i])
			{
				pa.v()[i] += this->particle_a[i] * (dt / 2);
			}
		}

		// Drift all the particles along their Jacobi Kepler ellipses
		// Can change false to true to use fixed iterations
		drift<true>(dt, pa.r(), pa.v(), begin, length, particle_dist, particle_energy, particle_vdotr, particle_mu, particle_mask);

		// find the accelerations of the heliocentric velocities
		helio_acc_particles<false>(pl, pa, begin, length, t, timestep_index);

		for (size_t i = begin; i < begin + length; i++)
		{
			if (!this->particle_mask[i])
			{
				pa.v()[i] += particle_a[i] * (dt / 2);
				pa.deathtime_index()[pa.id()[i]] = static_cast<uint32_t>(timestep_index + 1);
			}
		}
	}

	void WHIntegrator::gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length)
	{
		gather(particle_a, indices, begin, length);
	}

	void WHIntegrator::load_h0(const HostPlanetPhaseSpace& pl)
	{
		const Vf64_3* r_log = &pl.r_log().old;

		for (size_t index = 0; index < pl.r_log().len_old; index++)
		{
			f64_3 h0(0);
			f64_3 a_common(0);
			for (size_t i = 1; i < pl.n_alive(); i++)
			{
				f64_3 r = (*r_log)[pl.log_index_at<true>(index, i)];
				float64_t r2 = r.lensq();
				double inv_r3 = 1. / (std::sqrt(r2) * r2);

				if (i == 1)
				{
					h0 = -r * pl.m()[1] * inv_r3;
				}

				if (i > 1)
				{
					float64_t mfac = pl.m()[i] * inv_r3;
					a_common -= r * mfac;
				}
			}

			h0 += a_common;
		
			planet_h0_log.old[index] = h0;
		}

		planet_h0_log.len_old = pl.r_log().len_old;
	}

	void WHIntegrator::step_planets(HostPlanetPhaseSpace& pl, float64_t t, double dt, size_t timestep_index)
	{
		(void) t;

		double new_dt = dt;

		for (size_t i = 1; i < pl.n_alive(); i++)
		{
			pl.v()[i] += planet_a[i] * (new_dt / 2);
		}

		// Convert the heliocentric velocities to Jacobi velocities 
		sr::convert::helio_to_jacobi_v_planets(pl, planet_eta, planet_vj);

		for (size_t i = 1; i < pl.n_alive(); i++)
		{
			// Each Jacobi Kepler problem has a different central mass
			this->planet_mu[i] = pl.m()[0] * this->planet_eta[i] / this->planet_eta[i - 1];
			this->planet_mask[i] = 0;
		}

		// Drift all the particles along their Jacobi Kepler ellipses
		drift<false>(new_dt, this->planet_rj, this->planet_vj, 1, pl.n_alive() - 1, planet_dist, planet_energy, planet_vdotr, planet_mu, planet_mask);

		// convert Jacobi vectors to helio. ones for acceleration calc 
		sr::convert::jacobi_to_helio_planets(planet_eta, planet_rj, planet_vj, pl);

		// find the accelerations of the heliocentric velocities
		helio_acc_planets(pl, timestep_index);

		Vf64_3* r_log = &pl.r_log().old;
		Vf64_3* v_log = &pl.v_log().old;

		std::copy(pl.r().begin() + 1, pl.r().end(), r_log->begin() + (pl.n_alive() - 1) * timestep_index);
		std::copy(pl.v().begin() + 1, pl.v().end(), v_log->begin() + (pl.n_alive() - 1) * timestep_index);

		for (size_t i = 1; i < pl.n_alive(); i++)
		{
			pl.v()[i] += planet_a[i] * (new_dt / 2);
		}
	}
}
}
