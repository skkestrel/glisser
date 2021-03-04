#include "convert.h"
#include <cmath>
#include <iostream>

namespace sr
{
namespace convert
{
	void jacobi_to_helio_planets(const Vf64& eta, const Vf64_3& rj, const Vf64_3& vj, HostPlanetPhaseSpace& pl)
	{
		// sun at origin
		pl.r()[0] = f64_3(0);
		pl.v()[0] = f64_3(0);

		pl.r()[1] = rj[1];
		pl.v()[1] = vj[1];

		f64_3 rsum(0), vsum(0);

		for (size_t i = 2; i < pl.n(); i++)
		{
			rsum += rj[i-1] * (pl.m()[i - 1] / eta[i-1]);
			vsum += vj[i-1] * (pl.m()[i - 1] / eta[i-1]);

			pl.r()[i] = rj[i] + rsum;
			pl.v()[i] = vj[i] + vsum;
		}
	}

	/*
	void helio_to_jacobi_r_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p)
	{
		for (size_t i = 0; i < p.n(); i++)
		{
			p.rj[i] = p.r()[i] - pl.bary_r;
		}
	}

	void helio_to_jacobi_v_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p)
	{
		for (size_t i = 0; i < p.n(); i++)
		{
			p.vj[i] = p.v()[i] - pl.bary_v;
		}
	}
	*/

	void helio_to_jacobi_v_planets(const HostPlanetPhaseSpace& p, const Vf64& eta, Vf64_3& vj)
	{
		// COM
		vj[0] = f64_3(0);

		// same as heliocentric
		vj[1] = p.v()[1];

		// momentum sum
		f64_3 psum(0);

		for (size_t i = 2; i < p.n(); i++)
		{
			// velocity of interior COM
			f64_3 vsum;

			psum += p.v()[i - 1] * p.m()[i - 1];
			vsum = psum / eta[i - 1];

			vj[i] = p.v()[i] - vsum;
		}

		psum += p.v()[p.n() - 1] * p.m()[p.n() - 1];
		// p.bary_v = psum / eta[p.n() - 1];
	}

	void helio_to_jacobi_r_planets(const HostPlanetPhaseSpace& p, const Vf64& eta, Vf64_3& rj)
	{
		// pick origin at baricenter
		rj[0] = f64_3(0);
		// first jacobi coordinate is same as heliocentric
		rj[1] = p.r()[1];

		f64_3 sum(0);
		f64_3 bary;

		for (size_t i = 2; i < p.n(); i++)
		{
			sum += p.r()[i - 1] * p.m()[i - 1];
			bary = sum / eta[i - 1];

			rj[i] = p.r()[i] - bary;
		}

		sum += p.r()[p.n() - 1] * p.m()[p.n() - 1];
		// p.bary_r = sum / p.eta[p.n() - 1];
	}

	void find_barycenter(const Vf64_3& r, const Vf64_3& v, const Vf64& m, size_t n, f64_3& r_out, f64_3& v_out, double& m_out)
	{
		f64_3 rsum(0);
		f64_3 vsum(0);

		m_out = 0;

		for (size_t i = 0; i < n; i++)
		{
			rsum += r[i] * m[i];
			vsum += v[i] * m[i];

			m_out += m[i];
		}

		r_out = rsum / m_out;
		v_out = vsum / m_out;
	}

	void to_bary(HostData& hd)
	{
		f64_3 r, v;
		double m;
		find_barycenter(hd.planets.r(), hd.planets.v(), hd.planets.m(), hd.planets.n(), r, v, m);

		for (size_t i = 0; i < hd.planets.n(); i++)
		{
			hd.planets.r()[i] -= r;
			hd.planets.v()[i] -= v;
		}

		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			hd.particles.r()[i] -= r;
			hd.particles.v()[i] -= v;
		}
	}

	void to_helio(HostData& hd)
	{
		f64_3 r;
		f64_3 v;
		r = hd.planets.r()[0];
		v = hd.planets.v()[0];

		for (size_t i = 0; i < hd.planets.n(); i++)
		{
			hd.planets.r()[i] -= r;
			hd.planets.v()[i] -= v;
		}

		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			hd.particles.r()[i] -= r;
			hd.particles.v()[i] -= v;
		}
	}

	// ehie - see orbel/orbel_ehie.f in swift
	double ehie(double e, double m)
	{
		int iflag = 0;
		int nper = static_cast<int>(m / (2 * M_PI));
		m = m - nper * 2 * M_PI;
		if (m < 0) m += 2 * M_PI;
		if (m > M_PI)
		{
			m = 2 * M_PI - m;
			iflag = 1;
		}
		
		// make a first guess
		double x = std::pow(6 * m, 1. / 3) - m;

		// NMAX = 3
		for (size_t i = 0; i < 3; i++)
		{
			double sa = std::sin(x + m);
			double ca = std::cos(x + m);
			double esa = e * sa;
			double eca = e * ca;
			double f = x - esa;
			double fp = 1 - eca;
			double dx = -f / fp;
			dx = -f / (fp + dx * esa / 2);
			dx = -f / (fp + dx * (esa + eca * dx / 3) / 2);
			x = x + dx;
		}

		if (iflag)
		{
			return 2 * M_PI - m - x;
			m = 2 * M_PI - m;
		}
		else
		{
			return m + x;
		}
	}

	template<int n>
	double eget(double e, double m)
	{
		double sm = std::sin(m);
		double cm = std::cos(m);
		double x = m + e * sm * (1 + e * (cm + e * (1 - 1.5 * sm * sm)));

		for (size_t i = 0; i < n; i++)
		{
			double sx = std::sin(x);
			double cx = std::cos(x);
			double es = e*sx;
			double ec = e*cx;
			double f = x - es - m;
			double fp = 1 - ec;
			double fpp = es;
			double fppp = ec;
			double dx = -f / fp;
			dx = -f / (fp + dx * fpp / 2);
			dx = -f / (fp + dx * fpp / 2 + dx * dx * fppp / 6);

			x = x + dx;
		}
		return x;
	}
	
	double ehybrid(double e, double m)
	{
		if (e < 0.18) return eget<1>(e, m);
		else if (e < 0.8) return eget<2>(e, m);
		else return ehie(e, m);
	}

	void from_elements_M(double mu, double a, double e, double i, double capom, double om, double M, f64_3* r_, f64_3* v)
	{
		double E = ehybrid(e, M);
		double f = std::acos((std::cos(E) - e) / (1 - e * std::cos(E)));

		E = E - static_cast<int>(E / 2 / M_PI) * 2 * M_PI;
		if (E > M_PI) E -= 2 * M_PI;
		else if (E < -M_PI) E += 2 * M_PI;

		// fix the sign of f
		if (E < 0) f = -f;

		/*
		std::cout << "M = " << M << " f = " << f << std::endl;
		std::cout << "E - e sin E = " << E - e * std::sin(E) << std::endl;
		*/

		from_elements(mu, a, e, i, capom, om, f, r_, v);
	}

	void from_elements(double mu, double a, double e, double i, double capom, double om, double f, f64_3* r_, f64_3* v)
	{
		using namespace std;

		double& x = r_->x;
		double& y = r_->y;
		double& z = r_->z;
		double& vx = v->x;
		double& vy = v->y;
		double& vz = v->z;

		double prec;
		double u, r, xhat, yhat, zhat;
		double h, hx, hy, hz;
		double thx, thy, thz;		/* components of theta hat */
		double thdot,  rdot;

		/* user must set this */
		prec = 1.0e-13;

		/* compute the r unit vector */
		u = om + f;
		xhat = cos(u)*cos(capom) - cos(i)*sin(capom)*sin(u);
		yhat = cos(u)*sin(capom) + cos(i)*cos(capom)*sin(u);
		zhat = sin(i)*sin(u);

		/* compute components of the specific angular momentum vector. Note
		   this is a unit vector, with magnitude in the variable h (below) 
		*/
		hx = sin(capom)*sin(i);
		hy = -cos(capom)*sin(i);
		hz = cos(i);

		/* compute the radius vector and magnitude of h */
		if( fabs( e - 1.0 ) > prec) {
		  r = a * (1.0 - e*e) / (1.0 + e*cos(f));
		  h = sqrt( mu*a*(1.0 - e*e) );
		} else {
		  h = sqrt( 2.0*mu*a );		/* a holds q for parabolic */  
		  r = a/( cos(f/2.0)*cos(f/2.0) );
		}
		/* This provides the position vector */
		x = r * xhat;
		y = r * yhat;
		z = r * zhat;

		/* compute the components of the unit vector theta hat */
		thx = hy * zhat - hz * yhat;
		thy = hz * xhat - hx * zhat;
		thz = hx * yhat - hy * xhat;

		/* obtain the velocity vector's components and calculate v */
		thdot =  h/(r*r);
		rdot  =  e*mu*sin(f)/h;

		vx = r * thdot * thx + rdot * xhat; 
		vy = r * thdot * thy + rdot * yhat; 
		vz = r * thdot * thz + rdot * zhat; 
	}

	void to_elements(double mu, f64_3 r, f64_3 v, int* esignout, double* aout, double* eout, double* iout, double* capomout, double* omout, double* fout)
	{
		using namespace std;

		double a, e, i, capom, om, f;
		int esign;

		double x = r.x, y = r.y, z = r.z;
		double vx = v.x, vy = v.y, vz = v.z;
		double pi,prec;
		double hsq,hx,hy,hz;
		double rr,xhat,yhat,zhat;
		double vsq,vdotr,fac;
		double Px,Py,Pz,modP,nx,ny,ecosw;
		double energy;

		pi = 2.0*asin(1.0);
		/* machine precision, user must set this */
		prec = 1.0e-13;

		/* compute the specific angular momentum */
		hx = y*vz - z*vy;
		hy = z*vx - x*vz;
		hz = x*vy - y*vx;
		hsq = hx*hx + hy*hy + hz*hz;

		/* As long as we are not on a radial orbit, compute elements */
		if ( hsq > prec )
		{

			/* compute the orbital inclination */
			i = acos( hz/sqrt(hsq) );

			/* compute the longitude of the ascending node */
			if( fabs( i ) < prec) {
				capom = 0.0;
			} 
			else if( fabs( pi - fabs( i) ) < prec ) {
				capom = 0.0;
			} 
			else {
				capom = atan2(hx, -hy);
			}

			/* compute some required quantities */
			vsq = vx*vx + vy*vy + vz*vz;
			vdotr =  x*vx + y*vy + z*vz;
			rr = sqrt(x*x + y*y + z*z);
			xhat = x/rr;
			yhat = y/rr;
			zhat = z/rr;
			nx = cos( capom);
			ny = sin( capom);

			/* compute the Hamilton vector and thus the eccentricity */
			fac = vsq * rr - mu;
			Px = fac * xhat - vdotr * vx;
			Py = fac * yhat - vdotr * vy;
			Pz = fac * zhat - vdotr * vz;
			modP = sqrt( Px*Px + Py*Py + Pz*Pz );
			e = modP / mu;

			/* compute the argument of pericenter */
			if (fabs(e) < prec) {
				om = 0.0;
			} 
			else {
				if ((i < prec) || (pi - i < prec)) {
					om = atan2(Py, Px);
				} else {
					ecosw = (nx * Px + ny * Py) / mu;
					om = acos(ecosw / e);
					if ( fabs(Pz) > prec ) {
						/* resolve sign ambiguity by sign of Pz  */
						om *= fabs(Pz)/Pz;
					}
				}
			}

			/* compute the orbital energy , and depending on its sign compute
			   the semimajor axis (or pericenter) and true anomaly      */
			energy = vsq/2.0 - mu/rr;
			if( fabs(energy) < prec) {
				esign = 0;		/* parabolic */
				a = 0.5 * hsq / mu;	/* actually PERICENTRIC DISTANCE */
				if ( fabs(vdotr) < prec ) {
					f = 0.0;
				} else {
					f = 2.0 * acos(sqrt(a/rr)) * vdotr/fabs(vdotr);
				} 
			} else if (energy > 0.0) {
				esign = 1;		/* hyperbolic */
				a = -0.5 * mu/energy;  /* will be negative */

				if ( fabs(vdotr) < prec ) {
					f = 0.0;
				} else {
					fac =  a * (1.0 - e * e)/rr - 1.0;
					f =  acos( fac/ e ) * vdotr/fabs(vdotr);
				} 
			} else {
				esign = -1;		/* elliptic */
				a = -0.5 * mu/energy;
				if ( fabs( e ) > prec ) {      
					if ( fabs(vdotr) < prec ) {
						if ( rr < a ) {		/* determine apside */
							f = 0.0;
						} else {
							f = pi;
						}
					} else {
						fac =  a * (1.0 - e * e)/rr - 1.0;
						f =  acos( fac/ e ) * vdotr/fabs(vdotr);

						// KZ
						if (std::isnan(f)) {
							if ( rr < a ) {		/* determine apside */
								f = 0.0;
							} else {
								f = pi;
							}
						}
					} 
				} else {                       /* compute circular cases */
					fac = (x * nx + y * ny)/rr;
					f = acos(fac);
					if ( fabs(z) > prec ) {
						/* resolve sign ambiguity by sign of z  */
						f *= fabs(z)/z;
					} else if ( (i < prec) || (pi - i < prec) ) {
						f = atan2(y,x) * cos( i);
					} 

				}
			}

		} else { 				/* PANIC: radial orbit */
			esign = 1;			/* call it hyperbolic */
			a = sqrt(x*x + y*y + z*z);
			e = HUGE_VAL;
			i = asin(z/sqrt(x*x + y*y + z*z) );	/* latitude above plane */
			capom = atan2(y,x);			/* azimuth */
			om = HUGE_VAL;
			f = HUGE_VAL;
		}

		if (esignout) *esignout = esign;
		if (aout) *aout = a;
		if (eout) *eout = e;
		if (iout) *iout = i;
		if (capomout) *capomout = capom;
		if (omout) *omout = om;
		if (fout) *fout = f;
	}

	double get_mean_anomaly(double e, double f)
	{
		double E = std::acos((e + std::cos(f)) / (1 + e * std::cos(f)));
		if (f < 0) E = -E;
		return E - e * std::sin(E);
	}

	double get_bary_mu(double center_mass, double pl_mass)
	{
		return std::pow(center_mass - pl_mass, 3) / std::pow(center_mass, 2);
	}
}
}
