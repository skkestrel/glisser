#include "convert.h"
#include <cmath>

void jacobi_to_helio_planets(const Vf64& eta, const Vf64_3& rj, const Vf64_3& vj, HostPlanetPhaseSpace& pl)
{
	// sun at origin
	pl.r[0] = f64_3(0);
	pl.v[0] = f64_3(0);

	pl.r[1] = rj[1];
	pl.v[1] = vj[1];

	f64_3 rsum(0), vsum(0);

	for (size_t i = 2; i < pl.n; i++)
	{
		rsum += rj[i-1] * (pl.m[i - 1] / eta[i-1]);
		vsum += vj[i-1] * (pl.m[i - 1] / eta[i-1]);

		pl.r[i] = rj[i] + rsum;
		pl.v[i] = vj[i] + vsum;
	}
}

/*
void helio_to_jacobi_r_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p)
{
	for (size_t i = 0; i < p.n; i++)
	{
		p.rj[i] = p.r[i] - pl.bary_r;
	}
}

void helio_to_jacobi_v_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p)
{
	for (size_t i = 0; i < p.n; i++)
	{
		p.vj[i] = p.v[i] - pl.bary_v;
	}
}
*/

void helio_to_jacobi_v_planets(const HostPlanetPhaseSpace& p, const Vf64& eta, Vf64_3& vj)
{
	// COM
	vj[0] = f64_3(0);

	// same as heliocentric
	vj[1] = p.v[1];

	// momentum sum
	f64_3 psum(0);

	for (size_t i = 2; i < p.n; i++)
	{
		// velocity of interior COM
		f64_3 vsum;

		psum += p.v[i - 1] * p.m[i - 1];
		vsum = psum / eta[i - 1];

		vj[i] = p.v[i] - vsum;
	}

	psum += p.v[p.n - 1] * p.m[p.n - 1];
	// p.bary_v = psum / eta[p.n - 1];
}

void helio_to_jacobi_r_planets(const HostPlanetPhaseSpace& p, const Vf64& eta, Vf64_3& rj)
{
	// pick origin at baricenter
	rj[0] = f64_3(0);

	// first jacobi coordinate is same as heliocentric
	rj[1] = p.r[1];

	f64_3 sum(0);
	f64_3 bary;

	for (size_t i = 2; i < p.n; i++)
	{
		sum += p.r[i - 1] * p.m[i - 1];
		bary = sum / eta[i - 1];

		rj[i] = p.r[i] - bary;
	}

	sum += p.r[p.n - 1] * p.m[p.n - 1];
	// p.bary_r = sum / p.eta[p.n - 1];
}

void find_barycenter(const Vf64_3& r, const Vf64_3& v, const Vf64& m, size_t n, f64_3& r_out, f64_3& v_out)
{
	f64_3 rsum(0);
	f64_3 vsum(0);

	double totalm = 0;

	for (size_t i = 0; i < n; i++)
	{
		rsum += r[i] * m[i];
		vsum += v[i] * m[i];

		totalm += m[i];
	}

	r_out = rsum / totalm;
	v_out = vsum / totalm;
}

void to_bary(HostData& hd)
{
	f64_3 r, v;
	find_barycenter(hd.planets.r, hd.planets.v, hd.planets.m, hd.planets.n, r, v);

	for (size_t i = 0; i < hd.planets.n; i++)
	{
		hd.planets.r[i] -= r;
		hd.planets.v[i] -= v;
	}

	for (size_t i = 0; i < hd.particles.n; i++)
	{
		hd.particles.r[i] -= r;
		hd.particles.v[i] -= v;
	}
}

void to_helio(HostData& hd)
{
	f64_3 r;
	f64_3 v;
	r = hd.planets.r[0];
	v = hd.planets.v[0];

	for (size_t i = 0; i < hd.planets.n; i++)
	{
		hd.planets.r[i] -= r;
		hd.planets.v[i] -= v;
	}

	for (size_t i = 0; i < hd.particles.n; i++)
	{
		hd.particles.r[i] -= r;
		hd.particles.v[i] -= v;
	}
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
