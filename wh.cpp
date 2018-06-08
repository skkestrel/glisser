#include "wh.h"
#include "convert.h"

#include <cmath>

void helio_acc_planets(HostPlanetPhaseSpace& p)
{
	Hvf64 inverse_helio_cubed(p.n), inverse_jacobi_cubed(p.n);

	for (size_t i = 1; i < p.n; i++)
	{
		f64 r2 = p.r.lensq();
		inverse_helio_cubed[i] = 1. / (std::sqrt(r2) * r2);
		f64 r2 = p.rj.lensq();
		inverse_jacobi_cubed[i] = 1. / (std::sqrt(r2) * r2);
        }
	
        /* compute common heliocentric acceleration */
	axs = ays = azs = 0.0;
	for(i=2; i<=npl; i++)    
	{
		mfac = m[i]*irh3[i];
		axs -= mfac*x[i];
		ays -= mfac*y[i];
		azs -= mfac*z[i];
        }

        /* Load this into all the arrays */
	for(i=1; i<=npl; i++)    
	{
		ax[i] = axs;
		ay[i] = ays;
		az[i] = azs;
        }
	
	/* Now do indirect acceleration ; note that planet 1 does not receive
           a contribution 
         */
	for(i=2; i<=npl; i++)    
	{
		ax[i] += m[0] * ( xj[i]*irj3[i] - x[i]*irh3[i] );
		ay[i] += m[0] * ( yj[i]*irj3[i] - y[i]*irh3[i] );
		az[i] += m[0] * ( zj[i]*irj3[i] - z[i]*irh3[i] );
        }
	
	/* next term ; again, first planet does not participate */
	axs = ays = azs = 0.0;
	for(i=2; i<=npl; i++)    
	{
		mfac = m[i]*m[0]*irj3[i]/eta[i-1];
		axs += mfac*xj[i];
		ays += mfac*yj[i];
		azs += mfac*zj[i];
		ax[i] += axs;
		ay[i] += ays;
		az[i] += azs;
        }

	/* Finally, incorporate the direct accelerations */
	for(i=1; i<=npl-1; i++)    
	{
	   for(j=i+1; j<=npl; j++)    
	   {
		dx = x[j] - x[i];
		dy = y[j] - y[i];
		dz = z[j] - z[i];
		r2 = dx*dx + dy*dy + dz*dz;
		irij3 = 1.0/ ( r2*sqrt(r2) );

		mfac = m[i]*irij3;
		ax[j] -= mfac*dx;
		ay[j] -= mfac*dy;
		az[j] -= mfac*dz;

		/* acc. on i is just negative, with m[j] instead */
		mfac = m[j]*irij3;
		ax[i] += mfac*dx;
		ay[i] += mfac*dy;
		az[i] += mfac*dz;
	   }
        }
}

void initialize(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa)
{

}

void step_planets(HostPlanetPhaseSpace& pl, double dt)
{
        /* 
	Kick velocities with non-Keplerian component of accelerations 
	for half a time step. This does not affect the Jacobi positions.
	The acceleration arrays should be already stored!
	*/

	for (size_t i = 1; i < pl.n; i++)
	{
		pl.v[i].x += pl.a[i].x * dt / 2;
		pl.v[i].y += pl.a[i].y * dt / 2;
		pl.v[i].z += pl.a[i].z * dt / 2;
	}

	// Convert the heliocentric velocities to Jacobi velocities 
	helio_to_jacobi_v(pl);

	// Drift all the particles along their Jacobi Kepler ellipses
	for (size_t i = 1; i < pl.n; i++)
	{
            // Each Jacobi Kepler problem has a different central mass

	    double mu = m[0] * pl.eta[i] / pl.eta[i - 1];
	    drift(mu,dt0,&xj[i],&yj[i],&zj[i],&vxj[i],&vyj[i],&vzj[i]);
        }

	// convert Jacobi vectors to helio. ones for acceleration calc 
	jacobi_to_helio(pl);

	// find the accelerations of the heliocentric velocities
	helio_interaction_planets(pl);

	for (size_t i=1; i < pl.n; i++)
	{
		pl.v[i].x += pl.a[i].x * dt / 2;
		pl.v[i].y += pl.a[i].y * dt / 2;
		pl.v[i].z += pl.a[i].z * dt / 2;
	}
}
