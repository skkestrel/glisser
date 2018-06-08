#include "convert.h"

// requires pl to have bary_r calculated
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

void helio_to_jacobi_r_planets(HostPlanetPhaseSpace& p)
{
	p.eta[0] = p.m[0];

	for (size_t i = 1; i < p.n; i++)
	{
		p.eta[i] = p.eta[i - 1] + p.m[i];
	}

	// compute jacobi coordinates

	// pick origin at baricenter
	p.rj[0].x = p.rj[0].y = p.rj[0].z = 0.0;

	// first jacobi coordinate is same as heliocentric
	p.rj[1] = p.r[1];

	f64_3 sum;
	f64_3 bary;

	sum.x = 0;
	sum.y = 0;
	sum.z = 0;

	for (size_t i = 2; i < p.n; i++)
	{
		sum.x += p.m[i - 1] * p.r[i - 1].x;
		sum.y += p.m[i - 1] * p.r[i - 1].y;
		sum.z += p.m[i - 1] * p.r[i - 1].z;

		bary.x = sum.x / eta[i - 1];
		bary.y = sum.y / eta[i - 1];
		bary.z = sum.z / eta[i - 1];

		p.xj[i] = p.x[i] - bary.x;
		p.yj[i] = p.y[i] - bary.y;
		p.zj[i] = p.z[i] - bary.z;
	}

	sum.x += p.m[p.n - 1] * p.r[p.n - 1].x;
	sum.y += p.m[p.n - 1] * p.r[p.n - 1].y;
	sum.z += p.m[p.n - 1] * p.r[p.n - 1].z;

	p.bary_r.x = sum.x / eta[p.n - 1];
	p.bary_r.y = sum.y / eta[p.n - 1];
	p.bary_r.z = sum.z / eta[p.n - 1];
}

void helio_to_jacobi_v_planets(HostPlanetPhaseSpace& p)
{
	// momentum sum
	f64_3 psum;

	// COM
	p.vj.x[0] = p.vj.y[0] = p.vj.z[0] = 0;

	// same as heliocentric
	p.vj[1] = p.v[1];

	psum.x = psum.y = psum.z = 0;

	for (size_t i = 2; i < p.n; i++)
	{
		// velocity of interior COM
		f64_3 vsum;

		psum.x += p.m[i - 1] * p.v[i - 1].x;
		psum.y += p.m[i - 1] * p.v[i - 1].y;
		psum.z += p.m[i - 1] * p.v[i - 1].z;

		vsum.x = psum.x / p.eta[i - 1];
		vsum.y = psum.y / p.eta[i - 1];
		vsum.z = psum.z / p.eta[i - 1];

		p.vj[i].x = p.v[i].x - vsum.x;
		p.vj[i].y = p.v[i].y - vsum.y;
		p.vj[i].z = p.v[i].z - vsum.z;
	}

	psum.x += p.m[p.n - 1] * p.v[p.n - 1].x;
	psum.y += p.m[p.n - 1] * p.v[p.n - 1].y;
	psum.z += p.m[p.n - 1] * p.v[p.n - 1].z;

	p.bary_v.x = psum.x / p.eta[i - 1];
	p.bary_v.y = psum.y / p.eta[i - 1];
	p.bary_v.z = psum.z / p.eta[i - 1];
}

void helio_to_jacobi_r(HostPlanetPhaseSpace& p)
{
	p.eta[0] = p.m[0];

	for (size_t i = 1; i < p.n; i++)
	{
		p.eta[i] = p.eta[i - 1] + p.m[i];
	}

	// compute jacobi coordinates

	// pick origin at baricenter
	p.rj[0].x = p.rj[0].y = p.rj[0].z = 0.0;

	// first jacobi coordinate is same as heliocentric
	p.rj[1] = p.r[1];

	f64_3 sum;
	f64_3 bary;

	sum.x = 0;
	sum.y = 0;
	sum.z = 0;

	for (size_t i = 2; i < p.n; i++)
	{
		sum.x += p.m[i - 1] * p.r[i - 1].x;
		sum.y += p.m[i - 1] * p.r[i - 1].y;
		sum.z += p.m[i - 1] * p.r[i - 1].z;

		bary.x = sum.x / eta[i - 1];
		bary.y = sum.y / eta[i - 1];
		bary.z = sum.z / eta[i - 1];

		p.xj[i] = p.x[i] - bary.x;
		p.yj[i] = p.y[i] - bary.y;
		p.zj[i] = p.z[i] - bary.z;
	}

	sum.x += p.m[p.n - 1] * p.r[p.n - 1].x;
	sum.y += p.m[p.n - 1] * p.r[p.n - 1].y;
	sum.z += p.m[p.n - 1] * p.r[p.n - 1].z;

	p.bary_r.x = sum.x / eta[p.n - 1];
	p.bary_r.y = sum.y / eta[p.n - 1];
	p.bary_r.z = sum.z / eta[p.n - 1];
}

void helio_to_bary(HostData& hd)
{
	f64_3 r;
	f64_3 v;

	r.x = r.y = r.z = 0;
	v.x = v.y = v.z = 0;
	double totalm = 0;

	for (size_t i = 0; i < hd.n_planet; i++)
	{
		r.x += hd.planet.r[i].x * hd.planet.m[i];
		r.y += hd.planet.r[i].y * hd.planet.m[i];
		r.z += hd.planet.r[i].z * hd.planet.m[i];

		v.x += hd.planet.v[i].x * hd.planet.m[i];
		v.y += hd.planet.v[i].y * hd.planet.m[i];
		v.z += hd.planet.v[i].z * hd.planet.m[i];

		totalm += hd.planet.m[i];
	}

	r.x /= totalm;
	r.y /= totalm;
	r.z /= totalm;
	v.x /= totalm;
	v.y /= totalm;
	v.z /= totalm;

	for (size_t i = 0; i < hd.n_planet; i++)
	{
		hd.planet.r[i].x -= r.x;
		hd.planet.r[i].y -= r.y;
		hd.planet.r[i].z -= r.z;

		hd.planet.v[i].x -= v.x;
		hd.planet.v[i].y -= v.y;
		hd.planet.v[i].z -= v.z;
	}

	for (size_t i = 0; i < hd.n_part; i++)
	{
		hd.particle.r[i].x -= r.x;
		hd.particle.r[i].y -= r.y;
		hd.particle.r[i].z -= r.z;

		hd.particle.v[i].x -= v.x;
		hd.particle.v[i].y -= v.y;
		hd.particle.v[i].z -= v.z;
	}
}

void bary_to_helio(HostData& hd)
{
	f64_3 r;
	f64_3 v;
	r = hd.planet.r[0];
	v = hd.planet.v[0];

	for (size_t i = 0; i < hd.n_planet; i++)
	{
		hd.planet.r[i].x -= r.x;
		hd.planet.r[i].y -= r.y;
		hd.planet.r[i].z -= r.z;

		hd.planet.v[i].x -= v.x;
		hd.planet.v[i].y -= v.y;
		hd.planet.v[i].z -= v.z;
	}

	hd.planet.bary_r.x = -r.x;
	hd.planet.bary_r.y = -r.y;
	hd.planet.bary_r.z = -r.z;

	hd.planet.bary_v.x = -v.x;
	hd.planet.bary_v.y = -v.y;
	hd.planet.bary_v.z = -v.z;

	for (size_t i = 0; i < hd.n_part; i++)
	{
		hd.particle.r[i].x -= r.x;
		hd.particle.r[i].y -= r.y;
		hd.particle.r[i].z -= r.z;

		hd.particle.v[i].x -= v.x;
		hd.particle.v[i].y -= v.y;
		hd.particle.v[i].z -= v.z;
	}
}
