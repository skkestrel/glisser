#include "data.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>

template<typename T>
void gather(std::vector<T>& values, const std::vector<size_t>& indices)
{
	std::vector<T> copy(values.begin(), values.end());
	for (size_t i = 0; i < values.size(); i++)
	{
		values[i] = copy[indices[i]];
	}
}

void HostParticlePhaseSpace::stable_partition_alive()
{
	std::vector<size_t> indices(n);
	std::iota(indices.begin(), indices.end(), 0);
	n_alive = std::stable_partition(indices.begin(), indices.end(), [this](size_t index)
			{ return deathflags[index] == 0; }) - indices.begin();

	gather(r, indices);
	gather(v, indices);
	gather(deathtime, indices);
	gather(deathflags, indices);
	gather(id, indices);
}

bool load_data(HostData& hd, std::string plin, std::string icsin, size_t tbsize, size_t ce_factor, size_t max_particle, bool readmomenta)
{
	std::ifstream plinfile(plin), icsinfile(icsin);

	size_t npl;
	plinfile >> npl;

	hd.planets = HostPlanetPhaseSpace(npl, tbsize, ce_factor);

	for (size_t i = 0; i < npl; i++)
	{
		plinfile >> hd.planets.m[i];
		plinfile >> hd.planets.r[i].x >> hd.planets.r[i].y >> hd.planets.r[i].z;
		plinfile >> hd.planets.v[i].x >> hd.planets.v[i].y >> hd.planets.v[i].z;

		if (readmomenta)
		{
			hd.planets.v[i].x /= hd.planets.m[i];
			hd.planets.v[i].y /= hd.planets.m[i];
			hd.planets.v[i].z /= hd.planets.m[i];
		}
	}

	size_t npart;
	icsinfile >> npart;
	if (max_particle > 0) npart = std::min(npart, max_particle);

	hd.particles = HostParticlePhaseSpace(npart);

	for (size_t i = 0; i < npart; i++)
	{
		icsinfile >> hd.particles.r[i].x >> hd.particles.r[i].y >> hd.particles.r[i].z;
		icsinfile >> hd.particles.v[i].x >> hd.particles.v[i].y >> hd.particles.v[i].z;

		std::string s;
		icsinfile >> s;
		if (!isdigit(s[0]))
		{
			icsinfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			hd.particles.deathtime[i] = 0;
			hd.particles.id[i] = i;
			hd.particles.deathflags[i] = 0;
		}
		else
		{
			hd.particles.deathtime[i] = std::stod(s);
			icsinfile >> hd.particles.deathflags[i] >> hd.particles.id[i];
		}
	}

	hd.particles.stable_partition_alive();

	return false;
}

void save_data(const HostData& hd, std::string plout, std::string icsout)
{
	std::ofstream ploutfile(plout), icsoutfile(icsout);

	ploutfile << hd.planets_snapshot.n << std::endl;
	ploutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.planets_snapshot.n; i++)
	{
		ploutfile << hd.planets_snapshot.m[i] << std::endl;
		ploutfile << hd.planets_snapshot.r[i].x << " " << hd.planets_snapshot.r[i].y << " " << hd.planets_snapshot.r[i].z << std::endl;
		ploutfile << hd.planets_snapshot.v[i].x << " " << hd.planets_snapshot.v[i].y << " " << hd.planets_snapshot.v[i].z << std::endl;
	}

	icsoutfile << hd.particles.n << std::endl;
	icsoutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		icsoutfile << hd.particles.r[i].x << " " << hd.particles.r[i].y << " " << hd.particles.r[i].z << std::endl;
		icsoutfile << hd.particles.v[i].x << " " << hd.particles.v[i].y << " " << hd.particles.v[i].z << std::endl;
		icsoutfile << hd.particles.deathtime[i] << " " << hd.particles.deathflags[i] << " " << hd.particles.id[i] << std::endl;
	}
}
