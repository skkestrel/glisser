#include "data.h"
#include <fstream>
#include <iomanip>
#include <limits>

bool load_data(HostData& hd, std::string plin, std::string icsin, size_t tbsize, size_t max_particle, bool readmomenta)
{
	hd.tbsize = tbsize;
	std::ifstream plinfile(plin), icsinfile(icsin);

	size_t npl;
	plinfile >> npl;

	hd.planets = HostPlanetPhaseSpace(npl, tbsize);

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
			hd.particles.flags[i] = 0;
		}
		else
		{
			hd.particles.deathtime[i] = std::stod(s);
			icsinfile >> hd.particles.flags[i] >> hd.particles.id[i];
		}
	}

	return false;
}

void save_data(const HostData& hd, const DeviceData& dd, std::string plout, std::string icsout, std::string infoout)
{
	std::ofstream ploutfile(plout), icsoutfile(icsout), infooutfile(infoout);

	ploutfile << hd.planets.n << std::endl;
	ploutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.planets.n; i++)
	{
		ploutfile << hd.planets.m[i] << std::endl;
		ploutfile << hd.planets.r[i].x << " " << hd.planets.r[i].y << " " << hd.planets.r[i].z << std::endl;
		ploutfile << hd.planets.v[i].x << " " << hd.planets.v[i].y << " " << hd.planets.v[i].z << std::endl;
	}

	icsoutfile << hd.particles.n << std::endl;
	icsoutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		icsoutfile << hd.particles.r[i].x << " " << hd.particles.r[i].y << " " << hd.particles.r[i].z << std::endl;
		icsoutfile << hd.particles.v[i].x << " " << hd.particles.v[i].y << " " << hd.particles.v[i].z << std::endl;
		icsoutfile << hd.particles.deathtime[i] << " " << hd.particles.flags[i] << " " << hd.particles.id[i] << std::endl;
	}

	infooutfile << std::setprecision(17);
	infooutfile << hd.t << " " << hd.particles.n_alive << std::endl;
}

void transfer_data(const HostData& hd, DeviceData& dd)
{
	dd.particles0 = DeviceParticlePhaseSpace(hd.particles.n);
	dd.particles1 = DeviceParticlePhaseSpace(hd.particles.n);

	dd.planets0 = DevicePlanetPhaseSpace(hd.planets.n, hd.tbsize);
	dd.planets1 = DevicePlanetPhaseSpace(hd.planets.n, hd.tbsize);

	/*
	dd.log_buffer_id = 0;
	dd.phase_space_id = 0;

	dd.planets.r_log0 = Dvf64_3(hd.planets.r_log.size());
	dd.planets.r_log1 = Dvf64_3(hd.planets.r_log.size());
	dd.planets.m = Dvf64(hd.n_planet);

	dd.ps0.r = dd.ps0.v = Dvf64_3(hd.n_part);
	dd.ps0.flags = Dvu8(hd.n_part);
	dd.ps0.deathtime = Dvu32(hd.n_part);
	dd.ps0.id = Dvu32(hd.n_part);

	dd.ps1.r = dd.ps1.v = Dvf64_3(hd.n_part);
	dd.ps1.flags = Dvu8(hd.n_part);
	dd.ps1.deathtime = Dvu32(hd.n_part);
	dd.ps1.id = Dvu32(hd.n_part);

	dd.gather_indices = Dvu32(hd.n_part);

	dd.particles.n_alive = 0;
	for (size_t i = 0; i < hd.n_part; i++)
	{
		if (~hd.particles.flags[i] & 0x0001)
		{
			dd.particles.n_alive++;
		}
	}

	dd.coefdt = Dvf64(hd.coefdt.size());

	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;
	thrust::copy(hd.particles.r.begin(), hd.particles.r.end(), ps.r.begin());
	thrust::copy(hd.particles.v.begin(), hd.particles.v.end(), ps.v.begin());

	thrust::copy(hd.planets.m.begin(), hd.planets.m.end(), dd.planets.m.begin());
	thrust::copy(hd.coefdt.begin(), hd.coefdt.end(), dd.coefdt.begin());

	thrust::copy(hd.id.begin(), hd.id.end(), ps.id.begin());
	thrust::copy(hd.deathtime.begin(), hd.deathtime.end(), ps.deathtime.begin());
	thrust::copy(hd.particles.flags.begin(), hd.particles.flags.end(), ps.flags.begin());
	*/
}

void recover_data(HostData& hd, const DeviceData& dd, cudaStream_t& stream)
{
	/*
	const DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;

	size_t n = hd.n_part;


	thrust::copy(thrust::cuda::par.on(stream), ps.r.begin(), ps.r.begin() + n, hd.particles.r.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.v.begin(), ps.v.begin() + n, hd.particles.v.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.flags.begin(), ps.flags.begin() + n, hd.particles.flags.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.deathtime.begin(), ps.deathtime.begin() + n, hd.particles.deathtime.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.id.begin(), ps.id.begin() + n, hd.particles.id.begin());
	*/

	/*
	// zip will crash the program

	auto iterator = thrust::make_zip_iterator(thrust::make_tuple(
				ps.r.begin(),
				ps.v.begin(),
				ps.flags.begin(), ps.deathtime.begin(), ps.id.begin()));
	thrust::copy(thrust::cuda::par.on(stream),
			iterator,
			iterator + n,
			thrust::make_zip_iterator(thrust::make_tuple(
					hd.particles.r.begin(), hd.particles.v.begin(),
					hd.particles.flags.begin(), hd.deathtime.begin(), hd.id.begin())));
	 */
}

