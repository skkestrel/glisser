#include "data.h"

bool load_data(HostData& hd, std::string plin, std::string icsin, size_t max_particle = 0, bool readmomenta = true)
{
	std::ifstream plinfile(plin), icsinfile(icsin);

	plinfile >> hd.n_planet;

	hd.r_planet = Hvf64_3(hd.n_planet);
	hd.v_planet = Hvf64_3(hd.n_planet);
	hd.m_planet = Hvf64(hd.n_planet);

	size_t lsize = hd.n_planet * (SIMPL_DEGREE - 1) * TIMEBLOCK_SIZE;
	hd.r_planet_log = Hvf64_3(lsize);

	for (size_t i = 0; i < hd.n_planet; i++)
	{
		plinfile >> hd.m_planet[i];
		plinfile >> hd.r_planet[i].x >> hd.r_planet[i].y >> hd.r_planet[i].z;
		plinfile >> hd.v_planet[i].x >> hd.v_planet[i].y >> hd.v_planet[i].z;

		if (readmomenta)
		{
			hd.v_planet[i].x /= hd.m_planet[i];
			hd.v_planet[i].y /= hd.m_planet[i];
			hd.v_planet[i].z /= hd.m_planet[i];
		}
	}

	icsinfile >> hd.n_part;
	if (max_particle > 0) hd.n_part = std::min(hd.n_part, max_particle);

	hd.r_part = Hvf64_3(hd.n_part);
	hd.v_part = Hvf64_3(hd.n_part);
	hd.part_flags = Hvu8(hd.n_part);
	hd.deathtime = Hvu32(hd.n_part);
	hd.id = Hvu32(hd.n_part);

	for (size_t i = 0; i < hd.n_part; i++)
	{
		icsinfile >> hd.r_part[i].x >> hd.r_part[i].y >> hd.r_part[i].z;
		icsinfile >> hd.v_part[i].x >> hd.v_part[i].y >> hd.v_part[i].z;

		std::string s;
		icsinfile >> s;
		if (!isdigit(s[0]))
		{
			icsinfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			hd.deathtime[i] = 0;
			hd.id[i] = i;
			hd.part_flags[i] = 0;
		}
		else
		{
			hd.deathtime[i] = std::stod(s);
			icsinfile >> hd.part_flags[i] >> hd.id[i];
		}
	}

	load_coefs(hd.coefdt, hd.dt);
	return false;
}

void save_data(const HostData& hd, const DeviceData& dd, std::string plout, std::string icsout, std::string infoout)
{
	std::ofstream ploutfile(plout), icsoutfile(icsout), infooutfile(infoout);

	ploutfile << hd.n_planet << std::endl;
	ploutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.n_planet; i++)
	{
		ploutfile << hd.m_planet[i] << std::endl;
		ploutfile << hd.r_planet[i].x << " " << hd.r_planet[i].y << " " << hd.r_planet[i].z << std::endl;
		ploutfile << hd.v_planet[i].x << " " << hd.v_planet[i].y << " " << hd.v_planet[i].z << std::endl;
	}

	icsoutfile << hd.n_part << std::endl;
	icsoutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.n_part; i++)
	{
		icsoutfile << hd.r_part[i].x << " " << hd.r_part[i].y << " " << hd.r_part[i].z << std::endl;
		icsoutfile << hd.v_part[i].x << " " << hd.v_part[i].y << " " << hd.v_part[i].z << std::endl;
		icsoutfile << hd.deathtime[i] << " " << hd.part_flags[i] << " " << hd.id[i] << std::endl;
	}

	infooutfile << std::setprecision(17);
	infooutfile << hd.t << " " << dd.n_part_alive << std::endl;
}

void transfer_data(const HostData& hd, DeviceData& dd)
{
	dd.log_buffer_id = 0;
	dd.phase_space_id = 0;

	dd.r_planet_log0 = Dvf64_3(hd.r_planet_log.size());
	dd.r_planet_log1 = Dvf64_3(hd.r_planet_log.size());
	dd.m_planet = Dvf64(hd.n_planet);

	dd.ps0.r = dd.ps0.v = Dvf64_3(hd.n_part);
	dd.ps0.flags = Dvu8(hd.n_part);
	dd.ps0.deathtime = Dvu32(hd.n_part);
	dd.ps0.id = Dvu32(hd.n_part);

	dd.ps1.r = dd.ps1.v = Dvf64_3(hd.n_part);
	dd.ps1.flags = Dvu8(hd.n_part);
	dd.ps1.deathtime = Dvu32(hd.n_part);
	dd.ps1.id = Dvu32(hd.n_part);

	dd.gather_indices = Dvu32(hd.n_part);

	dd.n_part_alive = 0;
	for (size_t i = 0; i < hd.n_part; i++)
	{
		if (~hd.part_flags[i] & 0x0001)
		{
			dd.n_part_alive++;
		}
	}

	dd.coefdt = Dvf64(hd.coefdt.size());

	DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;
	thrust::copy(hd.r_part.begin(), hd.r_part.end(), ps.r.begin());
	thrust::copy(hd.v_part.begin(), hd.v_part.end(), ps.v.begin());

	thrust::copy(hd.m_planet.begin(), hd.m_planet.end(), dd.m_planet.begin());
	thrust::copy(hd.coefdt.begin(), hd.coefdt.end(), dd.coefdt.begin());

	thrust::copy(hd.id.begin(), hd.id.end(), ps.id.begin());
	thrust::copy(hd.deathtime.begin(), hd.deathtime.end(), ps.deathtime.begin());
	thrust::copy(hd.part_flags.begin(), hd.part_flags.end(), ps.flags.begin());
}

void recover_data(HostData& hd, const DeviceData& dd, cudaStream_t& stream)
{
	const DevicePhaseSpace& ps = dd.phase_space_id % 2 ? dd.ps0 : dd.ps1;

	size_t n = hd.n_part;

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
					hd.r_part.begin(), hd.v_part.begin(),
					hd.part_flags.begin(), hd.deathtime.begin(), hd.id.begin())));
	 */

	thrust::copy(thrust::cuda::par.on(stream), ps.r.begin(), ps.r.begin() + n, hd.r_part.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.v.begin(), ps.v.begin() + n, hd.v_part.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.flags.begin(), ps.flags.begin() + n, hd.part_flags.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.deathtime.begin(), ps.deathtime.begin() + n, hd.deathtime.begin());
	thrust::copy(thrust::cuda::par.on(stream), ps.id.begin(), ps.id.begin() + n, hd.id.begin());
}

