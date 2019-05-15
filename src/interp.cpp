#include "interp.h"
#include <iostream>

namespace sr
{
namespace interp
{
	Interpolator::Interpolator()
	{
	}

	Interpolator::Interpolator(const sr::data::Configuration& config, sr::data::HostPlanetPhaseSpace& pl, std::string file)
		: input(file, std::ios_base::binary)
	{
		cur_ts = 0;
		user_dt = config.dt;

		pl = sr::data::HostPlanetPhaseSpace(config.interp_maxpl, config.tbsize);

		pl.m()[0] = sr::data::read_binary<float64_t>(input);
		sr::data::skip_binary(input, 32 - 8);

		aei0 = aei1 = oom0 = oom1 = aei_m1 = oom_m1 = daei = doom = aei_i = oom_i = Vf64_3(pl.n());

		t0 = std::numeric_limits<double>::infinity();

		t1 = sr::data::read_binary<float64_t>(input);
		t_m1 = std::numeric_limits<double>::quiet_NaN();
		npl1 = sr::data::read_binary<uint32_t>(input) + 1;

		pl.n_alive() = npl1;
		pl.n_alive_old() = npl1;

		sr::data::skip_binary(input, 32 - 8 - 4);

		if (t1 > config.t_0)
		{
			throw std::runtime_error("time t0 not found in datafile: planet history starts after t0");
		}

		idmap[0] = 0;

		for (size_t i = 1; i < npl1; i++)
		{
			uint32_t id = sr::data::read_binary<uint32_t>(input);
			size_t ind;
			if (idmap.find(id) == idmap.end())
			{
				size_t sz = idmap.size();
				idmap[id] = sz;
				ind = sz;
				pl.id()[ind] = id;
			}
			else
			{
				ind = idmap[id];
			}

			pl.m()[ind] = sr::data::read_binary<float32_t>(input);
			aei1[ind].x = sr::data::read_binary<float32_t>(input);
			aei1[ind].y = sr::data::read_binary<float32_t>(input);
			aei1[ind].z = sr::data::read_binary<float32_t>(input);
			oom1[ind].x = sr::data::read_binary<float32_t>(input);
			oom1[ind].y = sr::data::read_binary<float32_t>(input);
			oom1[ind].z = sr::data::read_binary<float32_t>(input);

			f64_3 r, v;
			sr::convert::from_elements_M(pl.m()[0] + pl.m()[ind], aei1[ind].x, aei1[ind].y, aei1[ind].z, oom1[ind].x, oom1[ind].y, oom1[ind].z, &r, &v);

			pl.r()[ind] = r;
			pl.v()[ind] = v;
		}
		pl.r()[0] = f64_3(0);
		pl.v()[0] = f64_3(0);
		aei1[0] = f64_3(0);
		oom1[0] = f64_3(0);
	}

	void Interpolator::fill(sr::data::HostPlanetPhaseSpace& pl, size_t nstep, double relative_t, double dt)
	{
		pl.n_alive_old() = pl.n_alive();

		// our n_alive doesn't include the sun
		pl.n_alive() = n_alive + 1;

		// the first timestep starts at t + dt
		relative_t += dt;

#warning TODO need to handle the case where planets disappear - need to be smart about array indices
		// currently cannot handle planets changing, and cannot handle planet indices switching around
		for (size_t i = 0; i < nstep; i++)
		{
			if (relative_t > t1 - t0 + 1e-8)
			{
				throw std::runtime_error("interpolation ran too long, timechunk exceeded lookup boundary");
			}

			pl.r()[0] = f64_3(0);
			pl.v()[0] = f64_3(0);
			for (size_t j = 0; j < n_alive; j++)
			{
				f64_3 aei = aei_i[j] + daei[j] * relative_t;
				f64_3 oom = oom_i[j] + doom[j] * relative_t;
				double gm = pl.m()[j] + pl.m()[0];

				f64_3 r, v;

				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);

				// offset for the sun
				pl.r()[j + 1] = r;
				pl.v()[j + 1] = v;
			}

			std::copy(pl.r().begin() + 1, pl.r().begin() + pl.n_alive(), pl.r_log().old.begin() + (pl.n_alive() - 1) * i);
			std::copy(pl.v().begin() + 1, pl.v().begin() + pl.n_alive(), pl.v_log().old.begin() + (pl.n_alive() - 1) * i);

			relative_t += dt;
		}

		pl.r_log().len_old = nstep;
		pl.v_log().len_old = nstep;
	}

	void Interpolator::next(sr::data::HostPlanetPhaseSpace& pl)
	{
		aei_m1 = aei0;
		oom_m1 = oom0;
		t_m1 = t0;

		npl0 = npl1;
		aei0 = aei1;
		oom0 = oom1;
		t0 = t1;
		rel_t = 0;

		t1 = sr::data::read_binary<float64_t>(input);

		if (!input)
		{
			throw EOSError();
		}

		npl1 = sr::data::read_binary<uint32_t>(input) + 1;

		double dt = t1 - t0;

		// find the closest number of timesteps available and use that as the real timestep
		// user_dt is the "suggestion"
		n_ts = std::max<size_t>(1, static_cast<size_t>(std::round(dt / user_dt)));
		eff_dt = dt / static_cast<double>(n_ts);

		sr::data::skip_binary(input, 32 - 8 - 4);

		for (size_t i = 1; i < pl.n(); i++)
		{
			// initialize with nans so we know which entries are empty
			aei1[ind] = f64_3(std::numeric_limits<double>::quiet_nan());
			oom1[ind] = f64_3(std::numeric_limits<double>::quiet_nan());
		}

		size_t interval_planet_index = 0;
		for (size_t i = 1; i < npl1; i++)
		{
			uint32_t id = sr::data::read_binary<uint32_t>(input);
			size_t ind;
			if (idmap.find(id) == idmap.end())
			{
				ind = (idmap[id] = idmap.size());
				pl.id()[ind] = id;
			}
			else
			{
				ind = idmap[id];
			}

			pl.m()[ind] = sr::data::read_binary<float32_t>(input);
			aei1[ind].x = sr::data::read_binary<float32_t>(input);
			aei1[ind].y = sr::data::read_binary<float32_t>(input);
			aei1[ind].z = sr::data::read_binary<float32_t>(input);
			oom1[ind].x = sr::data::read_binary<float32_t>(input);
			oom1[ind].y = sr::data::read_binary<float32_t>(input);
			oom1[ind].z = sr::data::read_binary<float32_t>(input);

			// if there was no entry for the planet at the beginning of the interval, skip
			if (std::isnan(aei0[ind].x))
			{
				continue;
			}

			// interval_planet_index keeps track of the array index of the planets that are alive
			// THROUGHOUT the current interval - that is, the planet must be alive both at the beginning
			// and end of the interval
			// we use the THROUGHOUT-ALIVE particles to fill in planet history later
			id[interval_planet_index] = id;
			aei_i[interval_planet_index] = aei0[ind];
			oom_i[interval_planet_index] = oom0[ind];

			daei[interval_planet_index] = (aei1[ind] - aei0[ind]) / dt;

			doom[interval_planet_index] = oom1[ind] - oom0[ind];
			if (doom[interval_planet_index].x > M_PI) doom[interval_planet_index].x -= 2 * M_PI;
			else if (doom[interval_planet_index].x < -M_PI) doom[interval_planet_index].x += 2 * M_PI;

			if (doom[interval_planet_index].y > M_PI) doom[interval_planet_index].y -= 2 * M_PI;
			else if (doom[interval_planet_index].x < -M_PI) doom[interval_planet_index].y += 2 * M_PI;

			doom[interval_planet_index] /= dt;

			// guess the mean motion frequency, a must be in AU and t in 

			double mmfreq = 2 * M_PI * std::sqrt(1 + pl.m()[ind] / pl.m()[0])
				* (std::pow(aei0[ind].x, -1.5) + std::pow(aei1[ind].x, -1.5)) / 2;

			// std::cout << ind << " predicted orbital period: " << 2 * M_PI / mmfreq << std::endl;
			// std::cout << ind << " freq: " << mmfreq << std::endl;

			double cmfin = oom0[ind].z + mmfreq * dt;
			cmfin = std::fmod(cmfin, 2 * M_PI);
			double corr = oom1[ind].z - cmfin;
			if (corr > M_PI) corr -= 2 * M_PI;
			else if (corr < -M_PI) corr += 2 * M_PI;

			// correct the frequency to match final mean motion

			// std::cout << ind << " correction: " << corr/dt << std::endl;

			mmfreq += corr / dt;

			doom[interval_planet_index].z = mmfreq;

			interval_planet_index++;
		}

		n_alive = interval_planet_index;
	}
}
}
