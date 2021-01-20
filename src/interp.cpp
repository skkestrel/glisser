#include "interp.h"
#include <iostream>
#include <iomanip>

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

		aei0 = aei1 = oom0 = oom1 = Vf64_3(pl.n());
		m0 = m1 = Vf64(pl.n());

		// the reduced series of arrays holds data that is aligned - that is, aei_i and aei_f hold the SAME planets - they are gauranteed
		// to be the same set of planets in the same order, unlike aei0 and aei1 which may hold different sets of planets
		// i.e. the reduced sets are contain planets in the intersection of aei0 and aei1
		reduced_m = reduced_m_old = Vf64(pl.n());
		reduced_daei = reduced_doom = Vf64_3(pl.n());
		reduced_aei_i = reduced_aei_f = reduced_aei_i_old = reduced_aei_f_old = Vf64_3(pl.n());
		reduced_oom_i = reduced_oom_f = reduced_oom_i_old = reduced_oom_f_old = Vf64_3(pl.n());

		reduced_ids = reduced_ids_old = Vu32(pl.n());

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

			m1[ind] = sr::data::read_binary<float32_t>(input);
			aei1[ind].x = sr::data::read_binary<float32_t>(input);
			aei1[ind].y = sr::data::read_binary<float32_t>(input);
			aei1[ind].z = sr::data::read_binary<float32_t>(input);
			oom1[ind].x = sr::data::read_binary<float32_t>(input);
			oom1[ind].y = sr::data::read_binary<float32_t>(input);
			oom1[ind].z = sr::data::read_binary<float32_t>(input);

			f64_3 r, v;
			sr::convert::from_elements_M(pl.m()[0] + m1[ind], aei1[ind].x, aei1[ind].y, aei1[ind].z, oom1[ind].x, oom1[ind].y, oom1[ind].z, &r, &v);

			pl.r()[i] = r;
			pl.v()[i] = v;
			pl.id()[i] = id;
			pl.m()[i] = m1[ind];
		}
		pl.r()[0] = f64_3(0);
		pl.v()[0] = f64_3(0);
		aei1[0] = f64_3(0);
		oom1[0] = f64_3(0);
	}

	void Interpolator::fill_one(sr::data::HostPlanetPhaseSpace& pl, double relative_t)
	{
		pl.n_alive() = n_alive + 1;

		for (size_t j = 0; j < n_alive; j++)
		{
			f64_3 aei = reduced_aei_i[j] + reduced_daei[j] * relative_t;
			f64_3 oom = reduced_oom_i[j] + reduced_doom[j] * relative_t;
			double gm = reduced_m[j] + pl.m()[0];

			f64_3 r, v;

			sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);

			// offset for the sun
			pl.r()[j + 1] = r;
			pl.v()[j + 1] = v;
		}
	}

	void Interpolator::fill(sr::data::HostPlanetPhaseSpace& pl, size_t nstep, double relative_t, double dt)
	{
		pl.n_alive_old() = pl.n_alive();

		// our n_alive doesn't include the sun
		pl.n_alive() = n_alive + 1;

		// the first timestep starts at t + dt
		relative_t += dt;

		for (size_t j = 0; j < n_alive; j++)
		{
			pl.id()[j + 1] = reduced_ids[j];
			pl.m()[j + 1] = reduced_m[j];
		}

		for (size_t i = 0; i < nstep; i++)
		{
			if (relative_t > t1 - t0 + 1e-8)
			{
				throw std::runtime_error("interpolation ran too long, timechunk exceeded lookup boundary");
			}

			pl.r()[0] = f64_3(0);
			pl.v()[0] = f64_3(0);

			// std::cout << "pt = " << relative_t << std::endl;
			for (size_t j = 0; j < n_alive; j++)
			{
				f64_3 aei = reduced_aei_i[j] + reduced_daei[j] * relative_t;
				f64_3 oom = reduced_oom_i[j] + reduced_doom[j] * relative_t;
				double gm = reduced_m[j] + pl.m()[0];

				f64_3 r, v;

				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);

				// offset for the sun
				pl.r()[j + 1] = r;
				pl.v()[j + 1] = v;

				// std::cout << std::setprecision(17) << reduced_ids[j] << ": " << r << " " << v << std::endl;
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
		t_m1 = t0;

		m0 = m1;
		npl0 = npl1;
		aei0 = aei1;
		oom0 = oom1;
		t0 = t1;
		rel_t = 0;

		reduced_m_old = reduced_m;

		reduced_aei_i_old = reduced_aei_i;
		reduced_aei_f_old = reduced_aei_f;

		reduced_oom_i_old = reduced_oom_i;
		reduced_oom_f_old = reduced_oom_f;

		reduced_ids_old = reduced_ids;

		n_alive_old = n_alive;

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
			aei1[i] = f64_3(std::numeric_limits<double>::quiet_NaN());
			oom1[i] = f64_3(std::numeric_limits<double>::quiet_NaN());
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
			reduced_ids[interval_planet_index] = id;

			// the mass is the mass at the end of the interval
			reduced_m[interval_planet_index] = m1[ind];

			reduced_aei_i[interval_planet_index] = aei0[ind];
			reduced_oom_i[interval_planet_index] = oom0[ind];
			reduced_aei_f[interval_planet_index] = aei1[ind];
			reduced_oom_f[interval_planet_index] = oom1[ind];

			reduced_daei[interval_planet_index] = (aei1[ind] - aei0[ind]) / dt;

			reduced_doom[interval_planet_index] = oom1[ind] - oom0[ind];
			if (reduced_doom[interval_planet_index].x > M_PI) reduced_doom[interval_planet_index].x -= 2 * M_PI;
			else if (reduced_doom[interval_planet_index].x < -M_PI) reduced_doom[interval_planet_index].x += 2 * M_PI;

			if (reduced_doom[interval_planet_index].y > M_PI) reduced_doom[interval_planet_index].y -= 2 * M_PI;
			else if (reduced_doom[interval_planet_index].y < -M_PI) reduced_doom[interval_planet_index].y += 2 * M_PI;

			reduced_doom[interval_planet_index] /= dt;

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

			mmfreq += corr / dt;

			reduced_doom[interval_planet_index].z = mmfreq;

			interval_planet_index++;
		}

		n_alive = interval_planet_index;
	}
}
}
