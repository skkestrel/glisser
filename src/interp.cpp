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

		input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		resolve_encounters = config.resolve_encounters;
		user_dt = config.dt;
		fast_factor = config.resolve_encounters ? config.wh_ce_n1 * config.wh_ce_n2 : 1;

		pl.m()[0] = sr::data::read_binary<float64_t>(input) / 365.24 / 365.24;
		sr::data::skip_binary(input, 32 - 8);

		// note: can't add planets in the middle
#warning TODO how to do dynamic planet array sizing?

		aei0 = aei1 = oom0 = oom1 = aei_m1 = oom_m1 = daei = doom = Vf64_3(pl.n());
		mmfreq = std::vector<double>(pl.n());

#warning TODO need to think about time units
		t1 = sr::data::read_binary<float64_t>(input) * 365.24;
		t0 = std::numeric_limits<double>::infinity();
		t_m1 = std::numeric_limits<double>::quiet_NaN();

		pl.n_alive() = sr::data::read_binary<uint32_t>(input) + 1;
#warning TODO set old n alive or new n alive?
		pl.n_alive_old() = pl.n_alive();

		sr::data::skip_binary(input, 32 - 8 - 4);

		if (t1 > config.t_0)
		{
			throw std::runtime_error("time t0 not found in datafile: planet history starts after t0");
		}

		idmap[0] = 0;

		for (size_t i = 1; i < pl.n_alive(); i++)
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

			pl.m()[ind] = sr::data::read_binary<float32_t>(input) / 365.24 / 365.24;
			aei1[ind].x = sr::data::read_binary<float32_t>(input);
			aei1[ind].y = sr::data::read_binary<float32_t>(input);
			aei1[ind].z = sr::data::read_binary<float32_t>(input);
			oom1[ind].x = sr::data::read_binary<float32_t>(input);
			oom1[ind].y = sr::data::read_binary<float32_t>(input);
			oom1[ind].z = sr::data::read_binary<float32_t>(input);

			f64_3 r, v;
			sr::convert::from_elements_M(pl.m()[1] + pl.m()[ind], aei1[ind].x, aei1[ind].y, aei1[ind].z, oom1[ind].x, oom1[ind].y, oom1[ind].z, &r, &v);

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
		// the first timestep starts at t + dt
		relative_t += dt;

#warning TODO need to handle the case where planets disappear - need to be smart about array indices
		// currently cannot handle planets changing, and cannot handle planet indices switching around
		for (size_t i = 0; i < nstep * fast_factor; i++)
		{
			if (relative_t > t1 - t0 + 1e-8)
			{
				throw std::runtime_error("interpolation ran too long, timechunk exceeded lookup boundary");
			}

			pl.r()[0] = f64_3(0);
			pl.v()[0] = f64_3(0);
			for (size_t j = 1; j < pl.n_alive(); j++)
			{
				f64_3 aei = aei0[j] + daei[j] * relative_t;
				f64_3 oom = oom0[j] + doom[j] * relative_t;
				double gm = pl.m()[j] + pl.m()[0];

				f64_3 r, v;

				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);

				pl.r()[j] = r;
				pl.v()[j] = v;
			}


			Vf64_3* r_log = &pl.r_log().slow_old;
			Vf64_3* v_log = &pl.v_log().slow_old;
			if (resolve_encounters)
			{
				r_log = &pl.r_log().old;
				v_log = &pl.v_log().old;
			}

			std::copy(pl.r().begin() + 1, pl.r().begin() + pl.n_alive(), r_log->begin() + (pl.n_alive() - 1) * i);
			std::copy(pl.v().begin() + 1, pl.v().begin() + pl.n_alive(), v_log->begin() + (pl.n_alive() - 1) * i);

			if (resolve_encounters)
			{
				if ((i + 1) % (fast_factor) == 0)
				{
					size_t slow_index = i / fast_factor;

					auto fast_begin = pl.r_log().old.begin() + i * (pl.n() - 1);
					std::copy(fast_begin, fast_begin + (pl.n() - 1), pl.r_log().slow_old.begin() + slow_index * (pl.n() - 1));

					fast_begin = pl.v_log().old.begin() + i * (pl.n() - 1);
					std::copy(fast_begin, fast_begin + (pl.n() - 1), pl.v_log().slow_old.begin() + slow_index * (pl.n() - 1));
				}
			}

			if (resolve_encounters)
			{
				relative_t += dt / static_cast<double>(fast_factor);
			}
			else
			{
				relative_t += dt;
			}
		}

		pl.r_log().len_old = nstep;
		pl.v_log().len_old = nstep;
	}

	void Interpolator::next(sr::data::HostPlanetPhaseSpace& pl)
	{
		cur_ts = 0;

		aei_m1 = aei0;
		oom_m1 = oom0;
		t_m1 = t0;

		aei0 = aei1;
		oom0 = oom1;
		t0 = t1;
		rel_t = 0;

		t1 = sr::data::read_binary<float64_t>(input) * 365.24;


		if (!input)
		{
			throw EOSError();
		}

		pl.n_alive_old() = pl.n_alive();
		pl.n_alive() = sr::data::read_binary<uint32_t>(input) + 1;


		double dt = t1 - t0;

		// find the closest number of timesteps available and use that as the real timestep
		// user_dt is the "suggestion"
		n_ts = static_cast<size_t>(std::round(dt / user_dt));
		eff_dt = dt / static_cast<double>(n_ts);

		sr::data::skip_binary(input, 32 - 8 - 4);

		for (size_t i = 1; i < pl.n_alive(); i++)
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

			pl.m()[ind] = sr::data::read_binary<float32_t>(input) / 365.24 / 365.24;
			aei1[ind].x = sr::data::read_binary<float32_t>(input);
			aei1[ind].y = sr::data::read_binary<float32_t>(input);
			aei1[ind].z = sr::data::read_binary<float32_t>(input);
			oom1[ind].x = sr::data::read_binary<float32_t>(input);
			oom1[ind].y = sr::data::read_binary<float32_t>(input);
			oom1[ind].z = sr::data::read_binary<float32_t>(input);

			daei[ind] = (aei1[ind] - aei0[ind]) / dt;

			doom[ind] = oom1[ind] - oom0[ind];
			if (doom[ind].x > M_PI) doom[ind].x -= 2 * M_PI;
			else if (doom[ind].x < -M_PI) doom[ind].x += 2 * M_PI;

			if (doom[ind].y > M_PI) doom[ind].y -= 2 * M_PI;
			else if (doom[ind].x < -M_PI) doom[ind].y += 2 * M_PI;

			doom[ind] /= dt;

			// guess the mean motion frequency, a must be in AU and t in 

			mmfreq[ind] = 2 * M_PI * std::sqrt(1 + pl.m()[ind] / pl.m()[0])
				* (std::pow(aei0[ind].x, -1.5) + std::pow(aei1[ind].x, -1.5)) / 2;

			// std::cout << ind << " predicted orbital period: " << 2 * M_PI / mmfreq[ind] << std::endl;
			// std::cout << ind << " freq: " << mmfreq[ind] << std::endl;


			// mmfreq thinks time is in years so convert to days
			mmfreq[ind] /= 365.24;

			double cmfin = oom0[ind].z + mmfreq[ind] * dt;
			cmfin = std::fmod(cmfin, 2 * M_PI);
			double corr = oom1[ind].z - cmfin;
			if (corr > M_PI) corr -= 2 * M_PI;
			else if (corr < -M_PI) corr += 2 * M_PI;

			// correct the frequency to match final mean motion

			// std::cout << ind << " correction: " << corr/dt << std::endl;

			mmfreq[ind] += corr / dt;

			doom[ind].z = mmfreq[ind];
		}
	}
}
}
