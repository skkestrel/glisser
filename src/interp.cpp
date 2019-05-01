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
		input.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		resolve_encounters = config.resolve_encounters;
		fast_factor = config.resolve_encounters ? config.wh_ce_n1 * config.wh_ce_n2 : 1;

		pl.m()[0] = sr::data::read_binary<float64_t>(input) / 365.24 / 365.24;
		sr::data::skip_binary(input, 32 - 8);

		// note: can't add planets in the middle
#warning TODO how to do dynamic planet array sizing? can't use pl.n() for interpolation

		aei0 = aei1 = oom0 = oom1 = daei = doom = Vf64_3(pl.n());
		mmfreq = std::vector<double>(pl.n());

#warning TODO need to fink about time units
		t0 = sr::data::read_binary<float64_t>(input) * 365.24;
		t1 = t0 - 1;


		pl.n_alive() = sr::data::read_binary<uint32_t>(input) + 1;
#warning TODO set old n alive or new n alive?
		pl.n_alive_old() = pl.n_alive();

		sr::data::skip_binary(input, 32 - 8 - 4);

		if (t0 > config.t_0)
		{
			throw std::runtime_error("time t0 not found in datafile");
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

	void Interpolator::fill(sr::data::HostPlanetPhaseSpace& pl, size_t nstep, double t, double dt)
	{
#warning TODO which dt should we use?
#warning TODO need to handle the case where planets disappear - need to be smart about array indices
		// currently cannot handle planets changing, and cannot handle planet indices switching around
		for (size_t i = 0; i < nstep * fast_factor; i++)
		{
			while (t > t1) next(pl);

			pl.r()[0] = f64_3(0);
			pl.v()[0] = f64_3(0);
			for (size_t j = 1; j < pl.n_alive(); j++)
			{
				f64_3 aei = aei0[j] + daei[j] * (t + dt - t0);
				f64_3 oom = oom0[j] + doom[j] * (t + dt - t0);
				double gm = pl.m()[j] + pl.m()[0];

				f64_3 r, v;

				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);

				/*
				std::cout << "s0: " << aei0[j] << " " << oom0[j] << std::endl;
				std::cout << "s1: " << aei1[j] << " " << oom1[j] << std::endl;
				std::cout << "st: " << aei << " " << oom << std::endl;
				std::cout << "rv: " << r << " " << v << std::endl;

				double a, e, I, O, o, f;

				sr::convert::to_elements(pl.m()[0] + pl.m()[j], r, v, nullptr, &a, &e, &I, &O, &o, &f);

				std::cout << "back: " << a << " " << e << " " << I << " " << O << " " << o << " " << f << std::endl;
				int x;
				std::cin >> x;
				*/


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
				t += dt / static_cast<double>(fast_factor);
			}
			else
			{
				t += dt;
			}
		}

		pl.r_log().len_old = nstep;
		pl.v_log().len_old = nstep;
	}

	void Interpolator::next(sr::data::HostPlanetPhaseSpace& pl)
	{
		aei0 = aei1;
		oom0 = oom1;
		t0 = t1;
		t1 = sr::data::read_binary<float64_t>(input) * 365.24;

		if (!input)
		{
			throw EOSError();
		}

		pl.n_alive_old() = pl.n_alive();
		pl.n_alive() = sr::data::read_binary<uint32_t>(input) + 1;


		double dt = t1 - t0;

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
