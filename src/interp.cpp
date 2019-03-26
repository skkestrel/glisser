#include "interp.h"

namespace sr
{
namespace interp
{
	Interpolator::Interpolator(const sr::data::Configuration& config, sr::data::HostPlanetPhaseSpace& pl, std::string file)
		: input(file, std::ios_base::binary)
	{
		resolve_encounters = config.resolve_encounters;
		fast_factor = config.wh_ce_n1 * config.wh_ce_n2;

		pl.m()[0] = sr::data::read_binary<float64_t>(input);
		sr::data::skip_binary(input, 32 - 8);

		// note: can't add planets in the middle

		// TODO time
		t0 = sr::data::read_binary<float64_t>(input) * 365.24;
		pl.n_alive() = sr::data::read_binary<uint32_t>(input) + 1;

		sr::data::skip_binary(input, 32 - 8 - 4);

		if (t0 > config.t_0)
		{
			throw std::runtime_error("time t0 not found in datafile");
		}

		idmap[0] = 0;

		for (size_t i = 0; i < pl.n(); i++)
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
			aei0[ind].x = sr::data::read_binary<float32_t>(input);
			aei0[ind].y = sr::data::read_binary<float32_t>(input);
			aei0[ind].z = sr::data::read_binary<float32_t>(input);
			oof0[ind].x = sr::data::read_binary<float32_t>(input);
			oof0[ind].y = sr::data::read_binary<float32_t>(input);
			oof0[ind].z = sr::data::read_binary<float32_t>(input);

			f64_3 r, v;
			sr::convert::from_elements_M(pl.m()[0] + pl.m()[ind], aei0[ind].x, aei0[ind].y, aei0[ind].z, oof0[ind].x, oof0[ind].y, oof0[ind].z, &r, &v);

			pl.r()[ind] = r;
			pl.v()[ind] = v;
		}
		aei0[0] = f64_3(0);
		oof0[0] = f64_3(0);
	}

	void Interpolator::fill(sr::data::HostPlanetPhaseSpace& pl, size_t nstep, double t, double dt)
	{
#warning TODO which dt should we use?
#warning TODO need to handle the case where planets disappear - need to be smart about array indices
		// currently cannot handle planets changing, and cannot handle planet indices switching around
		for (size_t i = 0; i < nstep * fast_factor; i++)
		{
			pl.r()[0] = f64_3(0);
			pl.v()[0] = f64_3(0);
			for (size_t j = 1; j < pl.n(); j++)
			{
				f64_3 aei = aei0[j] + daei[j] * (t + dt - t0);
				f64_3 oof = oof0[j] + doof[j] * (t + dt - t0);
				double gm = pl.m()[j] + pl.m()[0];

				f64_3 r, v;
				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oof.x, oof.y, oof.z, &r, &v);
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

			std::copy(pl.r().begin() + 1, pl.r().end(), r_log->begin() + (pl.n_alive() - 1) * i);
			std::copy(pl.v().begin() + 1, pl.v().end(), v_log->begin() + (pl.n_alive() - 1) * i);

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

			if (t > t0) next(pl);
		}

		pl.r_log().len_old = nstep;
		pl.v_log().len_old = nstep;
	}

	void Interpolator::next(sr::data::HostPlanetPhaseSpace& pl)
	{
		t1 = sr::data::read_binary<float64_t>(input) * 365.24;
		if (!input)
		{
			throw EOSError();
		}

		pl.n_alive_old() = pl.n_alive();
		pl.n_alive() = sr::data::read_binary<uint32_t>(input) + 1;

		double dt = t1 - t0;

		sr::data::skip_binary(input, 32 - 8 - 4);
		for (size_t i = 1; i < pl.n() + 1; i++)
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
			oof1[ind].x = sr::data::read_binary<float32_t>(input);
			oof1[ind].y = sr::data::read_binary<float32_t>(input);
			oof1[ind].z = sr::data::read_binary<float32_t>(input);

			f64_3 r, v;
			
			daei[ind] = (aei1[ind] - aei0[ind]) / dt;

			doof[ind] = oof1[ind] - oof0[ind];
			if (doof[ind].x > M_PI) doof[ind].x -= 2 * M_PI;
			else if (doof[ind].x < -M_PI) doof[ind].x += 2 * M_PI;

			if (doof[ind].y > M_PI) doof[ind].y -= 2 * M_PI;
			else if (doof[ind].x < -M_PI) doof[ind].y += 2 * M_PI;

			doof[ind] /= dt;

			// guess the mean motion frequency, a must be in AU and t in 

			mmfreq[ind] = 2 * M_PI * std::sqrt(1 + pl.m()[ind] / pl.m()[0])
				* (std::pow(aei0[ind].x, -1.5) + std::pow(aei1[ind].x, -1.5)) / 2;


			// mmfreq thinks time is in years so convert to days
			mmfreq[ind] /= 365.24;

			double cmfin = oof0[ind].z + mmfreq[ind] * dt;
			cmfin = std::fmod(cmfin, 2 * M_PI);
			double corr = oof1[ind].z - cmfin;
			if (corr > M_PI) corr -= 2 * M_PI;
			else if (corr < -M_PI) corr -= 2 * M_PI;

			// correct the frequency to match final mean motion

			mmfreq[ind] += corr / dt;
		}

		aei0 = aei1;
		oof0 = oof1;
		t0 = t1;
	}
}
}
