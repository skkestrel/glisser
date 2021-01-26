#include "interp.h"
#include <iostream>
#include <fstream>
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
		use_bary_interp = config.use_bary_interp;

		pl = sr::data::HostPlanetPhaseSpace(config.interp_maxpl, config.tbsize);

		pl.m()[0] = sr::data::read_binary<float64_t>(input);
		sr::data::skip_binary(input, 32 - 8);

		aei0 = aei1 = oom0 = oom1 = bary_rs = bary_vs = Vf64_3(pl.n());
		m0 = m1 = bary_cms = Vf64(pl.n());

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


		std::vector<size_t> ids;
		std::vector<size_t> inds;
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
			ids.push_back(id);
			inds.push_back(ind);
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

		// Using barycentric orbital elements for interpolation
		if(use_bary_interp)
		{
			Vf64_3 plr, plv;
			plr.push_back(pl.r()[0]);
			plv.push_back(pl.v()[0]);


			// Calculate center of mass for all massive bodies.
			center_mass = pl.m()[0];
			f64_3 center_r = pl.r()[0] * pl.m()[0];
			f64_3 center_v = pl.v()[0] * pl.m()[0];

			size_t i = 0;
			// temp_log << std::setprecision(9) << t0 << std::setprecision(7) << " " << aei1[4] << " " << oom1[4] << std::endl;
			for (auto ind : inds)
			{	
				f64_3 aei = aei1[ind];
				f64_3 oom = oom1[ind];
				double gm = pl.m()[ind] + pl.m()[0];
				f64_3 r, v;
				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);
				// temp_log << "Original Helio " << std::setprecision(9) << t1 << std::setprecision(14) << " " << r << " " << v << std::endl;

				center_mass += pl.m()[ind];
				center_r += r * pl.m()[ind];
				center_v += v * pl.m()[ind];

				plr.push_back(r);
				plv.push_back(v);
				i++;
			}
			center_r /= center_mass;
			center_v /= center_mass;
			// temp_log  << "Original center " << std::setprecision(9) << t1 << std::setprecision(14) << " " << center_r << " " << center_v  << std::endl;

			i = 1;
			for (auto ind : inds)
			{	
				double a, e, in, capom, om, f;
				sr::convert::to_elements(sr::convert::get_bary_mu(center_mass, pl.m()[ind]), plr[i] - center_r, plv[i] - center_v,
					nullptr, &a, &e, &in, &capom, &om, &f);

				// temp_log << "Original Bary " << std::setprecision(9) << t1 << std::setprecision(14) << " " 
				// << plr[i] - center_r << " " << plv[i] - center_v << std::endl;

				double mean = sr::convert::get_mean_anomaly(e, f);
				aei1[ind] = f64_3(a, e, in);
				oom1[ind] = f64_3(capom, om, mean);
				i++;
			}
		}
		pl.r()[0] = f64_3(0);
		pl.v()[0] = f64_3(0);
		aei1[0] = f64_3(0);
		oom1[0] = f64_3(0);
		center_mass = 0;
		temp_log.open("temp_log.txt");
		
	}

	void Interpolator::fill_one(sr::data::HostPlanetPhaseSpace& pl, double relative_t)
	{
		pl.n_alive() = pl_alive + 1;

		if(use_bary_interp)
		{
			f64_3 center_r = f64_3(0);
			f64_3 center_v = f64_3(0);
			Vf64_3 plr_bary, plv_bary;
			// Recover the barycentric coords
			for (size_t j = 0; j < pl_alive; j++)
			{
				f64_3 aei_bary = reduced_aei_i[j] + reduced_daei[j] * relative_t;
				f64_3 oom_bary = reduced_oom_i[j] + reduced_doom[j] * relative_t;
				f64_3 r_bary, v_bary;
				
				sr::convert::from_elements_M(sr::convert::get_bary_mu(center_mass, reduced_m[j]), aei_bary.x, aei_bary.y, aei_bary.z, 
											oom_bary.x, oom_bary.y, oom_bary.z, &r_bary, &v_bary);
				plr_bary.push_back(r_bary);
				plv_bary.push_back(v_bary);
				center_r += r_bary * reduced_m[j];
				center_v += v_bary * reduced_m[j];
				// if (j==3) temp_log << std::setprecision(9) << t0 + relative_t << std::setprecision(10) << " " << aei_bary << " " << oom_bary << std::endl;
			}
			center_r /= pl.m()[0];
			center_v /= pl.m()[0];

			// Convert back to Heliocentric coords
			for (size_t j = 0; j < pl_alive; j++)
			{
				pl.r()[j + 1] = plr_bary[j] + center_r;
				pl.v()[j + 1] = plv_bary[j] + center_v;
				// temp_log << std::setprecision(9) << t0 + relative_t << std::setprecision(14) << " " << pl.r()[j + 1] << " " << pl.v()[j + 1] << std::endl;
			}
			
		}
		else
		{
			for (size_t j = 0; j < pl_alive; j++)
			{
				f64_3 aei = reduced_aei_i[j] + reduced_daei[j] * relative_t;
				f64_3 oom = reduced_oom_i[j] + reduced_doom[j] * relative_t;
				double gm = reduced_m[j] + pl.m()[0];
				f64_3 r, v;

				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);

				// offset for the sun
				pl.r()[j + 1] = r;
				pl.v()[j + 1] = v;
				// temp_log << std::setprecision(9) << t0 + relative_t << std::setprecision(14) << " " << pl.r()[j + 1] << " " << pl.v()[j + 1] << std::endl;
			}
		}
	}

	void Interpolator::fill(sr::data::HostPlanetPhaseSpace& pl, size_t nstep, double relative_t, double dt)
	{
		pl.n_alive_old() = pl.n_alive();

		// our n_alive doesn't include the sun
		pl.n_alive() = pl_alive + 1;

		// the first timestep starts at t + dt
		relative_t += dt;

		for (size_t j = 0; j < pl_alive; j++)
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

			if(use_bary_interp)
			{
				f64_3 center_r = f64_3(0);
				f64_3 center_v = f64_3(0);
				Vf64_3 plr_bary, plv_bary;
				// Recover the barycentric coords
				for (size_t j = 0; j < pl_alive; j++)
				{
					f64_3 aei_bary = reduced_aei_i[j] + reduced_daei[j] * relative_t;
					f64_3 oom_bary = reduced_oom_i[j] + reduced_doom[j] * relative_t;
					f64_3 r_bary, v_bary;
					
					sr::convert::from_elements_M(sr::convert::get_bary_mu(center_mass, reduced_m[j]), aei_bary.x, aei_bary.y, aei_bary.z, 
												oom_bary.x, oom_bary.y, oom_bary.z, &r_bary, &v_bary);
					plr_bary.push_back(r_bary);
					plv_bary.push_back(v_bary);
					center_r += r_bary * reduced_m[j];
					center_v += v_bary * reduced_m[j];
					// std::cout << std::setprecision(10) << bary_cms[j + 1] + reduced_m[j] << std::endl;
					// if (j==3) temp_log << std::setprecision(9) << t0 + relative_t << std::setprecision(10) << " " << aei_bary << " " << oom_bary << std::endl;
				}
				center_r /= pl.m()[0];
				center_v /= pl.m()[0];
				// temp_log  << "Recovered center " << std::setprecision(9) << t0 + relative_t << std::setprecision(14) << " " << center_r << " " << center_v  << std::endl;

				// Convert back to Heliocentric coords
				for (size_t j = 0; j < pl_alive; j++)
				{
					pl.r()[j + 1] = plr_bary[j] + center_r;
					pl.v()[j + 1] = plv_bary[j] + center_v;
					// temp_log << std::setprecision(9) << t0 + relative_t << std::setprecision(14) << " " << pl.r()[j + 1] << " " << pl.v()[j + 1] << std::endl;
				}
				
			}
			else 
			{
				for (size_t j = 0; j < pl_alive; j++)
				{
					f64_3 aei = reduced_aei_i[j] + reduced_daei[j] * relative_t;
					f64_3 oom = reduced_oom_i[j] + reduced_doom[j] * relative_t;
					f64_3 r, v;
					sr::convert::from_elements_M(pl.m()[0] + reduced_m[j], aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);
					pl.r()[j + 1] = r;
					pl.v()[j + 1] = v;
					// temp_log << std::setprecision(9) << t0 + relative_t << std::setprecision(14) << " " << pl.r()[j + 1] << " " << pl.v()[j + 1] << std::endl;
				}
			// if (j==3) temp_log << std::setprecision(9) << t0 + relative_t << std::setprecision(7) << " " << aei << " " << oom << std::endl;
			// offset for the sun
			}

			
			// temp_log << t0 + relative_t << " " << two_body_energy1 << " " << two_body_energy2 << std::endl;

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

		pl_alive_old = pl_alive;

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

		std::vector<size_t> ids;
		std::vector<size_t> inds;

		pl_alive = 0;
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
			ids.push_back(id);
			inds.push_back(ind);

			// Read heliocentric orbital elements of all planets.
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
			pl_alive++;
		}

		pl.n_alive() = pl_alive + 1;

		
		// Using barycentric orbital elements for interpolation
		if(use_bary_interp)
		{
			Vf64_3 plr, plv;
			plr.push_back(pl.r()[0]);
			plv.push_back(pl.v()[0]);


			// Calculate center of mass for all massive bodies.
			center_mass = pl.m()[0];
			f64_3 center_r = pl.r()[0] * pl.m()[0];
			f64_3 center_v = pl.v()[0] * pl.m()[0];

			size_t i = 0;
			// temp_log << std::setprecision(9) << t0 << std::setprecision(7) << " " << aei1[4] << " " << oom1[4] << std::endl;
			for (auto ind : inds)
			{	
				f64_3 aei = aei1[ind];
				f64_3 oom = oom1[ind];
				double gm = pl.m()[ind] + pl.m()[0];
				f64_3 r, v;
				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);
				// temp_log << "Original Helio " << std::setprecision(9) << t1 << std::setprecision(14) << " " << r << " " << v << std::endl;

				center_mass += pl.m()[ind];
				center_r += r * pl.m()[ind];
				center_v += v * pl.m()[ind];

				plr.push_back(r);
				plv.push_back(v);
				i++;
			}
			center_r /= center_mass;
			center_v /= center_mass;
			// temp_log  << "Original center " << std::setprecision(9) << t1 << std::setprecision(14) << " " << center_r << " " << center_v  << std::endl;

			i = 1;
			for (auto ind : inds)
			{	
				double a, e, in, capom, om, f;
				sr::convert::to_elements(sr::convert::get_bary_mu(center_mass, pl.m()[ind]), plr[i] - center_r, plv[i] - center_v,
					nullptr, &a, &e, &in, &capom, &om, &f);

				// temp_log << "Original Bary " << std::setprecision(9) << t1 << std::setprecision(14) << " " 
				// << plr[i] - center_r << " " << plv[i] - center_v << std::endl;

				double mean = sr::convert::get_mean_anomaly(e, f);
				aei1[ind] = f64_3(a, e, in);
				oom1[ind] = f64_3(capom, om, mean);
				i++;
			}
			// temp_log << std::setprecision(9) << t0 << std::setprecision(7) << " " << aei1[4] << " " << oom1[4] << std::endl;
		}

		size_t i = 0;
		for (auto id : ids)
		{

			size_t ind = inds[i];
			// i keeps track of the array index of the planets that are alive
			// THROUGHOUT the current interval - that is, the planet must be alive both at the beginning
			// and end of the interval
			// we use the THROUGHOUT-ALIVE particles to fill in planet history later
			reduced_ids[i] = id;

			// the mass is the mass at the end of the interval
			reduced_m[i] = m1[ind];

			reduced_aei_i[i] = aei0[ind];
			reduced_oom_i[i] = oom0[ind];
			reduced_aei_f[i] = aei1[ind];
			reduced_oom_f[i] = oom1[ind];

			reduced_daei[i] = (aei1[ind] - aei0[ind]) / dt;

			reduced_doom[i] = oom1[ind] - oom0[ind];
			if (reduced_doom[i].x > M_PI) reduced_doom[i].x -= 2 * M_PI;
			else if (reduced_doom[i].x < -M_PI) reduced_doom[i].x += 2 * M_PI;


			// FIXED IT!
			if (reduced_doom[i].y > M_PI) reduced_doom[i].y -= 2 * M_PI;
			else if (reduced_doom[i].y < -M_PI) reduced_doom[i].y += 2 * M_PI;

			reduced_doom[i] /= dt;

			// guess the mean motion frequency, a must be in AU and t in years

			double mmfreq = std::sqrt(pl.m()[ind] + pl.m()[0]) 
							* (std::pow(aei0[ind].x, -1.5) + std::pow(aei1[ind].x, -1.5)) / 2;
			// mmfreq -= (reduced_doom[i].x + reduced_doom[i].y);
			// std::cout << ind << " predicted orbital period: " << 2 * M_PI / mmfreq << std::endl;
			// std::cout << ind << " freq: " << mmfreq << std::endl;
			// std::cout << ind << " dt: " << dt << std::endl;
			// std::cout << ind << " oom0[ind].z: " << oom0[ind].z << std::endl;
			// std::cout << ind << " oom1[ind].z: " << oom1[ind].z << std::endl;

			double cmfin = oom0[ind].z + mmfreq * dt;
			cmfin = std::fmod(cmfin, 2 * M_PI);
			double corr = oom1[ind].z - cmfin;
			// std::cout << ind << " corr: " << corr << std::endl;
			if (corr > M_PI) corr -= 2 * M_PI;
			else if (corr < -M_PI) corr += 2 * M_PI;

			// correct the frequency to match final mean motion

			mmfreq += corr / dt;

			// std::cout << ind << " corr: " << corr << std::endl;
			// std::cout << ind << " freq: " << mmfreq << std::endl;

			reduced_doom[i].z = mmfreq;
			i++;
		}
		pl_alive = i;
	}
}
}
