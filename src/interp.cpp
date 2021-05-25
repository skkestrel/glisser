#include "interp.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

namespace sr
{
namespace interp
{
	Interpolator::Interpolator()
	{
	}

	Interpolator::Interpolator(const sr::data::Configuration& config, sr::data::HostPlanetPhaseSpace& pl, std::string file, bool _binary_hist, bool _single_precision)
		: input(file), binary_hist(_binary_hist), single_precision(_single_precision)
	{
		// temp_log.open("temp_log_interp.txt");
		// temp_log << std::setprecision(17);
		cur_ts = 0;
		user_dt = config.dt;
		use_jacobi_interp = config.use_jacobi_interp;

		pl = sr::data::HostPlanetPhaseSpace(config.interp_maxpl, config.tbsize);

		aei0 = aei1 = oom0 = oom1 = bary_rs = bary_vs = Vf64_3(pl.n());
		m0 = m1 = bary_cms = Vf64(pl.n());
		rplsq0 = rplsq1 = Vf64(pl.n());

		// the reduced series of arrays holds data that is aligned - that is, aei_i and aei_f hold the SAME planets - they are gauranteed
		// to be the same set of planets in the same order, unlike aei0 and aei1 which may hold different sets of planets
		// i.e. the reduced sets are contain planets in the intersection of aei0 and aei1
		reduced_ids = reduced_ids_old = Vu32(pl.n());
		reduced_m = reduced_m_old = Vf64(pl.n());
		reduced_rplsq = reduced_rplsq_old = Vf64(pl.n());
		reduced_daei = reduced_doom = Vf64_3(pl.n());
		reduced_aei_i = reduced_aei_f = reduced_aei_i_old = reduced_aei_f_old = Vf64_3(pl.n());
		reduced_oom_i = reduced_oom_f = reduced_oom_i_old = reduced_oom_f_old = Vf64_3(pl.n());
		

		jacobi_aei_i = jacobi_aei_f = jacobi_aei_i_old = jacobi_aei_f_old = Vf64_3(pl.n());
		jacobi_oom_i = jacobi_oom_f = jacobi_oom_i_old = jacobi_oom_f_old = Vf64_3(pl.n());

		planet_eta = Vf64(pl.n());
		planet_rj = planet_vj = Vf64_3(pl.n());	

		t0 = std::numeric_limits<double>::infinity(); 
		t_m1 = std::numeric_limits<double>::quiet_NaN();

		if (binary_hist)
		{
			t1 = sr::data::read_binary<float64_t>(input);
			pl.m()[0] = sr::data::read_binary<float64_t>(input);
			npl1 = sr::data::read_binary<uint32_t>(input) + 1; 
		}
		else
		{
			input >> t1 >> pl.m()[0] >> npl1;
			npl1 += 1;
		}
						
		pl.n_alive() = npl1;
		pl.n_alive_old() = npl1;

		if (t1 > config.t_0)
		{
			throw std::runtime_error("time t0 not found in datafile: planet history starts after t0");
		}

		idmap[0] = 0;


		std::vector<size_t> ids;
		std::vector<size_t> inds;
		for (size_t i = 1; i < npl1; i++)
		{
			uint32_t id;
			if (binary_hist)
			{
				id = sr::data::read_binary<uint32_t>(input);
			}
			else
			{
				input >> id;
			}

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

			if(binary_hist)
			{
				if (single_precision)
				{
					m1[ind] = sr::data::read_binary<float32_t>(input);
					rplsq1[ind] = sr::data::read_binary<float32_t>(input);
					aei1[ind].x = sr::data::read_binary<float32_t>(input);
					aei1[ind].y = sr::data::read_binary<float32_t>(input);
					aei1[ind].z = sr::data::read_binary<float32_t>(input);
					oom1[ind].x = sr::data::read_binary<float32_t>(input);
					oom1[ind].y = sr::data::read_binary<float32_t>(input);
					oom1[ind].z = sr::data::read_binary<float32_t>(input);
				}
				else
				{
					m1[ind] = sr::data::read_binary<double>(input);
					rplsq1[ind] = sr::data::read_binary<double>(input);
					aei1[ind].x = sr::data::read_binary<double>(input);
					aei1[ind].y = sr::data::read_binary<double>(input);
					aei1[ind].z = sr::data::read_binary<double>(input);
					oom1[ind].x = sr::data::read_binary<double>(input);
					oom1[ind].y = sr::data::read_binary<double>(input);
					oom1[ind].z = sr::data::read_binary<double>(input);
				}
			}
			else
			{
				input >> m1[ind] >> rplsq1[ind]>> aei1[ind].x >> aei1[ind].y >> aei1[ind].z
						  >> oom1[ind].x >> oom1[ind].y >> oom1[ind].z;
			}

			f64_3 r, v;
			sr::convert::from_elements_M(pl.m()[0] + m1[ind], aei1[ind].x, aei1[ind].y, aei1[ind].z, oom1[ind].x, oom1[ind].y, oom1[ind].z, &r, &v);

			pl.id()[i] = id;
			pl.m()[i] = m1[ind];
			pl.rplsq()[i] = rplsq1[ind];
			pl.r()[i] = r;
			pl.v()[i] = v;
		}

		// Using barycentric orbital elements for interpolation
		if(use_jacobi_interp) 
		{
			planet_eta[0] = pl.m()[0];
			for (size_t i = 1; i < pl.n_alive(); i++)
			{
				planet_eta[i] = planet_eta[i - 1] + pl.m()[i];
			}
			sr::convert::helio_to_jacobi_r_planets(pl, planet_eta, planet_rj);
			sr::convert::helio_to_jacobi_v_planets(pl, planet_eta, planet_vj);
			size_t i = 1;
			for (auto ind : inds)
			{	
				double a, e, in, capom, om, f;
				sr::convert::to_elements(planet_eta[i], planet_rj[i], planet_vj[i],
					nullptr, &a, &e, &in, &capom, &om, &f);

				aei1[ind] = f64_3(a, e, in);
				oom1[ind] = f64_3(capom, om, sr::convert::get_mean_anomaly(e, f));
				jacobi_aei_f[i-1] = aei1[ind];
				jacobi_oom_f[i-1] = oom1[ind];
				// temp_log  << t1 << " " << pl.r()[i] << std::endl;
				i++;
			}

		}
		pl.r()[0] = f64_3(0);
		pl.v()[0] = f64_3(0);
		aei1[0] = f64_3(0);
		oom1[0] = f64_3(0);
		center_mass = 0;
		
	}

	void Interpolator::fill_one(sr::data::HostPlanetPhaseSpace& pl, double relative_t)
	{
		pl.n_alive() = pl_alive + 1;

		if(use_jacobi_interp)
		{
			planet_rj[0] = f64_3(0);
			planet_vj[0] = f64_3(0);
			for (size_t j = 0; j < pl_alive; j++)
			{
				f64_3 aei_jacobi = reduced_aei_i[j] + reduced_daei[j] * relative_t;
				f64_3 oom_jacobi = reduced_oom_i[j] + reduced_doom[j] * relative_t;
				f64_3 r_jacobi, v_jacobi;
				
				sr::convert::from_elements_M(planet_eta[j + 1], aei_jacobi.x, aei_jacobi.y, aei_jacobi.z, 
											oom_jacobi.x, oom_jacobi.y, oom_jacobi.z, &r_jacobi, &v_jacobi);
				planet_rj[j + 1] = r_jacobi;
				planet_vj[j + 1] = v_jacobi;
				// if (j==pl_alive-1) temp_log << t0 + relative_t << " " << aei_jacobi << " " << oom_jacobi << std::endl;
			}

			sr::convert::jacobi_to_helio_planets(planet_eta, planet_rj, planet_vj, pl);
			// temp_log << t0 + relative_t << " " << pl.r()[pl_alive] << std::endl;
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
			pl.rplsq()[j + 1] = reduced_rplsq[j];
		}
		for (size_t i = 0; i < nstep; i++)
		{
			if (relative_t > t1 - t0 + 1e-8)
			{
				throw std::runtime_error("interpolation ran too long, timechunk exceeded lookup boundary");
			}

			pl.r()[0] = f64_3(0);
			pl.v()[0] = f64_3(0);

			if(use_jacobi_interp)
			{
				planet_rj[0] = f64_3(0);
				planet_vj[0] = f64_3(0);
				for (size_t j = 0; j < pl_alive; j++)
				{
					f64_3 aei_jacobi = reduced_aei_i[j] + reduced_daei[j] * relative_t;
					f64_3 oom_jacobi = reduced_oom_i[j] + reduced_doom[j] * relative_t;
					f64_3 r_jacobi, v_jacobi;
					
					sr::convert::from_elements_M(planet_eta[j + 1], aei_jacobi.x, aei_jacobi.y, aei_jacobi.z, 
												oom_jacobi.x, oom_jacobi.y, oom_jacobi.z, &r_jacobi, &v_jacobi);
					planet_rj[j + 1] = r_jacobi;
					planet_vj[j + 1] = v_jacobi;
					// if (j==pl_alive-1) temp_log << t0 + relative_t << " " << aei_jacobi << " " << oom_jacobi << std::endl;
				}

				sr::convert::jacobi_to_helio_planets(planet_eta, planet_rj, planet_vj, pl);
				// temp_log << t0 + relative_t << " " << pl.r()[pl_alive] << std::endl;
				
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
				}
			// offset for the sun
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
		rplsq0 = rplsq1;
		npl0 = npl1;
		aei0 = aei1;
		oom0 = oom1;
		t0 = t1;
		rel_t = 0;

		reduced_ids_old = reduced_ids;
		reduced_m_old = reduced_m;
		reduced_rplsq_old = reduced_rplsq;
		
		reduced_aei_i_old = reduced_aei_i;
		reduced_aei_f_old = reduced_aei_f;
		reduced_oom_i_old = reduced_oom_i;
		reduced_oom_f_old = reduced_oom_f;

		jacobi_aei_i_old = jacobi_aei_i;
		jacobi_aei_f_old = jacobi_aei_f;
		jacobi_oom_i_old = jacobi_oom_i;
		jacobi_oom_f_old = jacobi_oom_f;

		jacobi_aei_i = jacobi_aei_f;
		jacobi_oom_i = jacobi_oom_f;
		pl_alive_old = pl_alive;

		if (binary_hist)
		{
			t1 = sr::data::read_binary<float64_t>(input);
			pl.m()[0] = sr::data::read_binary<float64_t>(input);
			npl1 = sr::data::read_binary<uint32_t>(input) + 1; 
		}
		else
		{
			input >> t1 >> pl.m()[0] >> npl1;
			npl1 += 1;
		}

		if (!input)
		{
			throw EOSError();
		}


		double dt = t1 - t0;

		// find the closest number of timesteps available and use that as the real timestep
		// user_dt is the "suggestion"
		n_ts = std::max<size_t>(1, static_cast<size_t>(std::round(dt / user_dt)));
		eff_dt = dt / static_cast<double>(n_ts);


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
			uint32_t id;
			if (binary_hist)
			{
				id = sr::data::read_binary<uint32_t>(input);
			}
			else
			{
				input >> id;
			}
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
			if(binary_hist)
			{
				if (single_precision)
				{
					m1[ind] = sr::data::read_binary<float32_t>(input);
					rplsq1[ind] = sr::data::read_binary<float32_t>(input);
					aei1[ind].x = sr::data::read_binary<float32_t>(input);
					aei1[ind].y = sr::data::read_binary<float32_t>(input);
					aei1[ind].z = sr::data::read_binary<float32_t>(input);
					oom1[ind].x = sr::data::read_binary<float32_t>(input);
					oom1[ind].y = sr::data::read_binary<float32_t>(input);
					oom1[ind].z = sr::data::read_binary<float32_t>(input);
				}
				else
				{
					m1[ind] = sr::data::read_binary<double>(input);
					rplsq1[ind] = sr::data::read_binary<double>(input);
					aei1[ind].x = sr::data::read_binary<double>(input);
					aei1[ind].y = sr::data::read_binary<double>(input);
					aei1[ind].z = sr::data::read_binary<double>(input);
					oom1[ind].x = sr::data::read_binary<double>(input);
					oom1[ind].y = sr::data::read_binary<double>(input);
					oom1[ind].z = sr::data::read_binary<double>(input);
				}
			}
			else
			{
				input >> m1[ind] >> rplsq1[ind] >> aei1[ind].x >> aei1[ind].y >> aei1[ind].z
						  >> oom1[ind].x >> oom1[ind].y >> oom1[ind].z;
			}



			// if there was no entry for the planet at the beginning of the interval, skip
			if (std::isnan(aei0[ind].x))
			{
				continue;
			}
			pl_alive++;
			
		}

		pl.n_alive() = pl_alive + 1;

		// Using barycentric orbital elements for interpolation

		if(use_jacobi_interp)
		{
			size_t i = 0;
			for (auto ind : inds)
			{	
				f64_3 aei = aei1[ind];
				f64_3 oom = oom1[ind];
				double gm = pl.m()[ind] + pl.m()[0];
				f64_3 r, v;
				sr::convert::from_elements_M(gm, aei.x, aei.y, aei.z, oom.x, oom.y, oom.z, &r, &v);

				pl.r()[i + 1] = r;
				pl.v()[i + 1] = v;
				i++;
			}


			planet_eta[0] = pl.m()[0];
			for (i = 1; i < pl.n_alive(); i++)
			{
				planet_eta[i] = planet_eta[i - 1] + pl.m()[i];
			}
			sr::convert::helio_to_jacobi_r_planets(pl, planet_eta, planet_rj);
			sr::convert::helio_to_jacobi_v_planets(pl, planet_eta, planet_vj);

			i = 1;
			for (auto ind : inds)
			{	
				double a, e, in, capom, om, f;
				sr::convert::to_elements(planet_eta[i], planet_rj[i], planet_vj[i],
					nullptr, &a, &e, &in, &capom, &om, &f);

				aei1[ind] = f64_3(a, e, in);
				oom1[ind] = f64_3(capom, om, sr::convert::get_mean_anomaly(e, f));
				jacobi_aei_f[i-1] = aei1[ind];
				jacobi_oom_f[i-1] = oom1[ind];
				i++;
			}
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
			reduced_rplsq[i] = rplsq1[ind];

			reduced_aei_i[i] = aei0[ind];
			reduced_oom_i[i] = oom0[ind];
			reduced_aei_f[i] = aei1[ind];
			reduced_oom_f[i] = oom1[ind];

			reduced_daei[i] = (aei1[ind] - aei0[ind]) / dt;

			reduced_doom[i] = oom1[ind] - oom0[ind];
			double dcur = reduced_doom[i].x + reduced_doom[i].y;
			double lam0 = oom0[ind].x + oom0[ind].y + oom0[ind].z;
			double lam1 = oom1[ind].x + oom1[ind].y + oom1[ind].z;

			// Omega
			if (reduced_doom[i].x > M_PI) reduced_doom[i].x -= 2 * M_PI;
			else if (reduced_doom[i].x < -M_PI) reduced_doom[i].x += 2 * M_PI;
			double Omefreq = reduced_doom[i].x / dt;

			// Curly pi (longitude of apo)
			if (dcur > M_PI) dcur -= 2 * M_PI;
			else if (dcur < -M_PI) dcur += 2 * M_PI;
			double cirfreq = dcur / dt;

			// guess the mean motion frequency, a must be in AU and t in days

			double mmfreq = std::sqrt(pl.m()[ind] + pl.m()[0]) 
							* (std::pow(aei0[ind].x, -1.5) + std::pow(aei1[ind].x, -1.5)) / 2;
			// mmfreq -= (reduced_doom[i].x + reduced_doom[i].y);
			// std::cout << ind << " predicted orbital period: " << 2 * M_PI / mmfreq << std::endl;
			// std::cout << ind << " freq: " << mmfreq << std::endl;
			// std::cout << ind << " dt: " << dt << std::endl;
			// std::cout << ind << " oom0[ind].z: " << oom0[ind].z << std::endl;
			// std::cout << ind << " oom1[ind].z: " << oom1[ind].z << std::endl;

			double cmfin = lam0 + mmfreq * dt;
			cmfin = std::fmod(cmfin, 2 * M_PI);
			double corr = lam1 - cmfin;
			// std::cout << ind << " corr: " << corr << std::endl;
			if (corr > M_PI) corr -= 2 * M_PI;
			else if (corr < -M_PI) corr += 2 * M_PI;

			// correct the frequency to match final mean motion

			mmfreq += corr / dt;

			// std::cout << ind << " corr: " << corr << std::endl;
			// std::cout << ind << " freq: " << mmfreq << std::endl;

			reduced_doom[i].x = Omefreq;
			reduced_doom[i].y = cirfreq - Omefreq;
			reduced_doom[i].z = mmfreq - cirfreq;
			i++;
		}
		pl_alive = i;
	}
}
}
