#include "swift.h"
#include "convert.h"
#include "util.h"
#include "interp.h"
#include <cstring>
#include <libgen.h>
#include <fstream>
#include <sys/wait.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

namespace sr
{
namespace swift
{
	SwiftEncounterIntegrator::SwiftEncounterIntegrator(const sr::data::Configuration& config_, double t_, double dt_, size_t prev_tbsize_, size_t cur_tbsize_) : t(t_), dt(dt_), prev_tbsize(prev_tbsize_), cur_tbsize(cur_tbsize_), config(config_)
	{
	}

	void SwiftEncounterIntegrator::begin_integrate(const sr::data::HostPlanetPhaseSpace& pl, const sr::data::HostParticlePhaseSpace& pa, const sr::interp::Interpolator& interp, bool old)
	{
		if (_children.size() > 0)
		{
			throw std::runtime_error(".");
		}

		std::string swift_basename;

		{
			char copied_path[512];
			std::strcpy(copied_path, config.swift_path.c_str());
			swift_basename = std::string(basename(copied_path));
		}


		size_t encounter_start = pa.n_alive() - pa.n_encounter();
		size_t n_encounter = pa.n_encounter();

		// For example with swift_part_min = 10:
		// 1 particle -> 1 process
		// 10 particles -> 1 process
		// 19 particles -> 1 process
		// 20 particles -> 2 process
		// 29 particles -> 2 process
		size_t num_proc = std::min(n_encounter / config.swift_part_min, static_cast<size_t>(config.num_swift));
		if (num_proc == 0) num_proc = 1;

		for (unsigned int i = 0; i < num_proc; i++)
		{
			size_t chunk_begin = encounter_start + n_encounter * i / num_proc;
			size_t chunk_end = encounter_start + n_encounter * (i + 1) / num_proc;

			std::string mypid = std::to_string(::getpid());

			sr::util::make_dir("/tmp/glisse");

			std::string history_path = "/tmp/glisse/" + mypid + "hist" + std::to_string(i);
			std::string param_path = "/tmp/glisse/" + mypid + "param" + std::to_string(i);
			std::string tp_path = "/tmp/glisse/" + mypid + "tp" + std::to_string(i);
			std::string pl_path = "/tmp/glisse/" + mypid + "pl";

			write_planetary_history(pl, interp, history_path, old);
			write_param_in(param_path);
			write_tp_in(pa, chunk_begin, chunk_end, tp_path);
			write_pl_in(pl, pl_path);

			::pid_t pid = ::fork();
			if (pid == (pid_t) -1)
			{
				throw std::runtime_error("fork error");
			}
			else if (pid == 0)
			{
				// child process
				int null = ::open("/dev/null", O_WRONLY);
				::dup2(null, 1);
				::execl(config.swift_path.c_str(), swift_basename.c_str(), param_path.c_str(), pl_path.c_str(), tp_path.c_str(), nullptr);
				exit(-1);
			}
			else
			{
				// parent process
				_children.push_back(pid);
			}
		}
	}

	void SwiftEncounterIntegrator::end_integrate(sr::data::HostParticlePhaseSpace& pa)
	{
		if (_children.size() == 0)
		{
			throw std::runtime_error(".");
		}

		for (const pid_t& pid : _children)
		{
			int status;
			::waitpid(pid, &status, 0);
			if (WEXITSTATUS(status) != 0)
			{
				throw std::runtime_error("swift did not terminate normally");
			}
		}

		
		
		// TODO read particle data and insert back into pa
	}

	void SwiftEncounterIntegrator::write_param_in(std::string dest) const
	{
		std::ofstream file(dest);
		file << "0 " << static_cast<double>(prev_tbsize + cur_tbsize) * dt << " " << dt << std::endl;
		file << "999999 9999999" << std::endl;
		file << "F T F F T F" << std::endl;
		file << "0.5 500 200 -1 T" << std::endl;
		file << "/dev/null" << std::endl;
		file << "unknown" << std::endl;
	}

	void SwiftEncounterIntegrator::write_pl_in(const sr::data::HostPlanetPhaseSpace& pl, std::string dest) const
	{
		std::ofstream file(dest);
		file << pl.n_alive() << std::endl;

		for (size_t i = 0; i < pl.n_alive(); i++)
		{
			if (i == 0)
			{
				file << pl.m()[i] << std::endl;
			}
			else
			{
				file << pl.m()[i] << " " << config.cull_radius << std::endl;
			}

			if (i == 0 && (pl.r()[i].lensq() > 1e-16 || pl.v()[i].lensq() > 1e-16))
			{
				throw std::runtime_error("not in heliocentric coordinates");
			}

			file << pl.r()[i] << std::endl;
			file << pl.v()[i] << std::endl;
		}
	}

	void SwiftEncounterIntegrator::write_tp_in(const sr::data::HostParticlePhaseSpace& pa, size_t chunk_begin, size_t chunk_end, std::string dest) const
	{
		std::ofstream file(dest);
		file << chunk_end - chunk_begin << std::endl;
		for (size_t i = chunk_begin; i < chunk_end; i++)
		{
			file << pa.r()[i] << std::endl;
			file << pa.v()[i] << std::endl;
		
			for (size_t j = 0; j < config.swift_statlen; j++)
			{
				file << "0 ";
			}
			file << std::endl;
			for (size_t j = 0; j < config.swift_statlen; j++)
			{
				file << "0.0 ";
			}
			file << std::endl;
		}
	}

	// this function handles both two-step integrations, and initial one-step integrations
	// no matter whether the integration is one or two step, the lookup is always bounded by interp.t0 and interp.t1
	// since the swift integratio never crosses the lookup boundary
	void SwiftEncounterIntegrator::write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, const sr::interp::Interpolator& interp, std::string dest, bool old) const
	{
		std::ofstream file(dest, std::ios_base::binary);
		// TODO how to deal with changing planet number?

		size_t npl = old ? pl.n_alive_old() : pl.n_alive();

		if (old)
		{
			sr::data::write_binary(file, static_cast<double>(interp.t_m1));
			sr::data::write_binary(file, static_cast<uint32_t>(npl));

			for (size_t i = 1; i < npl; i++)
			{
				sr::data::write_binary(file, static_cast<uint32_t>(pl.id()[i]));
				sr::data::write_binary(file, static_cast<float>(pl.m()[i]));
				sr::data::write_binary(file, static_cast<float>(interp.aei_m1[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.aei_m1[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.aei_m1[i].z));
				sr::data::write_binary(file, static_cast<float>(interp.oom_m1[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.oom_m1[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.oom_m1[i].z));
			}
		}

		sr::data::write_binary(file, static_cast<double>(interp.t0));
		sr::data::write_binary(file, static_cast<uint32_t>(npl));

		for (size_t i = 1; i < npl; i++)
		{
			sr::data::write_binary(file, static_cast<uint32_t>(pl.id()[i]));
			sr::data::write_binary(file, static_cast<float>(pl.m()[i]));
			sr::data::write_binary(file, static_cast<float>(interp.aei0[i].x));
			sr::data::write_binary(file, static_cast<float>(interp.aei0[i].y));
			sr::data::write_binary(file, static_cast<float>(interp.aei0[i].z));
			sr::data::write_binary(file, static_cast<float>(interp.oom0[i].x));
			sr::data::write_binary(file, static_cast<float>(interp.oom0[i].y));
			sr::data::write_binary(file, static_cast<float>(interp.oom0[i].z));
		}

		if (!old)
		{
			sr::data::write_binary(file, static_cast<double>(interp.t1));
			sr::data::write_binary(file, static_cast<uint32_t>(npl));
			for (size_t i = 1; i < npl; i++)
			{
				sr::data::write_binary(file, static_cast<uint32_t>(pl.id()[i]));
				sr::data::write_binary(file, static_cast<float>(pl.m()[i]));
				sr::data::write_binary(file, static_cast<float>(interp.aei1[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.aei1[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.aei1[i].z));
				sr::data::write_binary(file, static_cast<float>(interp.oom1[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.oom1[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.oom1[i].z));
			}
		}
	}

	void SwiftEncounterIntegrator::write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, std::string dest) const
	{
		// write three planet locations to interpolate: beginning of previous chunk, beginning of current chunk, end of current chunk?

		std::ofstream file(dest, std::ios_base::binary);
		size_t npl = pl.n_alive();

		if (prev_tbsize != 0)
		{
			// TODO need to get the planet position at t - prev_tbsize * dt, as the earliest time that the log holds is t - prev_tbsize * dt + dt
			sr::data::write_binary(file, t - static_cast<double>(prev_tbsize - 1) * dt);
			sr::data::write_binary(file, static_cast<uint32_t>(npl));

			for (size_t i = 1; i < npl; i++)
			{
				double a, e, I, O, o, m;
				sr::convert::to_elements(pl.m()[0], pl.r_log().get<true>()[pl.log_index_at<true>(0, i)],
						pl.v_log().get<true>()[pl.log_index_at<true>(0, i)],
						nullptr, &a, &e, &I, &O, &o, &m);
				sr::data::write_binary(file, static_cast<uint32_t>(pl.id()[i]));
				sr::data::write_binary(file, static_cast<float>(pl.m()[i]));
				sr::data::write_binary(file, static_cast<float>(a));
				sr::data::write_binary(file, static_cast<float>(e));
				sr::data::write_binary(file, static_cast<float>(I));
				sr::data::write_binary(file, static_cast<float>(O));
				sr::data::write_binary(file, static_cast<float>(o));
				sr::data::write_binary(file, static_cast<float>(m));
			}
		}

		sr::data::write_binary(file, static_cast<double>(t));
		sr::data::write_binary(file, static_cast<uint32_t>(npl));
		for (size_t i = 1; i < npl; i++)
		{
			double a, e, I, O, o, m;
			sr::convert::to_elements(pl.m()[0], pl.r_log().get<true>()[pl.log_index_at<true>(prev_tbsize - 1, i)],
					pl.v_log().get<true>()[pl.log_index_at<true>(prev_tbsize - 1, i)],
					nullptr, &a, &e, &I, &O, &o, &m);
			sr::data::write_binary(file, static_cast<uint32_t>(pl.id()[i]));
			sr::data::write_binary(file, static_cast<float>(pl.m()[i]));
			sr::data::write_binary(file, static_cast<float>(a));
			sr::data::write_binary(file, static_cast<float>(e));
			sr::data::write_binary(file, static_cast<float>(I));
			sr::data::write_binary(file, static_cast<float>(O));
			sr::data::write_binary(file, static_cast<float>(o));
			sr::data::write_binary(file, static_cast<float>(m));
		}

		npl = pl.n_alive();
		sr::data::write_binary(file, t + static_cast<double>(cur_tbsize) * dt);
		sr::data::write_binary(file, static_cast<uint32_t>(npl));
		for (size_t i = 1; i < npl; i++)
		{
			double a, e, I, O, o, m;
			sr::convert::to_elements(pl.m()[0], pl.r_log().get<false>()[pl.log_index_at<false>(cur_tbsize - 1, i)],
					pl.v_log().get<false>()[pl.log_index_at<false>(cur_tbsize - 1, i)],
					nullptr, &a, &e, &I, &O, &o, &m);
			sr::data::write_binary(file, static_cast<uint32_t>(pl.id()[i]));
			sr::data::write_binary(file, static_cast<float>(pl.m()[i]));
			sr::data::write_binary(file, static_cast<float>(a));
			sr::data::write_binary(file, static_cast<float>(e));
			sr::data::write_binary(file, static_cast<float>(I));
			sr::data::write_binary(file, static_cast<float>(O));
			sr::data::write_binary(file, static_cast<float>(o));
			sr::data::write_binary(file, static_cast<float>(m));
		}
	}
}
}
