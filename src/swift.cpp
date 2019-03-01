#include "swift.h"
#include "convert.h"
#include <cstring>
#include <libgen.h>
#include <fstream>
#include <sys/wait.h>

namespace sr
{
namespace swift
{
	SwiftEncounterIntegrator::SwiftEncounterIntegrator(const sr::data::Configuration& config_, double t_, size_t prev_tbsize_, size_t cur_tbsize_) : t(t_), config(config_), prev_tbsize(prev_tbsize_), cur_tbsize(cur_tbsize_)
	{
	}

	void SwiftEncounterIntegrator::begin_integrate(const sr::data::HostPlanetPhaseSpace& pl, const sr::data::HostParticlePhaseSpace& pa)
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

		for (unsigned int i = 0; i < config.num_thread; i++)
		{
			size_t chunk_begin = encounter_start + n_encounter * i / config.num_thread;
			size_t chunk_end = encounter_start + n_encounter * (i + 1) / config.num_thread;

			std::string history_path = "/tmp/hist" + std::to_string(i);
			std::string param_path = "/tmp/param" + std::to_string(i);
			std::string tp_path = "/tmp/tp" + std::to_string(i);
			std::string pl_path = "/tmp/pl";

			// write_planetary_history(pl, history_path);
			write_param_in(param_path);
			write_tp_in(pa, chunk_begin, chunk_end, tp_path);
			write_pl_in(pl, pl_path);

			pid_t pid = fork();
			if (pid == (pid_t) -1)
			{
				throw std::runtime_error("fork error");
			}
			else if (pid == 0)
			{
				// child process
				execl(config.swift_path.c_str(), swift_basename.c_str(), param_path.c_str(), tp_path.c_str(), pl_path.c_str(), nullptr);
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
			waitpid(pid, nullptr, 0);
		}

		// TODO read particle data and insert back into pa
	}

	void SwiftEncounterIntegrator::write_param_in(std::string dest) const
	{
		std::ofstream file(dest);
		file << "0 " << t << " " << config.dt << std::endl;
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
			file << pl.m()[i] << std::endl;

			if (i == 0 && (pl.r()[i].lensq() != 0 || pl.v()[i].lensq() != 0))
			{
				throw reinterpret_cast<void*>(1);
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
		}
		
		file << "0 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
		file << "0 0 0 0 0" << std::endl;
		file << "0 0 0 0 0" << std::endl;
		file << "0 0 0" << std::endl;
	}

	void SwiftEncounterIntegrator::write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, double time, std::string dest) const
	{
		// write two (three?) planet locations: beginning, beginning, end?
		// need to handle planet mergers as a special case - may need to break on a merger

		std::ofstream file(dest, std::ios_base::binary);
		int npl = pl.n_alive_old();

		// TODO use the right time here
		sr::data::write_binary<double>(file, t - prev_tbsize * config.dt);
		sr::data::write_binary<uint32_t>(file, npl);

		for (size_t i = 1; i < npl; i++)
		{
			double a, e, I, O, o, m;
			sr::convert::to_elements(pl.m()[0], pl.r_log().get<false, true>()[pl.log_index_at<true>(0, i)],
					pl.v_log().get<false, true>()[pl.log_index_at<true>(0, i)],
					nullptr, &a, &e, &I, &O, &o, &m);
			sr::data::write_binary<uint32_t>(file, pl.id()[i]);
			sr::data::write_binary<float>(file, pl.m()[i]);
			sr::data::write_binary<float>(file, a);
			sr::data::write_binary<float>(file, e);
			sr::data::write_binary<float>(file, I);
			sr::data::write_binary<float>(file, O);
			sr::data::write_binary<float>(file, o);
			sr::data::write_binary<float>(file, m);
		}

		npl = pl.n_alive();
		// TODO use right time
		sr::data::write_binary<double>(file, t);
		sr::data::write_binary<uint32_t>(file, npl);
		for (size_t i = 1; i < npl; i++)
		{
			double a, e, I, O, o, m;
			sr::convert::to_elements(pl.m()[0], pl.r_log().get<false, false>()[pl.log_index_at<false>(0, i)],
					pl.v_log().get<false, false>()[pl.log_index_at<false>(0, i)],
					nullptr, &a, &e, &I, &O, &o, &m);
			sr::data::write_binary<uint32_t>(file, pl.id()[i]);
			sr::data::write_binary<float>(file, pl.m()[i]);
			sr::data::write_binary<float>(file, a);
			sr::data::write_binary<float>(file, e);
			sr::data::write_binary<float>(file, I);
			sr::data::write_binary<float>(file, O);
			sr::data::write_binary<float>(file, o);
			sr::data::write_binary<float>(file, m);
		}

		npl = pl.n_alive();
		sr::data::write_binary<double>(file, t + (cur_tbsize - 1) * config.tbsize);
		sr::data::write_binary<uint32_t>(file, npl);
		for (size_t i = 1; i < npl; i++)
		{
			double a, e, I, O, o, m;
			sr::convert::to_elements(pl.m()[0], pl.r_log().get<false, false>()[pl.log_index_at<false>(cur_tbsize - 1, i)],
					pl.v_log().get<false, false>()[pl.log_index_at<false>(cur_tbsize - 1, i)],
					nullptr, &a, &e, &I, &O, &o, &m);
			sr::data::write_binary<uint32_t>(file, pl.id()[i]);
			sr::data::write_binary<float>(file, pl.m()[i]);
			sr::data::write_binary<float>(file, a);
			sr::data::write_binary<float>(file, e);
			sr::data::write_binary<float>(file, I);
			sr::data::write_binary<float>(file, O);
			sr::data::write_binary<float>(file, o);
			sr::data::write_binary<float>(file, m);
		}
	}
}
}
