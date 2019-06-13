#include "swift.h"
#include "convert.h"
#include "util.h"
#include "interp.h"

#include <iomanip>
#include <cstring>
#include <libgen.h>
#include <fstream>
#include <iostream>
#include <sys/wait.h>
#include <ext/stdio_filebuf.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

namespace sr
{
namespace swift
{
	SwiftEncounterIntegrator::SwiftEncounterIntegrator() { }

	SwiftEncounterIntegrator::SwiftEncounterIntegrator(
			const sr::data::Configuration& config,
			size_t npl)
	{
		swift_statlen = config.swift_statlen;
		swift_path = config.swift_path;
		outfolder = config.outfolder;
		swift_part_min = config.swift_part_min;
		num_swift = config.num_swift;
		outer_radius = config.outer_radius;

		for (unsigned int i = 0; i < swift_statlen; i++)
		{
			istat.push_back(std::vector<int32_t>(npl));
			rstat.push_back(std::vector<double>(npl));
		}
	}

	void SwiftEncounterIntegrator::begin_integrate(
			const sr::data::HostPlanetPhaseSpace& pl,
			const sr::data::HostParticlePhaseSpace& pa,
			const sr::interp::Interpolator& interp,
			bool old,
			double _t,
			double _rel_t,
			double _dt,
			size_t _prev_tbsize,
			size_t _cur_tbsize)
	{
		ASSERT(_children.size() == 0, "")

		t = _t;
		rel_t = _rel_t;
		dt = _dt;
		prev_tbsize = _prev_tbsize;
		cur_tbsize = _cur_tbsize;

		std::string datapath = sr::util::joinpath(outfolder, "swift_data");
		sr::util::make_dir(datapath);
		std::string swift_basename;

		{
			char copied_path[512];
			std::strcpy(copied_path, swift_path.c_str());
			swift_basename = std::string(basename(copied_path));
		}


		size_t encounter_start = pa.n_alive() - pa.n_encounter();
		size_t n_encounter = pa.n_encounter();

		// Set process count
		// For example with swift_part_min = 10:
		// 1 particle -> 1 process
		// 10 particles -> 1 process
		// 19 particles -> 1 process
		// 20 particles -> 2 process
		// 29 particles -> 2 process
		size_t num_proc = std::min(n_encounter / swift_part_min, static_cast<size_t>(num_swift));
		if (num_proc == 0) num_proc = 1;

		std::string history_path = sr::util::joinpath(datapath, "hist");

		write_planetary_history(pl, interp, history_path, old);

		for (unsigned int i = 0; i < num_proc; i++)
		{
			size_t chunk_begin = encounter_start + n_encounter * i / num_proc;
			size_t chunk_end = encounter_start + n_encounter * (i + 1) / num_proc;

			std::string param_path = sr::util::joinpath(datapath, "param" + std::to_string(i));
			std::string tp_path = sr::util::joinpath(datapath, "tp" + std::to_string(i));
			std::string tpout_path = sr::util::joinpath(datapath, "tpout" + std::to_string(i));

			write_param_in(param_path);
			write_tp_in(pa, chunk_begin, chunk_end, tp_path);

			int pipefd[2];

			if (::pipe(pipefd))
			{
				throw std::runtime_error("pipe error");
			}

			::pid_t pid = ::fork();
			if (pid == (::pid_t) -1)
			{
				throw std::runtime_error("fork error");
			}
			else if (pid == 0)
			{
				// child process
				// copy pipe write end into stdout
				::dup2(pipefd[1], 1);
				::close(pipefd[1]);
				::close(pipefd[0]);

				::execl(swift_path.c_str(), swift_basename.c_str(), param_path.c_str(), history_path.c_str(), tp_path.c_str(), tpout_path.c_str(), nullptr);
				::exit(-1);
			}


			// parent process
			::close(pipefd[1]);
			_children.push_back(ChildProcess(pid, pipefd[0], tpout_path, chunk_begin, chunk_end));
		}
	}

	void SwiftEncounterIntegrator::end_integrate(sr::data::HostParticlePhaseSpace& pa)
	{
		if (_children.size() == 0)
		{
			throw std::runtime_error(".");
		}

		for (const ChildProcess& child : _children)
		{
			int status;
			::waitpid(child.pid, &status, 0);

			// construct buffer from file descriptor
			__gnu_cxx::stdio_filebuf<char> filebuf(child.piper, std::ios::in);
			std::istream is(&filebuf);

			// std::cout << "encounter for (" << prev_tbsize << ", " << cur_tbsize << ")" << std::endl;
			// std::cout << std::string(std::istreambuf_iterator<char>(is), {}) << std::endl;


			// std::cin.get();


			if (WEXITSTATUS(status) != 0)
			{
				throw std::runtime_error("swift did not terminate normally");
			}

			std::ifstream output(child.tpout);

			size_t npa;

			output >> npa;

			ASSERT(npa == child.chunk_end - child.chunk_begin, "")

			for (size_t i = 0; i < npa; i++)
			{
				double x, y, z, vx, vy, vz;

				output >> x >> y >> z >> vx >> vy >> vz;

				pa.r()[i + child.chunk_begin] = f64_3(x, y, z);
				pa.v()[i + child.chunk_begin] = f64_3(vx, vy, vz);

				uint32_t id = pa.id()[i + child.chunk_begin];
				size_t statindex;
				if (statmap.count(id) > 0)
				{
					statindex = statmap[id];
				}
				else
				{
					statindex = statmap[id] = statmap.size();
				}

				for (size_t j = 0; j < swift_statlen; j++)
				{
					output >> istat[j][statindex];
				}

				for (size_t j = 0; j < swift_statlen; j++)
				{
					output >> rstat[j][statindex];
				}

				uint16_t deathflags = 0x0000;
				// particle is dead

				if (istat[0][statindex] == 1)
				{
					// std::cout << "die" << std::endl;

					if (istat[1][statindex] == -1)
					{
						deathflags = 0x04;
					}
					else if (istat[1][statindex] == -2 || istat[1][statindex] == -3)
					{
						deathflags = 0x02;
					}
					else if (istat[1][statindex] == -4 || istat[1][statindex] == 1)
					{
						deathflags = 0x80;
					}
					else
					{
						// first planet is index 0 here, or index 2 in swift
						uint32_t planet_index = istat[1][statindex] - 2;

						deathflags = static_cast<uint16_t>(planet_id_list[planet_index] << 8) | 0x80;
					}
				}
				// particle is alive
				else
				{
					/*
						

					// still in encounter!
					if (istat[1][statindex] != 0)
					{
						uint32_t swift_planet_index = std::abs(istat[1][statindex]);

						// sun
						if (swift_planet_index == 1)
						{
							deathflags = 0x01;
						}
						else
						{
							uint32_t planet_index = swift_planet_index - 2;
							deathflags = static_cast<uint16_t>(planet_id_list[planet_index] << 8) | 0x01;
						}
					}
					else
					{
						// clear flags
						deathflags = 0;
					}

					*/


					// note: if we try to read the encounter flag from swift and replicate it here, the particle is going
					// to get handled on the very next encounter timechunk, which we don't want

					// isntead, set the encounter flag to 0 and let the GPU detect that the particle is in encounter using hill radius
					deathflags = 0;
				}

				// std::cout << istat[1][statindex] << std::endl;
				// std::cout << rstat[3][statindex] << " " << rstat[4][statindex] << std::endl;

				pa.deathflags()[i + child.chunk_begin] = deathflags;
			}
		}

		_children.clear();
	}

	void SwiftEncounterIntegrator::write_param_in(std::string dest) const
	{
		std::ofstream file(dest);

		file << rel_t - static_cast<double>(prev_tbsize) * dt << " " << rel_t + static_cast<double>(cur_tbsize) * dt << " " << dt << std::endl;
		file << std::setprecision(17) << std::endl;
		file << "999999 9999999" << std::endl;
		file << "F T F F T F" << std::endl;
		file << "0.5 " << outer_radius << " -1. -1. T" << std::endl;
		file << "/dev/null" << std::endl;
		file << "unknown" << std::endl;
		file << std::flush;
	}

	void SwiftEncounterIntegrator::write_tp_in(const sr::data::HostParticlePhaseSpace& pa, size_t chunk_begin, size_t chunk_end, std::string dest) const
	{
		std::ofstream file(dest);
		file << chunk_end - chunk_begin << std::endl;

		for (size_t i = chunk_begin; i < chunk_end; i++)
		{
			file << std::setprecision(17) << std::endl;
			file << pa.r()[i] << std::endl;
			file << pa.v()[i] << std::endl;
		
			for (size_t j = 0; j < swift_statlen; j++)
			{
				file << "0 ";
			}
			file << std::endl;
			for (size_t j = 0; j < swift_statlen; j++)
			{
				file << "0.0 ";
			}
			file << std::endl;
		}
	}

	// this function handles both two-step integrations, and initial one-step integrations
	// no matter whether the integration is one or two step, the lookup is always bounded by interp.t0 and interp.t1
	// since the swift integratio never crosses the lookup boundary
	void SwiftEncounterIntegrator::write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, const sr::interp::Interpolator& interp, std::string dest, bool old)
	{
		std::ofstream file(dest, std::ios_base::binary);

		sr::data::write_binary(file, static_cast<double>(pl.m()[0]));
		sr::data::pad_binary(file, 32 - 8);

		if (old)
		{
			planet_id_list = interp.reduced_ids_old;

			size_t npl = pl.n_alive_old();

			// interp n_alive doesn't include the sun
			ASSERT(interp.n_alive_old + 1 == pl.n_alive_old(), "")

			sr::data::write_binary(file, static_cast<double>(0));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, 32 - 8 - 4);

			// reduced ones use a d
			for (size_t i = 0; i < npl - 1; i++)
			{
				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids_old[i] + 1));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_m_old[i]));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_i_old[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_i_old[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_i_old[i].z));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_i_old[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_i_old[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_i_old[i].z));
			}

			sr::data::write_binary(file, static_cast<double>(interp.t0 - interp.t_m1));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, 32 - 8 - 4);

			for (size_t i = 0; i < npl - 1; i++)
			{
				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids_old[i] + 1));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_m_old[i]));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_f_old[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_f_old[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_f_old[i].z));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_f_old[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_f_old[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_f_old[i].z));
			}

		}
		else
		{
			planet_id_list = interp.reduced_ids;

			size_t npl = pl.n_alive();
			ASSERT(interp.n_alive + 1 == pl.n_alive(), "")

			sr::data::write_binary(file, static_cast<double>(0));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, 32 - 8 - 4);

			for (size_t i = 0; i < npl - 1; i++)
			{
				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids[i] + 1));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_m[i]));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_i[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_i[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_i[i].z));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_i[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_i[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_i[i].z));
			}

			sr::data::write_binary(file, static_cast<double>(interp.t1 - interp.t0));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, 32 - 8 - 4);

			for (size_t i = 0; i < npl - 1; i++)
			{
				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids[i] + 1));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_m[i]));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_f[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_f[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_aei_f[i].z));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_f[i].x));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_f[i].y));
				sr::data::write_binary(file, static_cast<float>(interp.reduced_oom_f[i].z));
			}
		}
	}
}
}
