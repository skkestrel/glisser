#include "swift.h"
#include "convert.h"
#include "util.h"
#include "interp.h"

#include <iomanip>
#include <chrono>
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
	using hrclock = std::chrono::high_resolution_clock;

	SwiftEncounterIntegrator::SwiftEncounterIntegrator() { }

	SwiftEncounterIntegrator::SwiftEncounterIntegrator(
			const sr::data::Configuration& config,
			size_t npa)
	{
		swift_statlen = config.swift_statlen;
		swift_path = config.swift_path;
		outfolder = config.outfolder;
		swift_part_min = config.swift_part_min;
		num_swift = config.num_swift;
		inner_bound = config.inner_bound;
		outer_bound = config.outer_bound;

		for (unsigned int i = 0; i < swift_statlen; i++)
		{
			istat.push_back(std::vector<int32_t>(npa));
			rstat.push_back(std::vector<double>(npa));
		}

		// temp_log.open("temp_log_swift.txt");
	}

	std::pair<double, double> SwiftEncounterIntegrator::begin_integrate(
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
		auto s_clock = hrclock::now();
		ASSERT(_children.size() == 0, "")

		t = _t;
		rel_t = _rel_t;
		dt = _dt;
		prev_tbsize = _prev_tbsize;
		cur_tbsize = _cur_tbsize;
		old_log = old;

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

		
		double iotime = 0;

		auto clock = hrclock::now();
		write_planetary_history(pl, interp, history_path);

		std::chrono::duration<double, std::milli> millis = hrclock::now() - clock;
		iotime += millis.count();

		for (unsigned int i = 0; i < num_proc; i++)
		{
			size_t chunk_begin = encounter_start + n_encounter * i / num_proc;
			size_t chunk_end = encounter_start + n_encounter * (i + 1) / num_proc;

			std::string param_path = sr::util::joinpath(datapath, "param" + std::to_string(i));
			std::string tp_path = sr::util::joinpath(datapath, "tp" + std::to_string(i));
			std::string tpout_path = sr::util::joinpath(datapath, "tpout" + std::to_string(i));

			clock = hrclock::now();
			write_param_in(param_path);
			write_tp_in(pa, chunk_begin, chunk_end, tp_path);

			millis = hrclock::now() - clock;
			iotime += millis.count();

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

		millis = hrclock::now() - s_clock;
		double totaltime = millis.count();
		return std::make_pair(iotime, totaltime);
	}

	std::pair<double, double> SwiftEncounterIntegrator::end_integrate(sr::data::HostParticlePhaseSpace& pa)
	{
		double iotime = 0;
		double waittime = 0;

		if (_children.size() == 0)
		{
			throw std::runtime_error(".");
		}

		for (const ChildProcess& child : _children)
		{
			auto clock = hrclock::now();

			int status;
			::waitpid(child.pid, &status, 0);

			// construct buffer from file descriptor
			// __gnu_cxx::stdio_filebuf<char> filebuf(child.piper, std::ios::in);
			// std::istream is(&filebuf);
			::close(child.piper);

			std::chrono::duration<double, std::milli> millis = hrclock::now() - clock;
			waittime += millis.count();


			// std::cin.get();

			if (WEXITSTATUS(status) != 0)
			{
				throw std::runtime_error("swift did not terminate normally");
			}


			clock = hrclock::now();

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
					size_t statsiz = statmap.size();
					statmap[id] = statsiz;
					statindex = statsiz;
				}

				for (size_t j = 0; j < swift_statlen; j++)
				{
					int val;
					output >> val;
					if (j < 3)
					{
						istat[j][statindex] = val;
					}
					else
					{
						// total # of encounters with planet n
						istat[j][statindex] += val;
					}
				}

				for (size_t j = 0; j < swift_statlen; j++)
				{
					double real;
					output >> real;
					if (j < 3) rstat[j][statindex] = real;
					else
					{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
						if (rstat[j][statindex] == 0.0)
						{
							rstat[j][statindex] = real;
						}
						else if (real != 0.0)
						{
							// closest approach distance to a planet
							rstat[j][statindex] = std::min(rstat[j][statindex], real);
						}
#pragma GCC diagnostic pop
					}
				}

				uint16_t deathflags = 0x0000;
				// particle is dead

				if (istat[0][statindex] == 1)
				{
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
						// first planet is index 1 here, or index 2 in swift
						uint32_t planet_index = istat[1][statindex] - 1;

						deathflags = static_cast<uint16_t>(planet_id_list[planet_index] << 8) | 0x80;
					}
				}
				// particle is alive
				else
				{
					// still in encounter!
					if (istat[1][statindex] != 0)
					{
						uint32_t swift_planet_index = std::abs(istat[1][statindex]);
						// if we're in inner region, decrease the istat count of that planet by one since we are ending in an encounter
						if (istat[1][statindex] < 0)
						{
							istat[swift_planet_index + 1][statindex] -= 1;
						}

						// even if we're in encounter, don't udpate the particle encounter flag
						// instead, set the encounter flag to 0 and let the GPU detect that the particle is in encounter using hill radius
					}
					deathflags = 0;
				}

				pa.deathflags()[i + child.chunk_begin] = deathflags;
			}

			millis = hrclock::now() - clock;
			iotime += millis.count();
		}

		_children.clear();
		return std::make_pair(iotime, waittime);
	}

	void SwiftEncounterIntegrator::write_stat(std::string dest) const
	{
		std::ofstream file(dest);
		for (std::pair<uint32_t, size_t> pair : statmap)
		{
			file << std::fixed << std::setprecision(17) << pair.first << std::endl;
			for (size_t j = 0; j < swift_statlen; j++)
			{
				file << istat[j][pair.second] << " ";
			}
			file << std::endl;
			for (size_t j = 0; j < swift_statlen; j++)
			{
				file << rstat[j][pair.second] << " ";
			}
			file << std::endl;
		}
	}

	void SwiftEncounterIntegrator::write_param_in(std::string dest) const
	{
		std::ofstream file(dest);

		file << rel_t - static_cast<double>(prev_tbsize) * dt << " " << rel_t + static_cast<double>(cur_tbsize) * dt << " " << dt << std::endl;
		file << std::setprecision(17);
		file << "99999999 99999999" << std::endl;
		file << "F T F F T F" << std::endl;
		file << inner_bound << " " << outer_bound << " -1. -1. T" << std::endl;
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
			file << std::setprecision(17);
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
	// since the swift integrator never crosses the lookup boundary
	void SwiftEncounterIntegrator::write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, const sr::interp::Interpolator& interp, std::string dest)
	{
		std::ofstream file(dest, std::ios_base::binary);

		size_t binary_chunk_size = 60;
		sr::data::write_binary(file, static_cast<double>(pl.m()[0]));
		sr::data::pad_binary(file, binary_chunk_size - 8);
		// temp_log << std::setprecision(17);
		// temp_log << "solar mass: "<< pl.m()[0] << std::endl;

		if (old_log)
		{
			// temp_log << "old_log" << std::endl;
			planet_id_list = interp.reduced_ids_old;

			size_t npl = pl.n_alive_old();

			// interp n_alive doesn't include the sun
			ASSERT(interp.pl_alive_old + 1 == pl.n_alive_old(), "")

			sr::data::write_binary(file, static_cast<double>(0));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, binary_chunk_size - 8 - 4);

			// temp_log << "planets alive: "<< npl - 1 << std::endl;

			// reduced ones use a d
			for (size_t i = 0; i < npl - 1; i++)
			{
				// temp_log << "planet: " << interp.reduced_ids_old[i] + 1 << std::endl;
				// temp_log << interp.reduced_m_old[i] << std::endl;
				// temp_log << interp.jacobi_aei_i_old[i] << ' ' << interp.jacobi_oom_i_old[i] << std::endl;

				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids_old[i] + 1));
				sr::data::write_binary(file, static_cast<double>(interp.reduced_m_old[i]));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_i_old[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_i_old[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_i_old[i].z));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_i_old[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_i_old[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_i_old[i].z));
			}

			// temp_log << "t0: " << interp.t0 << " t_m1: " << interp.t_m1 << std::endl;
			sr::data::write_binary(file, static_cast<double>(interp.t0 - interp.t_m1));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, binary_chunk_size - 8 - 4);

			for (size_t i = 0; i < npl - 1; i++)
			{
				// temp_log << "planet: " << interp.reduced_ids_old[i] + 1 << std::endl;
				// temp_log << interp.reduced_m_old[i] << std::endl;
				// temp_log << interp.jacobi_aei_f_old[i] << ' ' << interp.jacobi_oom_f_old[i] << std::endl;

				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids_old[i] + 1));
				sr::data::write_binary(file, static_cast<double>(interp.reduced_m_old[i]));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_f_old[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_f_old[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_f_old[i].z));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_f_old[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_f_old[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_f_old[i].z));
			}

		}
		else
		{
			// temp_log << "new_log" << std::endl;
			planet_id_list = interp.reduced_ids;

			size_t npl = pl.n_alive();
			ASSERT(interp.pl_alive + 1 == pl.n_alive(), "")

			sr::data::write_binary(file, static_cast<double>(0));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, binary_chunk_size - 8 - 4);
			// temp_log << "planets alive: "<< npl - 1 << std::endl;

			for (size_t i = 0; i < npl - 1; i++)
			{
				// temp_log << "planet: " << interp.reduced_ids[i] + 1 << std::endl;
				// temp_log << interp.reduced_m[i] << std::endl;
				// temp_log << interp.jacobi_aei_i[i] << ' ' << interp.jacobi_oom_i[i] << std::endl;
				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids[i] + 1));
				sr::data::write_binary(file, static_cast<double>(interp.reduced_m[i]));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_i[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_i[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_i[i].z));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_i[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_i[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_i[i].z));
			}
			// temp_log << "t0: " << interp.t0 << " t_m1: " << interp.t_m1 << std::endl;
			sr::data::write_binary(file, static_cast<double>(interp.t1 - interp.t0));
			sr::data::write_binary(file, static_cast<uint32_t>(npl - 1));
			sr::data::pad_binary(file, binary_chunk_size - 8 - 4);

			for (size_t i = 0; i < npl - 1; i++)
			{
				// temp_log << "planet: " << interp.reduced_ids[i] + 1 << std::endl;
				// temp_log << interp.reduced_m[i] << std::endl;
				// temp_log << interp.jacobi_aei_f[i] << ' ' << interp.jacobi_oom_f[i] << std::endl;
				sr::data::write_binary(file, static_cast<uint32_t>(interp.reduced_ids[i] + 1));
				sr::data::write_binary(file, static_cast<double>(interp.reduced_m[i]));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_f[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_f[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_aei_f[i].z));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_f[i].x));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_f[i].y));
				sr::data::write_binary(file, static_cast<double>(interp.jacobi_oom_f[i].z));
			}
		}
	}
}
}
