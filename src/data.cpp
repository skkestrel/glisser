#include "data.h"
#include "util.h"
#include "convert.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace sr
{
namespace data
{
	// particles having close encounters with any planets (sun excluded)
	size_t stable_partition_alive_indices(const Vu16& flags, size_t begin, size_t length, std::unique_ptr<Vs>* indices)
	{
		auto new_indices = std::make_unique<Vs>(length);
		std::iota(new_indices->begin(), new_indices->end(), 0);

		auto val = std::stable_partition(new_indices->begin(), new_indices->end(), [begin, &flags](size_t index)
				{ return ((flags[index + begin] & 0x00FE) == 0) && (flags[index + begin] != 0x0001); }) - new_indices->begin() + begin;

		*indices = std::move(new_indices);
		return val;
	}

	// particles with non-zero low byte are unflagged.
	size_t stable_partition_unflagged_indices(const Vu16& flags, size_t begin, size_t length, std::unique_ptr<Vs>* indices)
	{
		auto new_indices = std::make_unique<Vs>(length);
		std::iota(new_indices->begin(), new_indices->end(), 0);

		auto val = std::stable_partition(new_indices->begin(), new_indices->end(), [begin, &flags](size_t index)
				{ return flags[index + begin] == 0; }) - new_indices->begin() + begin;

		*indices = std::move(new_indices);
		return val;
	}

	void HostParticleSnapshot::gather(const Vs& indices, size_t begin, size_t length)
	{
		sr::data::gather(r, indices, begin, length);
		sr::data::gather(v, indices, begin, length);
		sr::data::gather(id, indices, begin, length);
	}

	void HostParticleSnapshot::resize(size_t length)
	{
		r.resize(length);
		v.resize(length);
		id.resize(length);
		n_alive = std::min(length, n_alive);
		n = length;
	}

	void HostParticleSnapshot::filter(const Vs& filter, HostParticleSnapshot& out) const
	{
		out = HostParticleSnapshot(filter.size());
		out.n_alive = 0;

		size_t index = 0;

		for (size_t i : filter)
		{
			out.id[index] = id[i];
			out.r[index] = r[i];
			out.v[index] = v[i];

			if (i < n_alive)
			{
				out.n_alive++;
			}

			index++;
		}
	}

	void HostParticlePhaseSpace::gather(const Vs& indices, size_t begin, size_t length)
	{
		base.gather(indices, begin, length);
		sr::data::gather(deathflags(), indices, begin, length);
		sr::data::gather(deathtime_index(), indices, begin, length);
	}

	void HostParticlePhaseSpace::filter(const Vs& filter, HostParticlePhaseSpace& out) const
	{
		base.filter(filter, out.base);
		out._n_encounter = 0;
		out._deathflags = Vu16(out.base.n);
		out._deathtime = Mu32f64();
		out._deathtime_index = Vu32(out.base.n);

		size_t index = 0;

		for (size_t i : filter)
		{
			out._deathtime[base.id[i]] = _deathtime.at(base.id[i]);
			out._deathflags[index] = _deathflags[i];
			out._deathtime_index[index] = _deathtime_index[i];

			index++;
		}
	}

	std::unique_ptr<Vs> HostParticlePhaseSpace::stable_partition_unflagged(size_t begin, size_t length)
	{
		std::unique_ptr<Vs> indices;
		n_alive() = stable_partition_unflagged_indices(deathflags(), begin, length, &indices);
		
		this->gather(*indices, begin, length);

		return indices;
	}

	std::unique_ptr<Vs> HostParticlePhaseSpace::stable_partition_alive(size_t begin, size_t length)
	{
		std::unique_ptr<Vs> indices;
		n_alive() = stable_partition_alive_indices(deathflags(), begin, length, &indices);
		
		this->gather(*indices, begin, length);

		return indices;
	}

	std::unique_ptr<Vs> HostParticleSnapshot::sort_by_id(size_t begin, size_t length)
	{
		auto new_indices = std::make_unique<Vs>(length);

		std::iota(new_indices->begin(), new_indices->end(), 0);
		std::sort(new_indices->begin(), new_indices->end(), [begin, this](size_t a, size_t b)
				{ return id[a + begin] < id[b + begin]; });

		this->gather(*new_indices, 0, length);

		return new_indices;
	}

	std::unique_ptr<Vs> HostParticlePhaseSpace::sort_by_id(size_t begin, size_t length)
	{
		auto indices = base.sort_by_id(begin, length);

		this->gather(*indices, 0, length);
		return indices;
	}

	Configuration::Configuration()
	{
		t_0 = 0;
		t_f = 1e4;
		dt = 1;
		tbsize = 1024;
		write_bary_track = true;
		write_ele_track = true;
		hill_factor = 3;
		inner_bound = 0.005;
		interp_maxpl = 9;
		max_kep = 10;

		outer_bound = 3000;
		resync_every = 1;
		planet_hist_every = 0;
		swift_statlen = 13;
		log_every = 10;
		track_every = 0;
		split_track_file = 0;
		num_swift = 1;
		swift_part_min = 10;

		interp_planets = false;
		use_jacobi_interp = true;
		planet_history_file = "";

		dump_every = 1000;
		max_particle = static_cast<uint32_t>(-1);
		resolve_encounters = false;
		big_g = 1;
		icsin = "";
		plin = "";
		hybridin = "state.in";
		hybridout = "state.out";
		readsplit = 0;
		writesplit = 0;
		writebinary = 0;
		readbinary = 0;
		outfolder = "output/";
		swift_path = "";

		read_single_hist = true;
		read_binary_hist = true;
		write_single_hist = true;
		write_binary_hist = true;

		write_encounter_log = false;
		diagnostic_mode = false;
	}

	Configuration Configuration::output_config() const
	{
// #pragma GCC warning "TODO"
		Configuration new_conf = *this;
		return new_conf;
	}

	void read_configuration(std::istream& in, Configuration* out)
	{
		size_t linenum = 0;
		std::string line;
		while (std::getline(in, line))
		{
			linenum++;
			if (line.length() == 0) continue;
			if (line[0] == '#') continue;
			while(line[0] == ' ' or line[0] == '\t')
			{
				line.erase(0, 1);
			}

			size_t split = line.find(' ');
			if (split == std::string::npos)
			{
				std::ostringstream ss;
				ss << "Unrecognized line " << linenum << ": " << line;
				throw std::runtime_error(ss.str());
			}

			std::string first = line.substr(0, split);
			std::string second = line.substr(split + 1, line.length() - split - 1);

			try
			{
				// Directories
				if (first == "Input-File")
					out->hybridin = second;
				else if (first == "Output-Folder")
					out->outfolder = second;
				else if (first == "Output-File")
					out->hybridout = second;
				else if (first == "Read-Binary-Input")
					out->readbinary = std::stoi(second) != 0;
				else if (first == "Write-Binary-Output")
					out->writebinary = std::stoi(second) != 0;

				// Read Planet History
				else if (first == "Read-Planet-History")
					out->interp_planets = std::stoi(second) != 0;
					// If "Read-Planet-History" is set to 0, the following parameters will be ignored.
					else if (first == "Planet-History-File")
						out->planet_history_file = second;
					else if (first == "Planet-History-Max-Body-Count")
						out->interp_maxpl = std::stou(second);
					else if (first == "Read-Binary-History")
						out->read_binary_hist = std::stoi(second) != 0;
					else if (first == "Read-Single-Precision-History")
						out->read_single_hist = std::stoi(second) != 0;
					// else if (first == "Use-Jacobi-Interpolation")
					// 	out->use_jacobi_interp = std::stoi(second) != 0;

				// Time
				else if (first == "Initial-Time")
					out->t_0 = std::stod(second);
				else if (first == "Final-Time")
					out->t_f = std::stod(second);
				else if (first == "Time-Step")
					out->dt = std::stod(second);
				else if (first == "Time-Block-Size")
					out->tbsize = std::stou(second);

				// Write Intervals
				else if (first == "Resync-Interval")
					out->resync_every = std::stou(second);	
				else if (first == "Log-Interval")
					out->log_every = std::stou(second);
				else if (first == "Dump-Interval")
					out->dump_every = std::stou(second);
				else if (first == "Track-Interval")
					out->track_every = std::stou(second);
					// If "Track-Interval" is set to 0, the following parameters will be ignored.
					else if (first == "Write-Single-Precision-Track")
						out->write_single_track = std::stoi(second) != 0;
					else if (first == "Write-Elements-Track")
						out->write_ele_track = std::stoi(second) != 0;
					else if (first == "Write-Bary-Track")
						out->write_bary_track = std::stoi(second) != 0;
					else if (first == "Split-Track-File")
					out->split_track_file = std::stou(second);
				else if (first == "Planet-History-Interval")
					out->planet_hist_every = std::stou(second);
					// If "Planet-History-Interval" is set to 0, the following parameters will be ignored.
					else if (first == "Write-Binary-History")
						out->write_binary_hist = std::stoi(second) != 0;
					else if (first == "Write-Single-Precision-History")
						out->write_single_hist = std::stoi(second) != 0;
				else if (first == "Write-Encounter-Log")
					out->write_encounter_log = std::stoi(second) != 0;
		
				// Close Encounters
				else if (first == "Resolve-Encounters")
					out->resolve_encounters = std::stoi(second) != 0;
					// If "Resolve-Encounters" is set to 0, the following parameters will be ignored.
					else if (first == "Swift-Path")
						out->swift_path = second;
					else if (first == "Swift-Process-Count")
						out->num_swift = std::stou(second);
					else if (first == "Swift-Process-Min-Particle-Count")
						out->swift_part_min = std::stou(second);
					else if (first == "Swift-Status-Length")
						out->swift_statlen = std::stou(second);

				// Other Parameters
				else if (first == "Big-G")
					out->big_g = std::stod(second);
				else if (first == "Max-Kepler-Iterations")
					out->max_kep = std::stou(second);
				else if (first == "Hill-Radius-Factor")
					out->hill_factor = std::stod(second);
				else if (first == "Particle-Inner-Boundary")
					out->inner_bound = std::stod(second);
				else if (first == "Particle-Outer-Boundary")
					out->outer_bound = std::stod(second);
				else if (first == "Particle-Count-Limit")
					out->max_particle = std::stou(second);
				else if (first == "Diagnositc-Mode")
					out->diagnostic_mode = std::stoi(second) != 0;
				else
					throw std::invalid_argument("bad");
			}
			catch (std::invalid_argument)
			{
				std::ostringstream ss;
				ss << "Could not parse line " << linenum << ": " << line;
				throw std::runtime_error(ss.str());
			}
		}

		if (out->outfolder.empty())
		{
			std::cerr << "Warning: Output-Folder was not specified, defaulting to ./" << std::endl;
			out->outfolder = "./";
		}
		if (out->resync_every != 1 && out->resolve_encounters)
		{
			throw std::runtime_error("Error: Resync-Interval must be 1 when Resolve-Encounters is enabled");
		}
		if (!out->interp_planets && out->resolve_encounters)
		{
			throw std::runtime_error("Error: Cannot resolve encounters if not using an interpolated planetary history file");
		}
	}

	void write_configuration(std::ostream& outstream, const Configuration& out)
	{
		outstream << std::setprecision(17);
		outstream << "# Directories" << std::endl;
		if (out.hybridin != "state.in") outstream << "Input-File " << out.hybridin << std::endl;
		outstream << "Output-Folder " << out.outfolder << std::endl;
		if (out.hybridout != "state.out") outstream << "Output-File " << out.hybridout << std::endl;

		if (out.readbinary) outstream << "Read-Binary-Input " << out.readbinary << std::endl;
		if (out.writebinary) outstream << "Write-Binary-Output " << out.writebinary << std::endl;
		outstream << std::endl;
		
		if (out.interp_planets) 
		{
			outstream << "# Read Planet History" << std::endl;
			outstream << "Read-Planet-History " << out.interp_planets << std::endl;
			outstream << "    Planet-History-File " << out.planet_history_file << std::endl;
			outstream << "    Planet-History-Max-Body-Count " << out.interp_maxpl << std::endl;
			outstream << "    Read-Binary-History " << out.read_binary_hist << std::endl;
			outstream << "    Read-Single-Precision-History " << out.read_single_hist << std::endl;
			// outstream << "    Use-Jacobi-Interpolation " << out.use_jacobi_interp << std::endl;
			outstream << std::endl;
		}
		
		outstream << "# Time" << std::endl;
		outstream << "Initial-Time " << out.t_0 << std::endl;
		outstream << "Final-Time " << out.t_f << std::endl;
		outstream << "Time-Step " << out.dt << std::endl;
		outstream << "Time-Block-Size " << out.tbsize << std::endl;
		outstream << std::endl;
		
		outstream << "# Write Intervals" << std::endl;
		outstream << "Resync-Interval " << out.resync_every << std::endl;
		outstream << "Log-Interval " << out.log_every << std::endl;
		outstream << "Dump-Interval " << out.dump_every << std::endl;
		outstream << "Track-Interval " << out.track_every << std::endl;
		if (out.track_every) 
		{
			outstream << "    Write-Single-Precision-Track " << out.write_single_track << std::endl;
			outstream << "    Write-Elements-Track " << out.write_ele_track << std::endl;
			outstream << "    Write-Bary-Track " << out.write_bary_track << std::endl;
			if (out.split_track_file) outstream << "    Split-Track-File " << out.split_track_file << std::endl;
		}
		if (out.planet_hist_every) 
		{
			outstream << "Planet-History-Interval " << out.planet_hist_every << std::endl;
			outstream << "    Write-Binary-History " << out.write_binary_hist << std::endl;
			outstream << "    Write-Single-Precision-History " << out.write_single_hist << std::endl;
		}
		outstream << "Write-Encounter-Log " << out.write_encounter_log << std::endl;
		outstream << std::endl;

		if (out.resolve_encounters) 
		{
			outstream << "# Close Encounters" << std::endl;
			outstream << "Resolve-Encounters " << out.resolve_encounters << std::endl;
			outstream << "    Swift-Path " << out.swift_path << std::endl;
			outstream << "    Swift-Process-Count " << out.num_swift << std::endl;
			outstream << "    Swift-Process-Min-Particle-Count " << out.swift_part_min << std::endl;
			outstream << "    Swift-Status-Length " << out.swift_statlen << std::endl;
			outstream << std::endl;
		}

		outstream << "# Other Parameters" << std::endl;
		if (std::fabs(out.big_g - 1.0) > 1e-14) outstream << "Big-G " << out.big_g << std::endl;
		outstream << "Max-Kepler-Iterations " << out.max_kep << std::endl;
		outstream << "Hill-Radius-Factor " << out.hill_factor << std::endl;
		outstream << "Particle-Inner-Boundary " << out.inner_bound << std::endl;
		outstream << "Particle-Outer-Boundary " << out.outer_bound << std::endl;
		outstream << "Particle-Count-Limit " << out.max_particle << std::endl;		
		if (out.diagnostic_mode) outstream << "Diagnositc-Mode " << out.diagnostic_mode << std::endl;	
	}	

	/* 	
	bool load_planet_data(HostPlanetPhaseSpace& pl, const Configuration& config, std::istream& plin)
	{
		size_t npl;
		plin >> npl;

		pl = HostPlanetPhaseSpace(npl, config.tbsize);

		for (size_t i = 0; i < npl; i++)
		{
			plin >> pl.m()[i];
			plin >> pl.r()[i].x >> pl.r()[i].y >> pl.r()[i].z;
			plin >> pl.v()[i].x >> pl.v()[i].y >> pl.v()[i].z;

			pl.id()[i] = static_cast<uint32_t>(i);
		}

		return false;
	}

	bool load_data_nohybrid(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config, std::istream& plin, std::istream& icsin)
	{
		load_planet_data(pl, config, plin);

		size_t npart;
		icsin >> npart;
		npart = std::min(npart, static_cast<size_t>(config.max_particle));

		pa = HostParticlePhaseSpace(npart);

		for (uint32_t i = 0; i < npart; i++)
		{
			icsin >> pa.r()[i].x >> pa.r()[i].y >> pa.r()[i].z;
			icsin >> pa.v()[i].x >> pa.v()[i].y >> pa.v()[i].z;

			std::string s;
			icsin >> s;
			if (!isdigit(s[0]))
			{
				icsin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				pa.deathtime_map()[i] = 0;
				pa.id()[i] = i;
				pa.deathflags()[i] = 0;
			}
			else
			{
				pa.deathtime_map()[i] = std::stof(s);
				icsin >> pa.deathflags()[i] >> pa.id()[i];
			}
		}

		return false;
	} 
	*/

	bool load_data_hybrid(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config, std::istream& in)
	{
		std::string s;

		std::getline(in, s);
		std::istringstream ss(s);
		size_t npl;
		ss >> npl;
		
		pl = HostPlanetPhaseSpace(npl, config.tbsize);

		// load the Sun
		std::getline(in, s);
		ss = std::istringstream(s);
		ss >> pl.m()[0];

		pl.id()[0] = 0;
		pl.rplsq()[0] = 0;
		pl.r()[0].x = 0;
		pl.r()[0].y = 0;
		pl.r()[0].z = 0;
		pl.v()[0].x = 0;
		pl.v()[0].y = 0;
		pl.v()[0].z = 0;
		
		// load planets
		for (size_t i = 1; i < npl; i++)
		{
			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> pl.id()[i];

			std::getline(in, s);
			double rpl;
			ss = std::istringstream(s);
			ss >> pl.m()[i] >> rpl;
			pl.rplsq()[i] = rpl * rpl;

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> pl.r()[i].x >> pl.r()[i].y >> pl.r()[i].z;

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> pl.v()[i].x >> pl.v()[i].y >> pl.v()[i].z;
		}


		// load particles
		std::getline(in, s);
		ss = std::istringstream(s);
		size_t npart;
		ss >> npart;
		npart = std::min(npart, static_cast<size_t>(config.max_particle));
		pa = HostParticlePhaseSpace(npart);

		for (size_t i = 0; i < npart; i++)
		{
			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> pa.r()[i].x >> pa.r()[i].y >> pa.r()[i].z;

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> pa.v()[i].x >> pa.v()[i].y >> pa.v()[i].z;

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> pa.id()[i] >> pa.deathflags()[i];

			ss >> pa.deathtime_map()[pa.id()[i]];
		}

		return false;
	}

	bool load_data_hybrid_binary(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config, std::istream& in)
	{
		uint64_t templl;
		read_binary<uint64_t>(in, templl);
		size_t npl = static_cast<size_t>(templl);
		pl = HostPlanetPhaseSpace(npl, config.tbsize);

		read_binary<double>(in, pl.m()[0]);

		pl.id()[0] = 0;
		pl.rplsq()[0] = 0;
		pl.r()[0].x = 0;
		pl.r()[0].y = 0;
		pl.r()[0].z = 0;
		pl.v()[0].x = 0;
		pl.v()[0].y = 0;
		pl.v()[0].z = 0;

		for (size_t i = 1; i < pl.n(); i++)
		{
			double rpl;
			read_binary<uint32_t>(in, pl.id()[i]);
			read_binary<double>(in, pl.m()[i]);
			read_binary<double>(in, rpl);
			read_binary<double>(in, pl.r()[i].x);
			read_binary<double>(in, pl.r()[i].y);
			read_binary<double>(in, pl.r()[i].z);
			read_binary<double>(in, pl.v()[i].x);
			read_binary<double>(in, pl.v()[i].y);
			read_binary<double>(in, pl.v()[i].z);
			pl.rplsq()[i] = rpl * rpl;
		}

		read_binary<uint64_t>(in, templl);
		size_t npart = std::min(static_cast<size_t>(templl), static_cast<size_t>(config.max_particle));
		pa = HostParticlePhaseSpace(npart);

		for (size_t i = 0; i < pa.n(); i++)
		{
			read_binary<uint32_t>(in, pa.id()[i]);
			read_binary<double>(in, pa.r()[i].x);
			read_binary<double>(in, pa.r()[i].y);
			read_binary<double>(in, pa.r()[i].z);
			read_binary<double>(in, pa.v()[i].x);
			read_binary<double>(in, pa.v()[i].y);
			read_binary<double>(in, pa.v()[i].z);
			read_binary<uint16_t>(in, pa.deathflags()[i]);
			float deathtime;

			read_binary<float>(in, deathtime);
			pa.deathtime_map()[pa.id()[i]] = deathtime;
		}

		return !in;
	}

	bool load_data(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config)
	{
		bool ret; // false if successful


		if (!sr::util::does_file_exist(config.hybridin))
		{
			std::ostringstream ss;
			ss << "Input file " << config.hybridin << " does not exist";
			throw std::runtime_error(ss.str());
		}

		if (config.readbinary)
		{
			std::ifstream in(config.hybridin, std::ios_base::binary);
			ret = load_data_hybrid_binary(pl, pa, config, in);
		}
		else
		{
			std::ifstream in(config.hybridin);
			ret = load_data_hybrid(pl, pa, config, in);
		}
		

		if (!ret)
		{
			for (size_t i = 0; i < pl.n(); i++)
			{
				pl.m()[i] *= config.big_g;
			}
		}

		return ret;
	}

	void save_data_hybrid_binary(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, const Configuration& config, std::ostream& out)
	{
		(void) config;
		write_binary(out, static_cast<uint64_t>(pl.n_alive));
		write_binary(out, static_cast<double>(pl.m[0]));
		for (size_t i = 1; i < pl.n_alive; i++)
		{

			write_binary(out, pl.id[i]);
			write_binary(out, pl.m[i]);
			write_binary(out, sqrt(pl.rplsq[i]));
			write_binary(out, pl.r[i].x);
			write_binary(out, pl.r[i].y);
			write_binary(out, pl.r[i].z);
			write_binary(out, pl.v[i].x);
			write_binary(out, pl.v[i].y);
			write_binary(out, pl.v[i].z);
		}

		write_binary(out, static_cast<uint64_t>(pa.n()));
		for (size_t i = 0; i < pa.n(); i++)
		{
			write_binary(out, pa.id()[i]);
			write_binary(out, pa.r()[i].x);
			write_binary(out, pa.r()[i].y);
			write_binary(out, pa.r()[i].z);
			write_binary(out, pa.v()[i].x);
			write_binary(out, pa.v()[i].y);
			write_binary(out, pa.v()[i].z);
			write_binary(out, pa.deathflags()[i]);
			write_binary(out, pa.deathtime_map().at(pa.id()[i]));
		}
	}

	void save_data_hybrid(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, const Configuration& config, std::ostream& out)
	{
		(void) config;
		out << pl.n_alive << std::endl;
		out << std::setprecision(17);

		out << pl.m[0] << std::endl;
		for (size_t i = 1; i < pl.n_alive; i++)
		{
			out << pl.id[i] << std::endl;
			out << pl.m[i] << " " << sqrt(pl.rplsq[i]) << std::endl;
			out << pl.r[i].x << " " << pl.r[i].y << " " << pl.r[i].z << std::endl;
			out << pl.v[i].x << " " << pl.v[i].y << " " << pl.v[i].z << std::endl;
		}

		out << pa.n() << std::endl;
		out << std::setprecision(17);
		for (size_t i = 0; i < pa.n(); i++)
		{
			out << pa.r()[i].x << " " << pa.r()[i].y << " " << pa.r()[i].z << std::endl;
			out << pa.v()[i].x << " " << pa.v()[i].y << " " << pa.v()[i].z << std::endl;
			out << pa.id()[i] << " " << pa.deathflags()[i] << " " << pa.deathtime_map().at(pa.id()[i]) << std::endl;
		}
	}

	/*
	void save_data_nohybrid(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, const Configuration& config, std::ostream& plout, std::ostream& icsout)
	{
		(void) config;
		plout << pl.n_alive << std::endl;
		plout << std::setprecision(17);
		for (size_t i = 0; i < pl.n_alive; i++)
		{
			double m = pl.m[i];
			plout << m << std::endl;

			plout << pl.r[i].x << " " << pl.r[i].y << " " << pl.r[i].z << std::endl;
			plout << pl.v[i].x << " " << pl.v[i].y << " " << pl.v[i].z << std::endl;
		}

		icsout << pa.n() << std::endl;
		icsout << std::setprecision(17);
		for (size_t i = 0; i < pa.n(); i++)
		{
			icsout << pa.r()[i].x << " " << pa.r()[i].y << " " << pa.r()[i].z << std::endl;
			icsout << pa.v()[i].x << " " << pa.v()[i].y << " " << pa.v()[i].z << std::endl;
			icsout << pa.deathtime_map().at(pa.id()[i]) << " " << pa.deathflags()[i] << " " << pa.id()[i] << std::endl;
		}
	}
	*/

	void save_data_swift(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, std::ostream& plout, std::ostream& icsout)
	{
		plout << pl.n_alive << std::endl;
		plout << std::setprecision(17);
		for (size_t i = 0; i < pl.n_alive; i++)
		{
			double m = pl.m[i];
			plout << m << std::endl;

			plout << pl.r[i].x << " " << pl.r[i].y << " " << pl.r[i].z << std::endl;
			plout << pl.v[i].x << " " << pl.v[i].y << " " << pl.v[i].z << std::endl;
		}

		icsout << pa.n() << std::endl;
		icsout << std::setprecision(17);
		for (size_t i = 0; i < pa.n(); i++)
		{
			icsout << pa.r()[i].x << " " << pa.r()[i].y << " " << pa.r()[i].z << std::endl;
			icsout << pa.v()[i].x << " " << pa.v()[i].y << " " << pa.v()[i].z << std::endl;
			icsout << "0 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
			icsout << "0.0 0.0 0.0 0.0 0.0" << std::endl << "0.0 0.0 0.0 0.0 0.0" << std::endl << "0.0 0.0 0.0" << std::endl;
		}
	}

	void save_data(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, const Configuration& config, const std::string& outfile)
	{
		if (config.writebinary)
		{
			std::ofstream out(outfile, std::ios_base::binary);
			save_data_hybrid_binary(pl, pa, config, out);
		}
		else
		{
			std::ofstream out(outfile);
			save_data_hybrid(pl, pa, config, out);
		}
	}

	void write_summary(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, const Configuration& config, const float64_t current_time, 
	const float64_t elapsed_time,  const std::string& outfile)
	{
		Vs pl_death_count = Vs(pl.n, 0);
		Vs pl_encounter_count = Vs(pl.n, 0);
		size_t innerbound_count, outerbound_count, unbound_count, kepler_count, absorb_count, encounter_count;
		innerbound_count = outerbound_count = unbound_count = kepler_count = absorb_count = encounter_count = 0;
		for (size_t i = 0; i < pa.n(); i++)
		{
			auto deathflag = pa.deathflags()[i];
			if (deathflag == 0x01) innerbound_count++;
			else if (deathflag & 0x01)
			{
				auto idx = deathflag >> 8;
				pl_encounter_count[idx]++;
				encounter_count++;
			}
			else if (deathflag & 0x02) outerbound_count++;
			else if (deathflag & 0x08) unbound_count++;
			else if (deathflag & 0x04) kepler_count++;
			else if (deathflag & 0x80)
			{
				auto idx = deathflag >> 8;
				pl_death_count[idx]++;
				absorb_count++;
			}
		}
		std::ofstream out(outfile);
		size_t discard_count = pa.n() - pa.n_alive();

		out << "GLISSER integrated " << pa.n() << " particles with " << pl.n - 1 << " planets from " << config.t_0 << " to " << current_time << " days (or " << std::setprecision(4) << std::setw(5) << current_time/365.25 << " yr/" << current_time/365.25/1e3 << " kyr/" << current_time/365.25/1e6 << " Myr/" << current_time/365.25/1e9 << " Gyr)" << std::endl ;
		out << "The whole integration took " << std::setprecision(4) << std::setw(5) << elapsed_time << " min (or " << elapsed_time/60 << " hr/" 
		<< elapsed_time/60/24 << " days)" << std::endl << std::endl;


		out << std::setprecision(17) << pl.n_alive-1 << "/" << pl.n-1 << " (" << std::setprecision(4) << std::setw(5) << 
		       static_cast<double>(pl.n_alive-1) * 100.0 / static_cast<double>(pl.n-1) << "%)" << " planets are still alive" << std::endl;
		out << std::setprecision(17) << pa.n_alive() << "/" << pa.n() << " (" << std::setprecision(4) << std::setw(5) << 
		       static_cast<double>(pa.n_alive()) * 100.0 / static_cast<double>(pa.n()) <<"%)" << " particles are alive" << std::endl
			<< std::setprecision(17) << encounter_count << "/" << pa.n() << " (" << std::setprecision(4) << std::setw(5) << 
			   static_cast<double>(encounter_count) * 100.0 / static_cast<double>(pa.n()) <<"%)" << " particles are still in encounter" << std::endl;

		if (encounter_count)
		{
			out << std::setprecision(17) << "  Among " << encounter_count << " particles still in encounter:" << std::endl;
			for (size_t i = 1; i < pl.n; i++)
			{
				out << std::setprecision(17) << "    " << pl_encounter_count[i] << " (" << std::setprecision(4) << std::setw(5) << 
					static_cast<double>(pl_encounter_count[i])*100.0/static_cast<double>(encounter_count) << "%) in encounter with planet " << i << std::endl;
			}
		}
		
		out << std::endl;
		out << std::setprecision(17) << discard_count << "/" << pa.n() << " (" << std::setprecision(4) << std::setw(5) << 
		       static_cast<double>(discard_count)*100.0/static_cast<double>(pa.n()) <<"%)" << " particles are discarded" << std::endl; 

		if (discard_count)
		{
			out << std::setprecision(17) << "  Among " << discard_count << " discarded particles:" << std::endl;
			out << std::setprecision(17) << "    " << innerbound_count << " (" << std::setprecision(4) << std::setw(5) << 
			static_cast<double>(innerbound_count) * 100.0 / static_cast<double>(discard_count) << "%) entered inner bound" << std::endl;
			out << std::setprecision(17) << "    " << outerbound_count << " (" << std::setprecision(4) << std::setw(5) << 
			static_cast<double>(outerbound_count) * 100.0 / static_cast<double>(discard_count) << "%) exceeded outer bound" << std::endl;
			out << std::setprecision(17) << "    " << unbound_count << " (" << std::setprecision(4) << std::setw(5) << 
			static_cast<double>(unbound_count) * 100.0 / static_cast<double>(discard_count) << "%) orbits unbound" << std::endl;
			out << std::setprecision(17) << "    " << kepler_count << " (" << std::setprecision(4) << std::setw(5) << 
			static_cast<double>(kepler_count) * 100.0 / static_cast<double>(discard_count) << "%) kepler didn't converge" << std::endl;
			out << std::setprecision(17) << "    " << absorb_count << " (" << std::setprecision(4) << std::setw(5) << 
			static_cast<double>(absorb_count) * 100.0 / static_cast<double>(discard_count) << "%) absorbed by planets" << std::endl;

			if (absorb_count)
			{
				out << std::setprecision(17) << "    Among " << absorb_count << " particles absorbed by planets:" << std::endl;
				for (size_t i = 1; i < pl.n; i++)
				{
					out << std::setprecision(17) << "      " << pl_death_count[i] << " (" << std::setprecision(4) << std::setw(5) << 
						static_cast<double>(pl_death_count[i]) * 100.0 / static_cast<double>(absorb_count) << "%) absorbed by planet " << i << std::endl;
				}
			}
		}
	}

	void save_binary_track(std::ostream& trackout, const HostPlanetSnapshot& pl, const HostParticleSnapshot& pa, double time, bool to_elements, bool barycentric_track, bool single_precision)
	{

		// write time and solar mass
		sr::data::write_binary(trackout, static_cast<double>(time));
		sr::data::write_binary(trackout, static_cast<double>(pl.m[0]));


		// write planets
		if (pl.n_alive > 0)
		{
			sr::data::write_binary(trackout, static_cast<uint64_t>(pl.n_alive - 1));
		}
		else
		{
			sr::data::write_binary(trackout, static_cast<uint64_t>(0));
		}

		double center_mass = pl.m[0];
		f64_3 center_r = pl.r[0] * pl.m[0];
		f64_3 center_v = pl.v[0] * pl.m[0];
		if (barycentric_track)
		{
			sr::convert::find_barycenter(pl.r, pl.v, pl.m, pl.n_alive, center_r, center_v, center_mass);
		}

		if (to_elements) 
		{
			for (uint32_t i = 1; i < pl.n_alive; i++)
			{
			
				double a, e, in, capom, om, f;	
				if (barycentric_track) 
				{
					// sr::convert::get_bary_mu(center_mass, pl.m[i])
					sr::convert::to_elements(center_mass + pl.m[i], pl.r[i] - center_r, pl.v[i] - center_v, nullptr, &a, &e, &in, &capom, &om, &f);
				}
				else 
				{
					sr::convert::to_elements(pl.m[i] + pl.m[0], pl.r[i], pl.v[i],
					nullptr, &a, &e, &in, &capom, &om, &f);
				}

				sr::data::write_binary(trackout, static_cast<uint32_t>(pl.id[i]));

				if (single_precision)
				{
					sr::data::write_binary(trackout, static_cast<float>(pl.m[i]));
					sr::data::write_binary(trackout, static_cast<float>(a));
					sr::data::write_binary(trackout, static_cast<float>(e));
					sr::data::write_binary(trackout, static_cast<float>(in));
					sr::data::write_binary(trackout, static_cast<float>(capom));
					sr::data::write_binary(trackout, static_cast<float>(om));
					sr::data::write_binary(trackout, static_cast<float>(sr::convert::get_mean_anomaly(e, f)));
				}
				else
				{
					sr::data::write_binary(trackout, static_cast<double>(pl.m[i]));
					sr::data::write_binary(trackout, static_cast<double>(a));
					sr::data::write_binary(trackout, static_cast<double>(e));
					sr::data::write_binary(trackout, static_cast<double>(in));
					sr::data::write_binary(trackout, static_cast<double>(capom));
					sr::data::write_binary(trackout, static_cast<double>(om));
					sr::data::write_binary(trackout, static_cast<double>(sr::convert::get_mean_anomaly(e, f)));
				}
			}
		}
		else
		{
			for (uint32_t i = 1; i < pl.n_alive; i++)
			{
				sr::data::write_binary(trackout, static_cast<uint32_t>(pl.id[i]));

				if (single_precision)
				{
					sr::data::write_binary(trackout, static_cast<float>(pl.m[i]));
					sr::data::write_binary(trackout, static_cast<float>(pl.r[i].x - center_r.x));
					sr::data::write_binary(trackout, static_cast<float>(pl.r[i].y - center_r.y));
					sr::data::write_binary(trackout, static_cast<float>(pl.r[i].z - center_r.z));
					sr::data::write_binary(trackout, static_cast<float>(pl.v[i].x - center_v.x));
					sr::data::write_binary(trackout, static_cast<float>(pl.v[i].y - center_v.y));
					sr::data::write_binary(trackout, static_cast<float>(pl.v[i].z - center_v.z));
				}
				else
				{
					sr::data::write_binary(trackout, static_cast<double>(pl.m[i]));
					sr::data::write_binary(trackout, static_cast<double>(pl.r[i].x - center_r.x));
					sr::data::write_binary(trackout, static_cast<double>(pl.r[i].y - center_r.y));
					sr::data::write_binary(trackout, static_cast<double>(pl.r[i].z - center_r.z));
					sr::data::write_binary(trackout, static_cast<double>(pl.v[i].x - center_v.x));
					sr::data::write_binary(trackout, static_cast<double>(pl.v[i].y - center_v.y));
					sr::data::write_binary(trackout, static_cast<double>(pl.v[i].z - center_v.z));
				}
			}
		}

		// write test particles
		sr::data::write_binary(trackout, static_cast<uint64_t>(pa.n_alive));
		for (uint32_t i = 0; i < pa.n_alive; i++)
		{
			if (to_elements)
			{

				double a, e, in, capom, om, f;
				sr::convert::to_elements(center_mass, pa.r[i] - center_r, pa.v[i] - center_v,
						nullptr, &a, &e, &in, &capom, &om, &f);

				sr::data::write_binary(trackout, static_cast<uint32_t>(pa.id[i]));

				if (single_precision)
				{
					sr::data::write_binary(trackout, static_cast<float>(a));
					sr::data::write_binary(trackout, static_cast<float>(e));
					sr::data::write_binary(trackout, static_cast<float>(in));
					sr::data::write_binary(trackout, static_cast<float>(capom));
					sr::data::write_binary(trackout, static_cast<float>(om));
					sr::data::write_binary(trackout, static_cast<float>(sr::convert::get_mean_anomaly(e, f)));
				}
				else
				{
					sr::data::write_binary(trackout, static_cast<double>(a));
					sr::data::write_binary(trackout, static_cast<double>(e));
					sr::data::write_binary(trackout, static_cast<double>(in));
					sr::data::write_binary(trackout, static_cast<double>(capom));
					sr::data::write_binary(trackout, static_cast<double>(om));
					sr::data::write_binary(trackout, static_cast<double>(sr::convert::get_mean_anomaly(e, f)));
				}
			}
			else
			{
				sr::data::write_binary(trackout, static_cast<uint32_t>(pa.id[i]));
				
				if (single_precision)
				{
					sr::data::write_binary(trackout, static_cast<float>(pa.r[i].x - center_r.x));
					sr::data::write_binary(trackout, static_cast<float>(pa.r[i].y - center_r.y));
					sr::data::write_binary(trackout, static_cast<float>(pa.r[i].z - center_r.z));
					sr::data::write_binary(trackout, static_cast<float>(pa.v[i].x - center_v.x));
					sr::data::write_binary(trackout, static_cast<float>(pa.v[i].y - center_v.y));
					sr::data::write_binary(trackout, static_cast<float>(pa.v[i].z - center_v.z));
				}
				else
				{
					sr::data::write_binary(trackout, static_cast<double>(pa.r[i].x - center_r.x));
					sr::data::write_binary(trackout, static_cast<double>(pa.r[i].y - center_r.y));
					sr::data::write_binary(trackout, static_cast<double>(pa.r[i].z - center_r.z));
					sr::data::write_binary(trackout, static_cast<double>(pa.v[i].x - center_v.x));
					sr::data::write_binary(trackout, static_cast<double>(pa.v[i].y - center_v.y));
					sr::data::write_binary(trackout, static_cast<double>(pa.v[i].z - center_v.z));
				}
			}
		}
		trackout.flush();
	}

	void save_txt_hist(std::ostream& histout, const HostPlanetSnapshot& pl, double time, bool single_precision)
	{
		int precision = (single_precision) ? 8 : 17;
		histout << std::setprecision(precision);
		histout << static_cast<double>(time) << " " << static_cast<double>(pl.m[0]) << " ";

		if (pl.n_alive > 0)
		{
			histout << static_cast<uint32_t>(pl.n_alive - 1) << std::endl;
		}
		else
		{
			histout << static_cast<uint32_t>(0) << std::endl;
		}

		// record finished
		for (uint32_t i = 1; i < pl.n_alive; i++)
		{
			double a, e, in, capom, om, f;
			sr::convert::to_elements(pl.m[i] + pl.m[0], pl.r[i], pl.v[i],
				nullptr, &a, &e, &in, &capom, &om, &f);

			histout << static_cast<uint32_t>(pl.id[i]) << " ";
			if (single_precision)
			{
				histout << static_cast<float>(pl.m[i]) << " " << static_cast<float>(pl.rplsq[i]) << " " << static_cast<float>(a) << " " 
				        << static_cast<float>(e) << " " << static_cast<float>(in) << " " << static_cast<float>(capom) << " "
						<< static_cast<float>(om) << " " << static_cast<float>(sr::convert::get_mean_anomaly(e, f)) << std::endl;

			}
			else
			{
				histout << static_cast<double>(pl.m[i]) << " " << static_cast<double>(pl.rplsq[i]) << " " << static_cast<double>(a) << " " 
				        << static_cast<double>(e) << " " << static_cast<double>(in) << " " << static_cast<double>(capom) << " "
						<< static_cast<double>(om) << " " << static_cast<double>(sr::convert::get_mean_anomaly(e, f)) << std::endl;
			}
		}
	}

	void save_binary_hist(std::ostream& histout, const HostPlanetSnapshot& pl, double time, bool single_precision)
	{

		sr::data::write_binary(histout, static_cast<double>(time));
		sr::data::write_binary(histout, static_cast<double>(pl.m[0]));

		if (pl.n_alive > 0)
		{
			sr::data::write_binary(histout, static_cast<uint32_t>(pl.n_alive - 1));
		}
		else
		{
			sr::data::write_binary(histout, static_cast<uint32_t>(0));
		}

		// record finished

		for (uint32_t i = 1; i < pl.n_alive; i++)
		{
			double a, e, in, capom, om, f;
			sr::convert::to_elements(pl.m[i] + pl.m[0], pl.r[i], pl.v[i],
				nullptr, &a, &e, &in, &capom, &om, &f);

			sr::data::write_binary(histout, static_cast<uint32_t>(pl.id[i]));
			if (single_precision)
			{
				sr::data::write_binary(histout, static_cast<float>(pl.m[i]));
				sr::data::write_binary(histout, static_cast<float>(pl.rplsq[i]));
				sr::data::write_binary(histout, static_cast<float>(a));
				sr::data::write_binary(histout, static_cast<float>(e));
				sr::data::write_binary(histout, static_cast<float>(in));
				sr::data::write_binary(histout, static_cast<float>(capom));
				sr::data::write_binary(histout, static_cast<float>(om));
				sr::data::write_binary(histout, static_cast<float>(sr::convert::get_mean_anomaly(e, f)));
			}
			else
			{
				sr::data::write_binary(histout, static_cast<double>(pl.m[i]));
				sr::data::write_binary(histout, static_cast<double>(pl.rplsq[i]));
				sr::data::write_binary(histout, static_cast<double>(a));
				sr::data::write_binary(histout, static_cast<double>(e));
				sr::data::write_binary(histout, static_cast<double>(in));
				sr::data::write_binary(histout, static_cast<double>(capom));
				sr::data::write_binary(histout, static_cast<double>(om));
				sr::data::write_binary(histout, static_cast<double>(sr::convert::get_mean_anomaly(e, f)));
			}
		}

		histout.flush();
	}

	void TrackReader::check_state(const State& expected)
	{
		if (state != expected)
		{
			throw std::runtime_error("invalid reader state");
		}
	}

	void TrackReader::read_time()
	{
		check_state(State::Start);

		sr::data::read_binary<double>(input, time);
		state = State::PlanetsBegin;
	}

	void TrackReader::begin_planets()
	{
		check_state(State::PlanetsBegin);

		planets = HostPlanetSnapshot();

		uint64_t templl;
		sr::data::read_binary<uint64_t>(input, templl);

		n_planets = static_cast<size_t>(templl);
		state = State::PlanetsEnd;
	}

	void TrackReader::read_planets(const Vu32* filter)
	{
		check_state(State::PlanetsEnd);
		std::streampos pos = input.tellg();

		// planets = HostPlanetSnapshot(static_cast<size_t>(n_planets));
		planets = HostPlanetSnapshot(static_cast<size_t>(filter ? filter->size() : n_planets));

		int index = 0;
		for (uint32_t i = 0; i < n_planets; i++)
		{
			uint32_t id;
			sr::data::read_binary<uint32_t>(input, id);

			if (!filter || std::find(filter->begin(), filter->end(), id) != filter->end())
			{
				planets.id[index] = id;

				float tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				planets.r[index].x = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				planets.r[index].y = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				planets.r[index].z = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				planets.v[index].x = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				planets.v[index].y = tempfloat;
				sr::data::read_binary<float>(input, tempfloat);
				planets.v[index].z = tempfloat;
				index++;
			}
			else
			{
				input.seekg(4 * 6, std::ios_base::cur);
			}
		}

		input.seekg(pos);
	}

	void TrackReader::end_planets()
	{
		check_state(State::PlanetsEnd);

		input.seekg(TRACK_PLANET_STRIDE * n_planets, std::ios_base::cur);
		state = State::ParticlesBegin;
	}

	void TrackReader::begin_particles()
	{
		check_state(State::ParticlesBegin);

		particles = HostParticleSnapshot();

		uint64_t templl;
		sr::data::read_binary<uint64_t>(input, templl);

		n_particles = static_cast<size_t>(templl);
		state = State::ParticlesEnd;
	}

	void TrackReader::read_particles()
	{
		check_state(State::ParticlesEnd);

		particles = HostParticleSnapshot(n_particles);
		std::streampos position = input.tellg();

		for (uint32_t i = 0; i < n_particles; i++)
		{
			sr::data::read_binary<uint32_t>(input, particles.id[i]);

			float tempfloat;
			sr::data::read_binary<float>(input, tempfloat);
			particles.r[i].x = tempfloat;
			sr::data::read_binary<float>(input, tempfloat);
			particles.r[i].y = tempfloat;
			sr::data::read_binary<float>(input, tempfloat);
			particles.r[i].z = tempfloat;
			sr::data::read_binary<float>(input, tempfloat);
			particles.v[i].x = tempfloat;
			sr::data::read_binary<float>(input, tempfloat);
			particles.v[i].y = tempfloat;
			sr::data::read_binary<float>(input, tempfloat);
			particles.v[i].z = tempfloat;
		}

		input.seekg(position);
	}

	bool TrackReader::read_particle(uint32_t id)
	{
		check_state(State::ParticlesEnd);

		std::streampos position = input.tellg();

		size_t result = bsearch_track(input, n_particles, id, TRACK_PARTICLE_STRIDE);

		if (result)
		{
			size_t newsize = particles.n + 1;
			particles.n = newsize;
			particles.n_alive = newsize;

			particles.r.resize(newsize);
			particles.v.resize(newsize);
			particles.id.resize(newsize);

			particles.r[newsize - 1].x = sr::data::read_binary<float>(input);
			particles.r[newsize - 1].y = sr::data::read_binary<float>(input);
			particles.r[newsize - 1].z = sr::data::read_binary<float>(input);
			particles.v[newsize - 1].x = sr::data::read_binary<float>(input);
			particles.v[newsize - 1].y = sr::data::read_binary<float>(input);
			particles.v[newsize - 1].z = sr::data::read_binary<float>(input);

			particles.id[newsize - 1] = id;
		}

		input.seekg(position);

		return result > 0;
	}

	void TrackReader::end_particles()
	{
		check_state(State::ParticlesEnd);
		input.seekg(TRACK_PARTICLE_STRIDE * n_particles, std::ios_base::cur);

		state = State::Finish;
	}

	size_t TrackReader::bsearch_track(std::istream& f, size_t npa, uint32_t partnum, size_t stride)
	{
		int64_t base = f.tellg();
		int64_t left = 0;
		int64_t right = npa - 1;

		size_t its = 0;
		while (left <= right)
		{
			its++;
			int64_t mid = left + (right - left) / 2;
			f.seekg(base + mid * stride);

			uint32_t id;
			sr::data::read_binary<uint32_t>(f, id);
			if (id == partnum)
			{
				return its;
			}
			else if (id > partnum)
			{
				right = mid - 1;
			}
			else
			{
				left = mid + 1;
			}
		}

		return 0;
	}

	void process_track(std::istream& input,
		const TrackReaderOptions& options,
		const std::function<void(HostPlanetSnapshot&, HostParticleSnapshot&, double)>& callback)
	{
		while (true)
		{
			TrackReader reader(input);
			bool skip = false;

			reader.read_time();
			if (!input) break;

			if (reader.time > options.max_time)
			{
				skip = true;
			}

			reader.begin_planets();
			if (!input) break;

			if (!skip)
			{
				reader.read_planets(options.take_all_planets ? nullptr : &options.planet_filter);
			}
			reader.end_planets();

			reader.begin_particles();
			if (!input) break;

			if (!skip)
			{
				if (options.take_all_particles)
				{
					reader.read_particles();
				}
				else
				{
					for (uint32_t id : options.particle_filter)
					{
						reader.read_particle(id);
					}
				}
			}
			reader.end_particles();

			if (!input) break;

			if (!skip)
			{
				callback(reader.planets, reader.particles, reader.time);
			}
		}
	}

	void read_tracks(const std::string& path,
		const TrackReaderOptions& options,
		const std::function<void(HostPlanetSnapshot&, HostParticleSnapshot&, double)>& callback)
	{
		sr::util::PathType pathtype = sr::util::get_path_type(path);

		if (pathtype == sr::util::PathType::Directory)
		{
			for (size_t i = 0; true; i++)
			{
				std::ostringstream ss;
				ss << path;
				if (path[path.size() - 1] != '/') ss << '/';
				ss << "track." << i << ".out";

				if (!options.silent)
				{
					std::cout << "Reading " << ss.str() << std::endl;
				}

				if (!sr::util::does_file_exist(ss.str()))
				{
					if (i == 0)
					{
						throw std::runtime_error("Directory not contain track files");
					}
					else
					{
						if (!options.silent)
						{
							std::cout << i << " files read" << std::endl;
						}
						return;
					}
				}

				std::ifstream input(ss.str(), std::ios_base::binary);

				process_track(input, options, callback);
			}
		}
		else if (pathtype == sr::util::PathType::File)
		{
			std::ifstream input(path, std::ios_base::binary);

			process_track(input, options, callback);
		}
		else
		{
			throw std::runtime_error("Path does not exist");
		}
	}

}
}
