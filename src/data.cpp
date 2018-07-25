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

namespace sr
{
namespace data
{
	size_t stable_partition_alive_indices(const std::vector<uint16_t>& flags, size_t begin, size_t length, std::unique_ptr<std::vector<size_t>>* indices)
	{
		auto new_indices = std::make_unique<std::vector<size_t>>(length);
		std::iota(new_indices->begin(), new_indices->end(), 0);

		auto val = std::stable_partition(new_indices->begin(), new_indices->end(), [begin, &flags](size_t index)
				{ return (flags[index + begin] & 0x00FE) == 0; }) - new_indices->begin() + begin;

		*indices = std::move(new_indices);
		return val;
	}

	size_t stable_partition_unflagged_indices(const std::vector<uint16_t>& flags, size_t begin, size_t length, std::unique_ptr<std::vector<size_t>>* indices)
	{
		auto new_indices = std::make_unique<std::vector<size_t>>(length);
		std::iota(new_indices->begin(), new_indices->end(), 0);

		auto val = std::stable_partition(new_indices->begin(), new_indices->end(), [begin, &flags](size_t index)
				{ return (flags[index + begin] & 0x00FF) == 0; }) - new_indices->begin() + begin;

		*indices = std::move(new_indices);
		return val;
	}

	void HostParticleSnapshot::gather(const std::vector<size_t>& indices, size_t begin, size_t length)
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

	void HostParticlePhaseSpace::gather(const std::vector<size_t>& indices, size_t begin, size_t length)
	{
		base.gather(indices, begin, length);
		sr::data::gather(deathtime(), indices, begin, length);
		sr::data::gather(deathflags(), indices, begin, length);

		if (deathtime_index().size() > 0)
		{
			sr::data::gather(deathtime_index(), indices, begin, length);
		}
	}

	std::unique_ptr<std::vector<size_t>> HostParticlePhaseSpace::stable_partition_unflagged(size_t begin, size_t length)
	{
		std::unique_ptr<std::vector<size_t>> indices;
		n_alive() = stable_partition_unflagged_indices(deathflags(), begin, length, &indices);
		
		this->gather(*indices, begin, length);

		return indices;
	}

	std::unique_ptr<std::vector<size_t>> HostParticlePhaseSpace::stable_partition_alive(size_t begin, size_t length)
	{
		std::unique_ptr<std::vector<size_t>> indices;
		n_alive() = stable_partition_alive_indices(deathflags(), begin, length, &indices);
		
		this->gather(*indices, begin, length);

		return indices;
	}

	std::unique_ptr<std::vector<size_t>> HostParticleSnapshot::sort_by_id(size_t begin, size_t length)
	{
		auto new_indices = std::make_unique<std::vector<size_t>>(length);

		std::iota(new_indices->begin(), new_indices->end(), 0);
		std::sort(new_indices->begin(), new_indices->end(), [begin, this](size_t a, size_t b)
				{ return id[a + begin] < id[b + begin]; });

		this->gather(*new_indices, 0, length);

		return new_indices;
	}

	std::unique_ptr<std::vector<size_t>> HostParticlePhaseSpace::sort_by_id(size_t begin, size_t length)
	{
		auto indices = base.sort_by_id(begin, length);

		this->gather(*indices, 0, length);
		return indices;
	}

	Configuration::Configuration()
	{
		num_thread = 4;
		t_0 = 0;
		t_f = 365e4;
		dt = 122;
		tbsize = 1024;
		wh_ce_n1 = 8;
		wh_ce_n2 = 4;
		wh_ce_r1 = 1;
		wh_ce_r2 = 3.5;
		cull_radius = 0.5;

		print_every = 10;
		energy_every = 1;
		track_every = 0;
		split_track_file = 0;

		use_gpu = true;

		dump_every = 1000;
		max_particle = static_cast<uint32_t>(-1);
		resolve_encounters = false;
		big_g = 1;
		icsin = "";
		plin = "";
		hybridin = "";
		hybridout = "state.out";
		readsplit = 0;
		writesplit = 0;
		writebinary = 0;
		readbinary = 0;
		outfolder = "output/";
		readmomenta = false;
		writemomenta = false;
	}

	Configuration Configuration::output_config() const
	{
#pragma GCC warning "TODO"
		return *this;
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
				if (first == "Initial-Time")
					out->t_0 = std::stod(second);
				else if (first == "Time-Step")
					out->dt = std::stod(second);
				else if (first == "Final-Time")
					out->t_f = std::stod(second);
				else if (first == "Time-Block-Size")
					out->tbsize = std::stou(second);
				else if (first == "Big-G")
					out->big_g = std::stod(second);
				else if (first == "Cull-Radius")
					out->cull_radius = std::stod(second);
				else if (first == "WH-Encounter-N1")
					out->wh_ce_n1 = std::stou(second);
				else if (first == "WH-Encounter-N2")
					out->wh_ce_n2 = std::stou(second);
				else if (first == "WH-Encounter-R1")
					out->wh_ce_r1 = std::stod(second);
				else if (first == "WH-Encounter-R2")
					out->wh_ce_r2 = std::stod(second);
				else if (first == "Enable-GPU")
					out->use_gpu = std::stoi(second) != 0;
				else if (first == "CPU-Thread-Count")
					out->num_thread = std::stou(second);
				else if (first == "Limit-Particle-Count")
					out->max_particle = std::stou(second);
				else if (first == "Log-Interval")
					out->print_every = std::stou(second);
				else if (first == "Status-Interval")
					out->energy_every = std::stou(second);
				else if (first == "Track-Interval")
					out->track_every = std::stou(second);
				else if (first == "Split-Track-File")
					out->split_track_file = std::stou(second);
				else if (first == "Dump-Interval")
					out->dump_every = std::stou(second);
				else if (first == "Resolve-Encounters")
					out->resolve_encounters = std::stoi(second) != 0;
				else if (first == "Write-Split-Output")
					out->writesplit = std::stoi(second) != 0;
				else if (first == "Read-Split-Input")
					out->readsplit = std::stoi(second) != 0;
				else if (first == "Write-Binary-Output")
					out->writebinary = std::stoi(second) != 0;
				else if (first == "Read-Binary-Input")
					out->readbinary = std::stoi(second) != 0;
				else if (first == "Input-File")
					out->hybridin = second;
				else if (first == "Output-File")
					out->hybridout = second;
				else if (first == "Particle-Input-File")
					out->icsin = second;
				else if (first == "Planet-Input-File")
					out->plin = second;
				else if (first == "Output-Folder")
					out->outfolder = second;
				else if (first == "Read-Input-Momenta")
					out->readmomenta = std::stoi(second) != 0;
				else if (first == "Write-Output-Momenta")
					out->writemomenta = std::stoi(second) != 0;
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
		if (out->cull_radius != 0.5 && out->resolve_encounters)
		{
			std::cerr << "Warning: Cull-Radius was set but Resolve-Encounters was also set: Planets will not be culled with Cull-Radius!" << std::endl;
		}
		if (out->readsplit && out->readbinary)
		{
			std::cerr << "Warning: Read-Split-Input was selected but Read-Binary-Input was also specified. Ignoring Read-Binary-Input" << std::endl;
		}
		if (out->writesplit && out->writebinary)
		{
			std::cerr << "Warning: Write-Split-Input was selected but Write-Binary-Input was also specified. Ignoring Write-Binary-Input" << std::endl;
		}

		if (!out->writesplit && out->hybridout == "")
		{
			throw std::runtime_error("Error: Output-File was not specified");
		}
		if (!out->writesplit && out->hybridin == "")
		{
			throw std::runtime_error("Error: Input-File was not specified");
		}
		if (out->readsplit && (out->plin == "" || out->icsin == ""))
		{
			throw std::runtime_error("Error: Read-Split-Input was selected but Particle-Input-File or Planet-Input-File were not specified");
		}
	}

	void write_configuration(std::ostream& outstream, const Configuration& out)
	{
		outstream << std::setprecision(17);
		outstream << "Initial-Time " << out.t_0 << std::endl;
		outstream << "Time-Step " << out.dt << std::endl;
		outstream << "Final-Time " << out.t_f << std::endl;
		outstream << "Time-Block-Size " << out.tbsize << std::endl;
		outstream << "Cull-Radius " << out.cull_radius << std::endl;
		outstream << "WH-Encounter-N1 " << out.wh_ce_n1 << std::endl;
		outstream << "WH-Encounter-N2 " << out.wh_ce_n2 << std::endl;
		outstream << "WH-Encounter-R1 " << out.wh_ce_r1 << std::endl;
		outstream << "WH-Encounter-R2 " << out.wh_ce_r2 << std::endl;
		outstream << "Enable-GPU " << out.use_gpu << std::endl;
		outstream << "CPU-Thread-Count" << out.num_thread << std::endl;
		outstream << "Limit-Particle-Count " << out.max_particle << std::endl;
		outstream << "Log-Interval " << out.print_every << std::endl;
		outstream << "Status-Interval " << out.energy_every << std::endl;
		outstream << "Track-Interval " << out.track_every << std::endl;
		outstream << "Split-Track-File " << out.split_track_file << std::endl;
		outstream << "Dump-Interval " << out.dump_every << std::endl;
		outstream << "Resolve-Encounters " << out.resolve_encounters << std::endl;
		outstream << "Write-Split-Output " << out.writesplit << std::endl;
		outstream << "Write-Binary-Output " << out.writebinary << std::endl;
		outstream << "Read-Split-Input " << out.readsplit << std::endl;
		outstream << "Read-Binary-Input " << out.readbinary << std::endl;
		outstream << "Input-File " << out.hybridin << std::endl;
		outstream << "Output-File " << out.hybridout << std::endl;
		outstream << "Particle-Input-File " << out.icsin << std::endl;
		outstream << "Planet-Input-File " << out.plin << std::endl;
		outstream << "Output-Folder " << out.outfolder << std::endl;
		outstream << "Read-Input-Momenta " << out.readmomenta << std::endl;
		outstream << "Write-Output-Momenta " << out.writemomenta << std::endl;
	}

	bool load_data_nohybrid(HostData& hd, const Configuration& config, std::istream& plin, std::istream& icsin)
	{
		size_t npl;
		plin >> npl;

		hd.planets = HostPlanetPhaseSpace(npl, config.tbsize, config.fast_timestep_factor());

		for (size_t i = 0; i < npl; i++)
		{
			plin >> hd.planets.m()[i];
			plin >> hd.planets.r()[i].x >> hd.planets.r()[i].y >> hd.planets.r()[i].z;
			plin >> hd.planets.v()[i].x >> hd.planets.v()[i].y >> hd.planets.v()[i].z;

			hd.planets.id()[i] = static_cast<uint32_t>(i);
		}

		size_t npart;
		icsin >> npart;
		npart = std::min(npart, static_cast<size_t>(config.max_particle));

		hd.particles = HostParticlePhaseSpace(npart, !config.use_gpu);

		for (size_t i = 0; i < npart; i++)
		{
			icsin >> hd.particles.r()[i].x >> hd.particles.r()[i].y >> hd.particles.r()[i].z;
			icsin >> hd.particles.v()[i].x >> hd.particles.v()[i].y >> hd.particles.v()[i].z;

			std::string s;
			icsin >> s;
			if (!isdigit(s[0]))
			{
				icsin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				hd.particles.deathtime()[i] = 0;
				hd.particles.id()[i] = static_cast<uint32_t>(i);
				hd.particles.deathflags()[i] = 0;
			}
			else
			{
				hd.particles.deathtime()[i] = std::stof(s);
				icsin >> hd.particles.deathflags()[i] >> hd.particles.id()[i];
			}
		}

		return false;
	}

	bool load_data_hybrid(HostData& hd, const Configuration& config, std::istream& in)
	{
		std::string s;

		std::getline(in, s);
		std::istringstream ss(s);
		size_t npl;
		ss >> npl;
		
		hd.planets = HostPlanetPhaseSpace(npl, config.tbsize, config.fast_timestep_factor());

		for (size_t i = 0; i < npl; i++)
		{
			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> hd.planets.m()[i];

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> hd.planets.r()[i].x >> hd.planets.r()[i].y >> hd.planets.r()[i].z;

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> hd.planets.v()[i].x >> hd.planets.v()[i].y >> hd.planets.v()[i].z;

			std::getline(in, s);
			ss = std::istringstream(s);
			
			ss >> hd.planets.id()[i];

			hd.planets.id()[i] = static_cast<uint32_t>(i);
		}

		std::getline(in, s);
		ss = std::istringstream(s);
		size_t npart;
		ss >> npart;
		npart = std::min(npart, static_cast<size_t>(config.max_particle));
		hd.particles = HostParticlePhaseSpace(npart, !config.use_gpu);

		for (size_t i = 0; i < npart; i++)
		{
			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> hd.particles.r()[i].x >> hd.particles.r()[i].y >> hd.particles.r()[i].z;

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> hd.particles.v()[i].x >> hd.particles.v()[i].y >> hd.particles.v()[i].z;

			std::getline(in, s);
			ss = std::istringstream(s);
			ss >> hd.particles.id()[i] >> hd.particles.deathflags()[i] >> hd.particles.deathtime()[i];
		}

		return false;
	}

	bool load_data_hybrid_binary(HostData& hd, const Configuration& config, std::istream& in)
	{
		uint64_t templl;
		read_binary<uint64_t>(in, templl);
		size_t npl = static_cast<size_t>(templl);
		hd.planets = HostPlanetPhaseSpace(npl, config.tbsize, config.fast_timestep_factor());

		for (size_t i = 0; i < hd.planets.n(); i++)
		{
			read_binary<uint32_t>(in, hd.planets.id()[i]);
			read_binary<double>(in, hd.planets.m()[i]);
			read_binary<double>(in, hd.planets.r()[i].x);
			read_binary<double>(in, hd.planets.r()[i].y);
			read_binary<double>(in, hd.planets.r()[i].z);
			read_binary<double>(in, hd.planets.v()[i].x);
			read_binary<double>(in, hd.planets.v()[i].y);
			read_binary<double>(in, hd.planets.v()[i].z);
		}

		read_binary<uint64_t>(in, templl);
		size_t npart = std::min(static_cast<size_t>(templl), static_cast<size_t>(config.max_particle));
		hd.particles = HostParticlePhaseSpace(npart, !config.use_gpu);

		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			read_binary<uint32_t>(in, hd.particles.id()[i]);
			read_binary<double>(in, hd.particles.r()[i].x);
			read_binary<double>(in, hd.particles.r()[i].y);
			read_binary<double>(in, hd.particles.r()[i].z);
			read_binary<double>(in, hd.particles.v()[i].x);
			read_binary<double>(in, hd.particles.v()[i].y);
			read_binary<double>(in, hd.particles.v()[i].z);
			read_binary<uint16_t>(in, hd.particles.deathflags()[i]);
			read_binary<float>(in, hd.particles.deathtime()[i]);
		}

		return !in;
	}

	bool load_data(HostData& hd, const Configuration& config)
	{
		bool ret;
		if (config.readsplit)
		{
			if (!sr::util::does_file_exist(config.plin))
			{
				std::ostringstream ss;
				ss << "Planet input file " << config.plin << " does not exist";
				throw std::runtime_error(ss.str());
			}
			if (!sr::util::does_file_exist(config.icsin))
			{
				std::ostringstream ss;
				ss << "Particle input file " << config.icsin << " does not exist";
				throw std::runtime_error(ss.str());
			}
			
			std::ifstream plinfile(config.plin), icsinfile(config.icsin);
			ret = load_data_nohybrid(hd, config, plinfile, icsinfile);
		}
		else
		{
			if (!sr::util::does_file_exist(config.hybridin))
			{
				std::ostringstream ss;
				ss << "Input file " << config.hybridin << " does not exist";
				throw std::runtime_error(ss.str());
			}

			if (config.readbinary)
			{
				std::ifstream in(config.hybridin, std::ios_base::binary);
				ret = load_data_hybrid_binary(hd, config, in);
			}
			else
			{
				std::ifstream in(config.hybridin);
				ret = load_data_hybrid(hd, config, in);
			}
		}

		if (!ret)
		{
			for (size_t i = 0; i < hd.planets.n(); i++)
			{
				if (config.readmomenta)
				{
					hd.planets.v()[i] /= hd.planets.m()[i];
				}

				hd.planets.m()[i] *= config.big_g;
			}

			hd.particles.stable_partition_alive(0, hd.particles.n());
		}

		return ret;
	}

	void save_data_hybrid_binary(const HostData& hd, const Configuration& config, std::ostream& out)
	{
		(void) config;
		write_binary(out, static_cast<uint64_t>(hd.planets.n_alive()));
		for (size_t i = 0; i < hd.planets.n_alive(); i++)
		{
			double m = hd.planets_snapshot.m[i];
			write_binary(out, hd.planets_snapshot.id[i]);
			write_binary(out, m);

			if (!config.writemomenta) m = 1;
			write_binary(out, hd.planets_snapshot.r[i].x);
			write_binary(out, hd.planets_snapshot.r[i].y);
			write_binary(out, hd.planets_snapshot.r[i].z);
			write_binary(out, hd.planets_snapshot.v[i].x * m);
			write_binary(out, hd.planets_snapshot.v[i].y * m);
			write_binary(out, hd.planets_snapshot.v[i].z * m);
		}

		write_binary(out, static_cast<uint64_t>(hd.particles.n()));
		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			write_binary(out, hd.particles.id()[i]);
			write_binary(out, hd.particles.r()[i].x);
			write_binary(out, hd.particles.r()[i].y);
			write_binary(out, hd.particles.r()[i].z);
			write_binary(out, hd.particles.v()[i].x);
			write_binary(out, hd.particles.v()[i].y);
			write_binary(out, hd.particles.v()[i].z);
			write_binary(out, hd.particles.deathflags()[i]);
			write_binary(out, hd.particles.deathtime()[i]);
		}
	}

	void save_data_hybrid(const HostData& hd, const Configuration& config, std::ostream& out)
	{
		(void) config;
		out << hd.planets_snapshot.n_alive << std::endl;
		out << std::setprecision(17);
		for (size_t i = 0; i < hd.planets_snapshot.n_alive; i++)
		{
			double m = hd.planets_snapshot.m[i];
			out << m << std::endl;
			if (!config.writemomenta) m = 1;
			out << hd.planets_snapshot.r[i].x << " " << hd.planets_snapshot.r[i].y << " " << hd.planets_snapshot.r[i].z << std::endl;
			out << hd.planets_snapshot.v[i].x * m << " " << hd.planets_snapshot.v[i].y * m << " " << hd.planets_snapshot.v[i].z * m << std::endl;
			out << hd.planets_snapshot.id[i] << std::endl;
		}

		out << hd.particles.n() << std::endl;
		out << std::setprecision(17);
		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			out << hd.particles.r()[i].x << " " << hd.particles.r()[i].y << " " << hd.particles.r()[i].z << std::endl;
			out << hd.particles.v()[i].x << " " << hd.particles.v()[i].y << " " << hd.particles.v()[i].z << std::endl;
			out << hd.particles.id()[i] << " " << hd.particles.deathflags()[i] << " " << hd.particles.deathtime()[i] << std::endl;
		}
	}

	void save_data_nohybrid(const HostData& hd, const Configuration& config, std::ostream& plout, std::ostream& icsout)
	{
		(void) config;
		plout << hd.planets_snapshot.n_alive << std::endl;
		plout << std::setprecision(17);
		for (size_t i = 0; i < hd.planets_snapshot.n_alive; i++)
		{
			double m = hd.planets_snapshot.m[i];
			plout << m << std::endl;

			if (!config.writemomenta) m = 1;
			plout << hd.planets_snapshot.r[i].x << " " << hd.planets_snapshot.r[i].y << " " << hd.planets_snapshot.r[i].z << std::endl;
			plout << hd.planets_snapshot.v[i].x * m << " " << hd.planets_snapshot.v[i].y * m << " " << hd.planets_snapshot.v[i].z * m << std::endl;
		}

		icsout << hd.particles.n() << std::endl;
		icsout << std::setprecision(17);
		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			icsout << hd.particles.r()[i].x << " " << hd.particles.r()[i].y << " " << hd.particles.r()[i].z << std::endl;
			icsout << hd.particles.v()[i].x << " " << hd.particles.v()[i].y << " " << hd.particles.v()[i].z << std::endl;
			icsout << hd.particles.deathtime()[i] << " " << hd.particles.deathflags()[i] << " " << hd.particles.id()[i] << std::endl;
		}
	}

	void save_data(const HostData& hd, const Configuration& config, const std::string& outfile, bool dump)
	{
		if (dump || !config.writesplit)
		{
			std::ostringstream ss;

			if (dump || (!config.writesplit && config.writebinary))
			{
				std::ofstream out(outfile, std::ios_base::binary);
				save_data_hybrid_binary(hd, config, out);
			}
			else
			{
				std::ofstream out(outfile);
				save_data_hybrid(hd, config, out);
			}
		}
		else
		{
			std::ofstream ploutfile(sr::util::joinpath(config.outfolder, "pl.out")), icsoutfile(sr::util::joinpath(config.outfolder, "ics.out"));
			save_data_nohybrid(hd, config, ploutfile, icsoutfile);
		}
	}

	void save_binary_track(std::ostream& trackout, const HostPlanetSnapshot& pl, const HostParticleSnapshot& pa, double time, bool to_elements)
	{
		sr::data::write_binary(trackout, static_cast<double>(time));
		sr::data::write_binary(trackout, static_cast<uint64_t>(pl.n_alive - 1));

		for (uint32_t i = 1; i < pl.n_alive; i++)
		{
			if (to_elements)
			{
				double a, e, in, capom, om, f;
				sr::convert::to_elements(pl.m[i] + pl.m[0], pl.r[i], pl.v[i],
					nullptr, &a, &e, &in, &capom, &om, &f);

				sr::data::write_binary(trackout, static_cast<uint32_t>(pl.id[i]));
				sr::data::write_binary(trackout, static_cast<float>(a));
				sr::data::write_binary(trackout, static_cast<float>(e));
				sr::data::write_binary(trackout, static_cast<float>(in));
				sr::data::write_binary(trackout, static_cast<float>(capom));
				sr::data::write_binary(trackout, static_cast<float>(om));
				sr::data::write_binary(trackout, static_cast<float>(f));
			}
			else
			{
				sr::data::write_binary(trackout, static_cast<uint32_t>(pl.id[i]));
				sr::data::write_binary(trackout, static_cast<float>(pl.r[i].x));
				sr::data::write_binary(trackout, static_cast<float>(pl.r[i].y));
				sr::data::write_binary(trackout, static_cast<float>(pl.r[i].z));
				sr::data::write_binary(trackout, static_cast<float>(pl.v[i].x));
				sr::data::write_binary(trackout, static_cast<float>(pl.v[i].y));
				sr::data::write_binary(trackout, static_cast<float>(pl.v[i].z));
			}
		}

		sr::data::write_binary(trackout, static_cast<uint64_t>(pa.n_alive));
		for (uint32_t i = 0; i < pa.n_alive; i++)
		{
			if (to_elements)
			{
				double a, e, in, capom, om, f;
				sr::convert::to_elements(pl.m[0], pa.r[i], pa.v[i],
						nullptr, &a, &e, &in, &capom, &om, &f);

				sr::data::write_binary(trackout, static_cast<uint32_t>(pa.id[i]));
				sr::data::write_binary(trackout, static_cast<float>(a));
				sr::data::write_binary(trackout, static_cast<float>(e));
				sr::data::write_binary(trackout, static_cast<float>(in));
				sr::data::write_binary(trackout, static_cast<float>(capom));
				sr::data::write_binary(trackout, static_cast<float>(om));
				sr::data::write_binary(trackout, static_cast<float>(f));
			}
			else
			{
				sr::data::write_binary(trackout, static_cast<uint32_t>(pa.id[i]));
				sr::data::write_binary(trackout, static_cast<float>(pa.r[i].x));
				sr::data::write_binary(trackout, static_cast<float>(pa.r[i].y));
				sr::data::write_binary(trackout, static_cast<float>(pa.r[i].z));
				sr::data::write_binary(trackout, static_cast<float>(pa.v[i].x));
				sr::data::write_binary(trackout, static_cast<float>(pa.v[i].y));
				sr::data::write_binary(trackout, static_cast<float>(pa.v[i].z));
			}
		}

		trackout.flush();
	}

	void load_binary_track(std::istream& trackin, HostPlanetSnapshot& pl, HostParticleSnapshot& pa, double& time, bool skipplanets, bool skipparticles)
	{
		sr::data::read_binary<double>(trackin, time);

		uint64_t templl;
		sr::data::read_binary<uint64_t>(trackin, templl);

		if (!trackin) return;

		if (skipplanets)
		{
			pl.n = static_cast<size_t>(templl);
			pl.n_alive = static_cast<size_t>(templl);

			trackin.seekg(TRACK_PLANET_STRIDE * pl.n, std::ios_base::cur);
		}
		else
		{
			pl = HostPlanetSnapshot(static_cast<size_t>(templl + 1));
			
			for (uint32_t i = 1; i < pl.n_alive; i++)
			{
				sr::data::read_binary<uint32_t>(trackin, pl.id[i]);

				float tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pl.r[i].x = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pl.r[i].y = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pl.r[i].z = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pl.v[i].x = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pl.v[i].y = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pl.v[i].z = tempfloat;
			}
		}

		sr::data::read_binary(trackin, templl);
		if (!trackin) return;

		if (skipparticles)
		{
			pa.n = static_cast<size_t>(templl);
			pa.n_alive = static_cast<size_t>(templl);

			trackin.seekg(TRACK_PARTICLE_STRIDE * pa.n, std::ios_base::cur);
		}
		else
		{
			pa = HostParticleSnapshot(static_cast<size_t>(templl));

			for (uint32_t i = 0; i < pa.n_alive; i++)
			{
				sr::data::read_binary<uint32_t>(trackin, pa.id[i]);

				float tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pa.r[i].x = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pa.r[i].y = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pa.r[i].z = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pa.v[i].x = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pa.v[i].y = tempfloat;
				sr::data::read_binary<float>(trackin, tempfloat);
				pa.v[i].z = tempfloat;
			}
		}
	}

	size_t bsearch_track(std::istream& f, size_t npa, uint32_t partnum, size_t stride)
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

	void process_track(std::istream& input, bool takeallparticles, const std::vector<uint32_t>& particles, bool killplanets, const std::function<void(HostPlanetSnapshot&, HostParticleSnapshot&, double)>& callback)
	{
		while (true)
		{
			double time;
			HostParticleSnapshot paout;
			HostPlanetSnapshot pl;

			load_binary_track(input, pl, paout, time, killplanets, !takeallparticles);

			if (!input) break;

			if (!takeallparticles)
			{
				input.seekg(-TRACK_PARTICLE_STRIDE * paout.n, std::ios_base::cur);

				for (size_t j = 0; j < particles.size(); j++)
				{
					int64_t position = input.tellg();
					size_t result = bsearch_track(input, paout.n, particles[j], TRACK_PARTICLE_STRIDE);

					if (result)
					{
						paout.r.resize(paout.r.size() + 1);
						paout.v.resize(paout.v.size() + 1);
						paout.id.resize(paout.id.size() + 1);

						paout.r[paout.r.size() - 1].x = sr::data::read_binary<float>(input);
						paout.r[paout.r.size() - 1].y = sr::data::read_binary<float>(input);
						paout.r[paout.r.size() - 1].z = sr::data::read_binary<float>(input);
						paout.v[paout.r.size() - 1].x = sr::data::read_binary<float>(input);
						paout.v[paout.r.size() - 1].y = sr::data::read_binary<float>(input);
						paout.v[paout.r.size() - 1].z = sr::data::read_binary<float>(input);

						paout.id[paout.id.size() - 1] = particles[j];
					}

					input.seekg(position);
				}

				input.seekg(TRACK_PARTICLE_STRIDE * paout.n, std::ios_base::cur);
				paout.n_alive = paout.n = paout.r.size();
			}

			callback(pl, paout, time);
		}
	}

	void read_tracks(const std::string& path, bool takeallparticles, const std::vector<uint32_t>& particle_filter, bool removeplanets,
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

				std::cout << "Reading " << ss.str() << std::endl;

				if (!sr::util::does_file_exist(ss.str()))
				{
					if (i == 0)
					{
						throw std::runtime_error("Directory not contain track files");
					}
					else
					{
						std::cout << i << " files read" << std::endl;
						return;
					}
				}

				std::ifstream input(ss.str(), std::ios_base::binary);

				process_track(input, takeallparticles, particle_filter, removeplanets, callback);
			}
		}
		else if (pathtype == sr::util::PathType::File)
		{
			std::ifstream input(path, std::ios_base::binary);

			process_track(input, takeallparticles, particle_filter, removeplanets, callback);
		}
		else
		{
			throw std::runtime_error("Path does not exist");
		}
	}
}
}
