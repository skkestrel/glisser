#include "data.h"
#include "util.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>

size_t stable_partition_alive_indices(const std::vector<uint16_t>& flags, size_t begin, size_t length, std::unique_ptr<std::vector<size_t>>* indices)
{
	auto new_indices = std::make_unique<std::vector<size_t>>(length);
	std::iota(new_indices->begin(), new_indices->end(), 0);

	auto val = std::stable_partition(new_indices->begin(), new_indices->end(), [begin, &flags](size_t index)
			{ return flags[index + begin] == 0; }) - new_indices->begin() + begin;

	if (indices)
	{
		*indices = std::move(new_indices);
	}

	return val;
}

std::unique_ptr<std::vector<size_t>> HostParticlePhaseSpace::stable_partition_alive(size_t begin, size_t length)
{
	if (length == static_cast<size_t>(-1))
	{
		length = n - begin;
	}

	std::unique_ptr<std::vector<size_t>> indices;
	n_alive = stable_partition_alive_indices(deathflags, begin, length, &indices);
	
	gather(r, *indices, begin, length);
	gather(v, *indices, begin, length);
	gather(deathtime, *indices, begin, length);
	gather(deathflags, *indices, begin, length);
	gather(id, *indices, begin, length);

	return indices;
}

Configuration::Configuration()
{
	t_0 = 0;
	t_f = 365e4;
	dt = 122;
	tbsize = 1024;
	wh_ce_n1 = 8;
	wh_ce_n2 = 4;
	wh_ce_r1 = 1;
	wh_ce_r2 = 3.5;

	print_every = 10;
	energy_every = 1;
	track_every = 0;
	dump_every = 1000;
	max_particle = static_cast<size_t>(-1);
	resolve_encounters = false;
	big_g = 1;
	icsin = "";
	plin = "";
	hybridin = "";
	hybridout = "state.out";
	readsplit = 0;
	writesplit = 0;
	dumpbinary = 1;
	trackbinary = 1;
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

bool read_configuration(std::istream& in, Configuration* out)
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
			std::cerr << "Unrecognized line " << linenum << std::endl;
			return true;
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
				out->tbsize = std::stoul(second);
			else if (first == "Big-G")
				out->big_g = std::stod(second);
			else if (first == "WH-Encounter-N1")
				out->wh_ce_n1 = std::stoull(second);
			else if (first == "WH-Encounter-N2")
				out->wh_ce_n2 = std::stoull(second);
			else if (first == "WH-Encounter-R1")
				out->wh_ce_r1 = std::stod(second);
			else if (first == "WH-Encounter-R2")
				out->wh_ce_r2 = std::stod(second);
			else if (first == "Limit-Particle-Count")
				out->max_particle = std::stoull(second);
			else if (first == "Log-Interval")
				out->print_every = std::stoull(second);
			else if (first == "Status-Interval")
				out->energy_every = std::stoull(second);
			else if (first == "Track-Interval")
				out->track_every = std::stoull(second);
			else if (first == "Dump-Interval")
				out->dump_every = std::stoull(second);
			else if (first == "Write-Binary-Dump")
				out->dumpbinary = std::stoi(second) != 0;
			else if (first == "Write-Binary-Track")
				out->trackbinary = std::stoi(second) != 0;
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
			std::cerr << "Unrecognized line " << linenum << std::endl;
			return true;
		}
	}

	if (out->outfolder.empty())
	{
		std::cerr << "Warning: Output-Folder was not specified, defaulting to ./" << std::endl;
		out->outfolder = "./";
	}
	if (out->dump_every == 0 && out->dumpbinary)
	{
		std::cerr << "Warning: Dumping is disabled but Write-Binary-Dump was specified" << std::endl;
	}
	if (out->track_every == 0 && out->trackbinary)
	{
		std::cerr << "Warning: Tracking is disabled but Write-Binary-Track was specified" << std::endl;
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
		std::cerr << "Error: Output-File was not specified" << std::endl;
		return true;
	}
	if (!out->writesplit && out->hybridin == "")
	{
		std::cerr << "Error: Input-File was not specified" << std::endl;
		return true;
	}
	if (out->readsplit && (out->plin == "" || out->icsin == ""))
	{
		std::cerr << "Error: Read-Split-Input was selected but Particle-Input-File or Planet-Input-File were not specified" << std::endl;
		return true;
	}

	return false;
}

void write_configuration(std::ostream& outstream, const Configuration& out)
{
	outstream << std::setprecision(17);
	outstream << "Initial-Time " << out.t_0 << std::endl;
	outstream << "Time-Step " << out.dt << std::endl;
	outstream << "Final-Time " << out.t_f << std::endl;
	outstream << "Time-Block-Size " << out.tbsize << std::endl;
	outstream << "WH-Encounter-N1 " << out.wh_ce_n1 << std::endl;
	outstream << "WH-Encounter-N2 " << out.wh_ce_n2 << std::endl;
	outstream << "WH-Encounter-R1 " << out.wh_ce_r1 << std::endl;
	outstream << "WH-Encounter-R2 " << out.wh_ce_r2 << std::endl;
	outstream << "Limit-Particle-Count " << out.max_particle << std::endl;
	outstream << "Log-Interval " << out.print_every << std::endl;
	outstream << "Status-Interval " << out.energy_every << std::endl;
	outstream << "Track-Interval " << out.track_every << std::endl;
	outstream << "Dump-Interval " << out.dump_every << std::endl;
	outstream << "Write-Binary-Track " << out.trackbinary << std::endl;
	outstream << "Write-Binary-Dump " << out.dumpbinary << std::endl;
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

std::string joinpath(const std::string& base, const std::string& file)
{
	return base + "/" + file;
}

bool load_data_nohybrid(HostData& hd, const Configuration& config, std::istream& plin, std::istream& icsin)
{
	size_t npl;
	plin >> npl;

	hd.planets = HostPlanetPhaseSpace(npl, config.tbsize, config.fast_timestep_factor());

	for (size_t i = 0; i < npl; i++)
	{
		plin >> hd.planets.m[i];
		plin >> hd.planets.r[i].x >> hd.planets.r[i].y >> hd.planets.r[i].z;
		plin >> hd.planets.v[i].x >> hd.planets.v[i].y >> hd.planets.v[i].z;

		hd.planets.id[i] = static_cast<uint32_t>(i);
	}

	size_t npart;
	icsin >> npart;
	npart = std::min(npart, config.max_particle);

	hd.particles = HostParticlePhaseSpace(npart);

	for (size_t i = 0; i < npart; i++)
	{
		icsin >> hd.particles.r[i].x >> hd.particles.r[i].y >> hd.particles.r[i].z;
		icsin >> hd.particles.v[i].x >> hd.particles.v[i].y >> hd.particles.v[i].z;

		std::string s;
		icsin >> s;
		if (!isdigit(s[0]))
		{
			icsin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			hd.particles.deathtime[i] = 0;
			hd.particles.id[i] = static_cast<uint32_t>(i);
			hd.particles.deathflags[i] = 0;
		}
		else
		{
			hd.particles.deathtime[i] = std::stof(s);
			icsin >> hd.particles.deathflags[i] >> hd.particles.id[i];
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
		ss >> hd.planets.m[i];

		std::getline(in, s);
		ss = std::istringstream(s);
		ss >> hd.planets.r[i].x >> hd.planets.r[i].y >> hd.planets.r[i].z;

		std::getline(in, s);
		ss = std::istringstream(s);
		ss >> hd.planets.v[i].x >> hd.planets.v[i].y >> hd.planets.v[i].z;

		std::getline(in, s);
		ss = std::istringstream(s);
		
		ss >> hd.planets.id[i];

		hd.planets.id[i] = static_cast<uint32_t>(i);
	}

	std::getline(in, s);
	ss = std::istringstream(s);
	size_t npart;
	ss >> npart;
	npart = std::min(npart, config.max_particle);

	hd.particles = HostParticlePhaseSpace(npart);

	for (size_t i = 0; i < npart; i++)
	{
		std::istringstream ss;
		std::string s;

		std::getline(in, s);
		ss = std::istringstream(s);
		ss >> hd.particles.r[i].x >> hd.particles.r[i].y >> hd.particles.r[i].z;

		std::getline(in, s);
		ss = std::istringstream(s);
		ss >> hd.particles.v[i].x >> hd.particles.v[i].y >> hd.particles.v[i].z;

		std::getline(in, s);
		ss = std::istringstream(s);
		ss >> hd.particles.id[i] >> hd.particles.deathflags[i] >> hd.particles.deathtime[i];
	}

	return false;
}

bool load_data_hybrid_binary(HostData& hd, const Configuration& config, std::istream& in)
{
	read_binary(in, hd.planets.n);
	hd.planets = HostPlanetPhaseSpace(hd.planets.n, config.tbsize, config.fast_timestep_factor());

	for (size_t i = 0; i < hd.planets.n; i++)
	{
		read_binary(in, hd.planets.id[i]);
		read_binary(in, hd.planets.m[i]);
		read_binary(in, hd.planets.r[i].x);
		read_binary(in, hd.planets.r[i].y);
		read_binary(in, hd.planets.r[i].z);
		read_binary(in, hd.planets.v[i].x);
		read_binary(in, hd.planets.v[i].y);
		read_binary(in, hd.planets.v[i].z);
	}

	read_binary(in, hd.particles.n);

	hd.particles.n = std::min(hd.particles.n, config.max_particle);
	hd.particles = HostParticlePhaseSpace(hd.particles.n);

	for (size_t i = 0; i < hd.particles.n; i++)
	{
		read_binary(in, hd.particles.id[i]);
		read_binary(in, hd.particles.r[i].x);
		read_binary(in, hd.particles.r[i].y);
		read_binary(in, hd.particles.r[i].z);
		read_binary(in, hd.particles.v[i].x);
		read_binary(in, hd.particles.v[i].y);
		read_binary(in, hd.particles.v[i].z);
		read_binary(in, hd.particles.deathflags[i]);
		read_binary(in, hd.particles.deathtime[i]);
	}

	return !in;
}

bool load_data(HostData& hd, const Configuration& config)
{
	bool ret;
	if (config.readsplit)
	{
		std::ifstream plinfile(config.plin), icsinfile(config.icsin);
		ret = load_data_nohybrid(hd, config, plinfile, icsinfile);
	}
	else
	{
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
		for (size_t i = 0; i < hd.planets.n; i++)
		{
			if (config.readmomenta)
			{
				hd.planets.v[i] /= hd.planets.m[i];
			}

			hd.planets.m[i] *= config.big_g;
		}

		hd.particles.stable_partition_alive();
	}

	return ret;
}

void save_data_hybrid_binary(const HostData& hd, const Configuration& config, std::ostream& out)
{
	(void) config;
	write_binary(out, hd.planets.n_alive);
	for (size_t i = 0; i < hd.planets.n_alive; i++)
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

	write_binary(out, hd.particles.n);
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		write_binary(out, hd.particles.id[i]);
		write_binary(out, hd.particles.r[i].x);
		write_binary(out, hd.particles.r[i].y);
		write_binary(out, hd.particles.r[i].z);
		write_binary(out, hd.particles.v[i].x);
		write_binary(out, hd.particles.v[i].y);
		write_binary(out, hd.particles.v[i].z);
		write_binary(out, hd.particles.deathflags[i]);
		write_binary(out, hd.particles.deathtime[i]);
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

	out << hd.particles.n << std::endl;
	out << std::setprecision(17);
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		out << hd.particles.r[i].x << " " << hd.particles.r[i].y << " " << hd.particles.r[i].z << std::endl;
		out << hd.particles.v[i].x << " " << hd.particles.v[i].y << " " << hd.particles.v[i].z << std::endl;
		out << hd.particles.id[i] << " " << hd.particles.deathflags[i] << " " << hd.particles.deathtime[i] << std::endl;
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

	icsout << hd.particles.n << std::endl;
	icsout << std::setprecision(17);
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		icsout << hd.particles.r[i].x << " " << hd.particles.r[i].y << " " << hd.particles.r[i].z << std::endl;
		icsout << hd.particles.v[i].x << " " << hd.particles.v[i].y << " " << hd.particles.v[i].z << std::endl;
		icsout << hd.particles.deathtime[i] << " " << hd.particles.deathflags[i] << " " << hd.particles.id[i] << std::endl;
	}
}

void save_data(const HostData& hd, const Configuration& config, const std::string& outfile, bool dump)
{
	if (dump || !config.writesplit)
	{
		std::ostringstream ss;

		if ((dump && config.dumpbinary) || (!config.writesplit && config.writebinary))
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
		std::ofstream ploutfile(joinpath(config.outfolder, "pl.out")), icsoutfile(joinpath(config.outfolder, "ics.out"));
		save_data_nohybrid(hd, config, ploutfile, icsoutfile);
	}
}
