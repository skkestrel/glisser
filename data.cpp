#include "data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>

template<typename T>
void gather(std::vector<T>& values, const std::vector<size_t>& indices, size_t begin, size_t length)
{
	std::vector<T> copy(values.begin() + begin, values.begin() + begin + length);
	for (size_t i = begin; i < begin + length; i++)
	{
		values[i] = copy[indices[i - begin]];
	}
}

void HostParticlePhaseSpace::stable_partition_alive(size_t begin, size_t length)
{
	if (length == static_cast<size_t>(-1))
	{
		length = n - begin;
	}

	std::vector<size_t> indices(length);
	std::iota(indices.begin(), indices.end(), 0);
	n_alive = std::stable_partition(indices.begin(), indices.end(), [begin, this](size_t index)
			{ return deathflags[index + begin] == 0; }) - indices.begin() + begin;

	gather(r, indices, begin, length);
	gather(v, indices, begin, length);
	gather(deathtime, indices, begin, length);
	gather(deathflags, indices, begin, length);
	gather(id, indices, begin, length);
}

Configuration::Configuration()
{
	t_0 = 0;
	t_f = 365e4;
	dt = 122;
	tbsize = 1024;
	ce_factor = 8;
	print_every = 10;
	track_every = 0;
	dump_every = 1000;
	max_particle = static_cast<size_t>(-1);
	resolve_encounters = false;
	icsin = "";
	plin = "";
	hybridin = "";
	readhybrid = 0;
	writehybrid = 0;
	dumpbinary = 1;
	trackbinary = 0;
	writehybridbinary = 0;
	readhybridbinary = 0;
	outfolder = "output/";
	readmomenta = false;
}


bool read_configuration(std::istream& in, Configuration* out)
{
	size_t linenum = 0;
	std::string line;
	while (std::getline(in, line))
	{
		linenum++;
		if (line.length() == 0)
		{
			continue;
		}

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
				out->tbsize = std::stoll(second);
			else if (first == "Encounter-Time-Factor")
				out->ce_factor = std::stoll(second);
			else if (first == "Limit-Particle-Count")
				out->max_particle = std::stoll(second);
			else if (first == "Status-Interval")
				out->print_every = std::stoll(second);
			else if (first == "Track-Interval")
				out->track_every = std::stoll(second);
			else if (first == "Dump-Interval")
				out->dump_every = std::stoll(second);
			else if (first == "Write-Binary-Dump")
				out->dumpbinary = std::stoi(second) != 0;
			else if (first == "Write-Binary-Track")
				out->trackbinary = std::stoi(second) != 0;
			else if (first == "Resolve-Encounters")
				out->resolve_encounters = std::stoi(second) != 0;
			else if (first == "Write-Hybrid-Output")
				out->writehybrid = std::stoi(second) != 0;
			else if (first == "Read-Hybrid-Input")
				out->readhybrid = std::stoi(second) != 0;
			else if (first == "Write-Hybrid-Binary-Output")
				out->writehybridbinary = std::stoi(second) != 0;
			else if (first == "Read-Hybrid-Binary-Input")
				out->readhybridbinary = std::stoi(second) != 0;
			else if (first == "Hybrid-Input-File")
				out->hybridin = second;
			else if (first == "Particle-Input-File")
				out->icsin = second;
			else if (first == "Planet-Input-File")
				out->plin = second;
			else if (first == "Output-Folder")
				out->outfolder = second;
			else if (first == "Read-Input-Momenta")
				out->readmomenta = std::stoi(second) != 0;
			else
				throw std::invalid_argument("bad");
		}
		catch (std::invalid_argument)
		{
			std::cerr << "Unrecognized line " << linenum << std::endl;
			return true;
		}
	}

	if (!out->resolve_encounters && out->ce_factor != 1)
	{
		std::cerr << "Warning: Resolve-Encounters was not selected but Encounter-Time-Factor is not 1" << std::endl;
		out->ce_factor = 1;
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
	if (!out->readhybrid && out->readhybridbinary)
	{
		std::cerr << "Warning: Read-Hybrid-Input was not selected but Read-Hybrid-Binary-Input was specified" << std::endl;
	}
	if (!out->writehybrid && out->writehybridbinary)
	{
		std::cerr << "Warning: Write-Hybrid-Input was not selected but Write-Hybrid-Binary-Input was specified" << std::endl;
	}

	if (out->readhybrid && out->hybridin == "")
	{
		std::cerr << "Error: Read-Hybrid-Input was selected but no hybrid input file was specified" << std::endl;
		return true;
	}
	if (!out->readhybrid && (out->plin == "" || out->icsin == ""))
	{
		std::cerr << "Error: Read-Hybrid-Input was not selected but no input files was specified" << std::endl;
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
	outstream << "Encounter-Time-Factor " << out.ce_factor << std::endl;
	outstream << "Limit-Particle-Count " << out.max_particle << std::endl;
	outstream << "Status-Interval " << out.print_every << std::endl;
	outstream << "Track-Interval " << out.track_every << std::endl;
	outstream << "Dump-Interval " << out.dump_every << std::endl;
	outstream << "Write-Binary-Track " << out.trackbinary << std::endl;
	outstream << "Write-Binary-Dump" << out.dumpbinary << std::endl;
	outstream << "Resolve-Encounters " << out.resolve_encounters << std::endl;
	outstream << "Write-Hybrid-Output " << out.writehybrid << std::endl;
	outstream << "Write-Hybrid-Binary-Output " << out.writehybridbinary << std::endl;
	outstream << "Read-Hybrid-Input " << out.readhybrid << std::endl;
	outstream << "Read-Hybrid-Binary-Input " << out.readhybridbinary << std::endl;
	outstream << "Hybrid-Input-File " << out.hybridin << std::endl;
	outstream << "Particle-Input-File " << out.icsin << std::endl;
	outstream << "Planet-Input-File " << out.plin << std::endl;
	outstream << "Output-Folder " << out.outfolder << std::endl;
	outstream << "Read-Input-Momenta " << out.readmomenta << std::endl;
}

std::string joinpath(const std::string& base, const std::string& file)
{
	return base + "/" + file;
}

bool load_data_nohybrid(HostData& hd, const Configuration& config)
{
	std::ifstream plinfile(config.plin), icsinfile(config.icsin);

	size_t npl;
	plinfile >> npl;

	hd.planets = HostPlanetPhaseSpace(npl, config.tbsize, config.ce_factor);

	for (size_t i = 0; i < npl; i++)
	{
		plinfile >> hd.planets.m[i];
		plinfile >> hd.planets.r[i].x >> hd.planets.r[i].y >> hd.planets.r[i].z;
		plinfile >> hd.planets.v[i].x >> hd.planets.v[i].y >> hd.planets.v[i].z;

		if (config.readmomenta)
		{
			hd.planets.v[i].x /= hd.planets.m[i];
			hd.planets.v[i].y /= hd.planets.m[i];
			hd.planets.v[i].z /= hd.planets.m[i];
		}

		hd.planets.id[i] = i;
	}

	size_t npart;
	icsinfile >> npart;
	npart = std::min(npart, config.max_particle);

	hd.particles = HostParticlePhaseSpace(npart);

	for (size_t i = 0; i < npart; i++)
	{
		icsinfile >> hd.particles.r[i].x >> hd.particles.r[i].y >> hd.particles.r[i].z;
		icsinfile >> hd.particles.v[i].x >> hd.particles.v[i].y >> hd.particles.v[i].z;

		std::string s;
		icsinfile >> s;
		if (!isdigit(s[0]))
		{
			icsinfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			hd.particles.deathtime[i] = 0;
			hd.particles.id[i] = i;
			hd.particles.deathflags[i] = 0;
		}
		else
		{
			hd.particles.deathtime[i] = std::stod(s);
			icsinfile >> hd.particles.deathflags[i] >> hd.particles.id[i];
		}
	}

	return false;
}

bool load_data_hybrid(HostData& hd, const Configuration& config)
{
	std::ifstream in(config.hybridin);

	size_t npl;
	in >> npl;

	hd.planets = HostPlanetPhaseSpace(npl, config.tbsize, config.ce_factor);

	for (size_t i = 0; i < npl; i++)
	{
		std::string s;
		std::getline(in, s);
		std::istringstream ss(s);
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

		if (config.readmomenta)
		{
			hd.planets.v[i] /= hd.planets.m[i];
		}

		hd.planets.id[i] = i;
	}

	size_t npart;
	in >> npart;
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

bool load_data_hybrid_binary(HostData& hd, const Configuration& config)
{
	std::ifstream in(config.hybridin, std::ios_base::binary);
	read_binary(in, hd.planets.n);
	hd.planets = HostPlanetPhaseSpace(hd.planets.n, config.tbsize, config.ce_factor);

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

		if (config.readmomenta)
		{
			hd.planets.v[i] /= hd.planets.m[i];
		}
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
	if (!config.readhybrid)
	{
		ret = load_data_nohybrid(hd, config);
	}
	else
	{
		if (config.readhybridbinary)
		{
			ret = load_data_hybrid_binary(hd, config);
		}
		else
		{
			ret = load_data_hybrid(hd, config);
		}
	}

	if (!ret) hd.particles.stable_partition_alive();
	return ret;
}

void save_data_hybrid_binary(const HostData& hd, std::ostream& out, const Configuration& config)
{
	(void) config;
	write_binary(out, hd.planets.n_alive);
	for (size_t i = 0; i < hd.planets.n_alive; i++)
	{
		write_binary(out, hd.planets.id[i]);
		write_binary(out, hd.planets.m[i]);
		write_binary(out, hd.planets.r[i].x);
		write_binary(out, hd.planets.r[i].y);
		write_binary(out, hd.planets.r[i].z);
		write_binary(out, hd.planets.v[i].x);
		write_binary(out, hd.planets.v[i].y);
		write_binary(out, hd.planets.v[i].z);
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

void save_data_hybrid(const HostData& hd, std::ostream& out, const Configuration& config)
{
	(void) config;

	out << hd.planets_snapshot.n_alive << std::endl;
	out << std::setprecision(17);
	for (size_t i = 0; i < hd.planets_snapshot.n_alive; i++)
	{
		out << hd.planets_snapshot.m[i] << std::endl;
		out << hd.planets_snapshot.r[i].x << " " << hd.planets_snapshot.r[i].y << " " << hd.planets_snapshot.r[i].z << std::endl;
		out << hd.planets_snapshot.v[i].x << " " << hd.planets_snapshot.v[i].y << " " << hd.planets_snapshot.v[i].z << std::endl;
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

void save_data_nohybrid(const HostData& hd, const Configuration& config)
{
	(void) config;
	std::ofstream ploutfile(joinpath(config.outfolder, "pl.out")), icsoutfile(joinpath(config.outfolder, "ics.out"));

	ploutfile << hd.planets_snapshot.n_alive << std::endl;
	ploutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.planets_snapshot.n_alive; i++)
	{
		ploutfile << hd.planets_snapshot.m[i] << std::endl;
		ploutfile << hd.planets_snapshot.r[i].x << " " << hd.planets_snapshot.r[i].y << " " << hd.planets_snapshot.r[i].z << std::endl;
		ploutfile << hd.planets_snapshot.v[i].x << " " << hd.planets_snapshot.v[i].y << " " << hd.planets_snapshot.v[i].z << std::endl;
	}

	icsoutfile << hd.particles.n << std::endl;
	icsoutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		icsoutfile << hd.particles.r[i].x << " " << hd.particles.r[i].y << " " << hd.particles.r[i].z << std::endl;
		icsoutfile << hd.particles.v[i].x << " " << hd.particles.v[i].y << " " << hd.particles.v[i].z << std::endl;
		icsoutfile << hd.particles.deathtime[i] << " " << hd.particles.deathflags[i] << " " << hd.particles.id[i] << std::endl;
	}
}

void save_data(const HostData& hd, const Configuration& config, bool dump, size_t dumpnum)
{
	if (dump || config.writehybrid)
	{
		std::ostringstream ss;

		if (dump)
		{
			ss << "dump/state." << dumpnum << ".out";
		}
		else
		{
			ss << "state.out";
		}

		if ((dump && config.dumpbinary) || (!dump && config.writehybridbinary))
		{
			std::ofstream out(joinpath(config.outfolder, ss.str()), std::ios_base::binary);
			save_data_hybrid_binary(hd, out, config);
		}
		else
		{
			std::ofstream out(joinpath(config.outfolder, ss.str()));
			save_data_hybrid(hd, out, config);
		}
	}

	save_data_nohybrid(hd, config);
}
