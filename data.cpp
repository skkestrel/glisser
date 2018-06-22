#include "data.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <numeric>

template<typename T>
void gather(std::vector<T>& values, const std::vector<size_t>& indices)
{
	std::vector<T> copy(values.begin(), values.end());
	for (size_t i = 0; i < values.size(); i++)
	{
		values[i] = copy[indices[i]];
	}
}

void HostParticlePhaseSpace::stable_partition_alive()
{
	std::vector<size_t> indices(n);
	std::iota(indices.begin(), indices.end(), 0);
	n_alive = std::stable_partition(indices.begin(), indices.end(), [this](size_t index)
			{ return deathflags[index] == 0; }) - indices.begin();

	gather(r, indices);
	gather(v, indices);
	gather(deathtime, indices);
	gather(deathflags, indices);
	gather(id, indices);
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
	max_particle = 0;
	resolve_encounters = false;
	icsin = "";
	plin = "";
	outfolder = "output/";
	readmomenta = false;
	enable_ascii_track = 0;
	enable_binary_track = 0;
}


bool read_configuration(std::istream& in, Configuration* out)
{
	size_t linenum = 0;
	std::string line;
	while (std::getline(in, line))
	{
		if (line.length() == 0)
		{
			linenum++;
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
			else if (first == "Encounter-Slowdown-Factor")
				out->ce_factor = std::stoll(second);
			else if (first == "Limit-Particle-Count")
				out->max_particle = std::stoll(second);
			else if (first == "Status-Interval")
				out->print_every = std::stoll(second);
			else if (first == "Track-Interval")
				out->track_every = std::stoll(second);
			else if (first == "Dump-Interval")
				out->dump_every = std::stoll(second);
			else if (first == "Resolve-Encounters")
				out->resolve_encounters = std::stoi(second) != 0;
			else if (first == "Enable-Ascii-Track")
				out->enable_ascii_track = std::stoi(second) != 0;
			else if (first == "Enable-Binary-Track")
				out->enable_binary_track = std::stoi(second) != 0;
			else if (first == "Use-Hybrid-Output")
				out->writehybrid = std::stoi(second) != 0;
			else if (first == "Use-Hybrid-Input")
				out->readhybrid = std::stoi(second) != 0;
			else if (first == "Hybrid-Input")
				out->hybridin = second;
			else if (first == "Particle-Input")
				out->icsin = second;
			else if (first == "Planet-Input")
				out->plin = second;
			else if (first == "Output-Folder")
				out->outfolder = second;
			else if (first == "Read-Momenta")
				out->readmomenta = std::stoi(second) != 0;
			else
				throw std::invalid_argument("bad");
		}
		catch (std::invalid_argument)
		{
			std::cerr << "Unrecognized line " << linenum << std::endl;
			return true;
		}

		linenum++;
	}

	if (!out->resolve_encounters)
	{
		out->ce_factor = 1;
	}

	if (out->outfolder.empty())
	{
		out->outfolder = "./";
	}

	if (out->readhybrid && out->hybridin == "")
	{
		std::cerr << "Use-Hybrid-Input was selected but no hybrid input file was specified" << std::endl;
		return true;
	}
	if (!out->readhybrid && (out->plin == "" || out->icsin == ""))
	{
		std::cerr << "Use-Hybrid-Input was not selected but no input files was specified" << std::endl;
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
	outstream << "Encounter-Slowdown-Factor " << out.ce_factor << std::endl;
	outstream << "Limit-Particle-Count " << out.max_particle << std::endl;
	outstream << "Status-Interval " << out.print_every << std::endl;
	outstream << "Enable-Ascii-Track " << out.enable_ascii_track << std::endl;
	outstream << "Enable-Binary-Track " << out.enable_binary_track << std::endl;
	outstream << "Track-Interval " << out.track_every << std::endl;
	outstream << "Dump-Interval " << out.dump_every << std::endl;
	outstream << "Resolve-Encounters " << out.resolve_encounters << std::endl;
	outstream << "Use-Hybrid-Output " << out.writehybrid << std::endl;
	outstream << "Use-Hybrid-Input " << out.readhybrid << std::endl;
	outstream << "Hybrid-Input " << out.hybridin << std::endl;
	outstream << "Particle-Input " << out.icsin << std::endl;
	outstream << "Planet-Input " << out.plin << std::endl;
	outstream << "Output-Folder " << out.outfolder << std::endl;
	outstream << "Read-Momenta " << out.readmomenta << std::endl;
}

std::string joinpath(const std::string& base, const std::string& file)
{
	return base + "/" + file;
}

bool load_data(HostData& hd, const Configuration& config)
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
	if (config.max_particle > 0) npart = std::min(npart, config.max_particle);

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

	hd.particles.stable_partition_alive();

	return false;
}

void save_data(const HostData& hd, const Configuration& config, bool dump, size_t dumpnum)
{
	std::ofstream ploutfile(joinpath(config.outfolder, "pl.out")), icsoutfile(joinpath(config.outfolder, "ics.out"));

	ploutfile << hd.planets_snapshot.n << std::endl;
	ploutfile << std::setprecision(17);
	for (size_t i = 0; i < hd.planets_snapshot.n; i++)
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
