#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

static const char USAGE[] = R"(import-swift
Usage:
    import-swift [options] <plin> <tpin> <outstate>

Options:
    -h, --help         Show this screen.
    -s <val>           Set the length of the istat and rstat arrays. [default: 13]
)";

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "import-swift");

	std::string pl = args["<plin>"].asString();
	std::string tp = args["<tpin>"].asString();

	std::string out = args["<outstate>"].asString();

	size_t statlen = std::stoi(args["-s"].asString());

	sr::data::HostData hd;

	std::ifstream plin(pl);
	std::ifstream tpin(tp);

	size_t ntp, npl;
	plin >> npl;
	std::cout << "npl " << npl << std::endl;
	hd.planets = sr::data::HostPlanetPhaseSpace(npl);

	for (size_t i = 0; i < npl; i++)
	{
		plin >> hd.planets.m()[i];

		// hd.planets.m()[i] /= 365.24 * 365.24;

		if (i == 0)
		{
			double r;
			plin >> r;
		}

		plin >> hd.planets.r()[i].x >> hd.planets.r()[i].y >> hd.planets.r()[i].z;
		plin >> hd.planets.v()[i].x >> hd.planets.v()[i].y >> hd.planets.v()[i].z;
		// hd.planets.v()[i] = hd.planets.v()[i] * (1 / 365.24);
		hd.planets.id()[i] = i;
	}

	tpin >> ntp;
	hd.particles = sr::data::HostParticlePhaseSpace(ntp);

	for (size_t i = 0; i < ntp; i++)
	{
		tpin >> hd.particles.r()[i].x >> hd.particles.r()[i].y >> hd.particles.r()[i].z;
		tpin >> hd.particles.v()[i].x >> hd.particles.v()[i].y >> hd.particles.v()[i].z;
		// hd.particles.v()[i] = hd.particles.v()[i] * (1 / 365.24);
		hd.particles.id()[i] = i;

		for (size_t j = 0; j < statlen; j++)
		{
			int istat;
			tpin >> istat;

			if (j == 0 && istat == 1)
			{
				hd.particles.deathflags()[i] = 0x81;
			}

			if (j == 1 && hd.particles.deathflags()[i])
			{
				if (istat == -1)
				{
					hd.particles.deathflags()[i] |= 0x0004;
				}
				if (istat == -2 || istat == -3)
				{
					hd.particles.deathflags()[i] |= 0x0002;
				}

				if (istat > 1)
				{
					hd.particles.deathflags()[i] |= ((istat - 1) << 8);
				}
			}
		}

		for (size_t j = 0; j < statlen; j++)
		{
			double rstat;
			tpin >> rstat;

			if (j == 0 && hd.particles.deathflags()[i])
			{
				hd.particles.deathtime_map()[i] = rstat;
			}
		}
	}

	save_data(hd.planets.base, hd.particles, sr::data::Configuration::create_dummy(), out);
}
