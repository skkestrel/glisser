#include "../src/data.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

static const char USAGE[] = R"(drive-res-amp
Usage:
    drive-res-amp [options] <input> <output>

Options:
    -h, --help                     Show this screen.
    -w <list>, --watch <list>      Take only the comma-separated list of particles, or "all" [default: none]
    -p <list>, --planets <list>    Take the comma-separated list of particles, or "none" [default: all]
    -l <n>, --split <n>            Split output every n bytes
    -s <n>, --skip <n>             Take every n time steps [default: 1]
    -t <t>, --tmax <t>             Take only up to given time
)";

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "prune-track");

	try
	{
		std::vector<uint32_t> particles;
		bool takeallparticles = false;
		if (args["--watch"].asString() == "all")
		{
			takeallparticles = true;
		}
		else if (args["--watch"].asString() == "none")
		{
			takeallparticles = false;
		}
		else
		{
			std::istringstream ss(args["--watch"].asString());
			std::string token;

			while (std::getline(ss, token, ','))
			{
				particles.push_back(static_cast<uint32_t>(std::stoul(token)));
			}
		}

		std::vector<uint32_t> planet_filter;
		bool takeallplanets = true;
		if (args["--planets"].asString() == "all")
		{
			takeallplanets = true;
		}
		else if (args["--planets"].asString() == "none")
		{
			takeallplanets = false;
		}
		else
		{
			std::istringstream ss(args["--planets"].asString());
			std::string token;

			while (std::getline(ss, token, ','))
			{
				planet_filter.push_back(static_cast<uint32_t>(std::stoul(token)));
			}
		}


		std::sort(particles.begin(), particles.end());

		int64_t splitbytes = 0;
		if (args["--split"])
		{
			splitbytes = args["split"].asLong();

			if (splitbytes < 0)
			{
				throw std::runtime_error("Cannot specify a negative integer for split");
			}
		}

		std::string inpath = args["<input>"].asString();
		std::string outpath = args["<output>"].asString();

		if (splitbytes > 0 && !sr::util::is_dir_empty(outpath))
		{
			std::cout << "Output folder is not empty! Do you want to continue?" << std::endl;
			std::cout << "Type \"Yes\" exactly as shown to continue: ";
		
			std::string s;
			std::getline(std::cin, s);

			if (s != "Yes") return -1;
		}

		std::ostringstream ss;

		if (splitbytes > 0)
		{
			sr::util::make_dir(outpath);

			ss << outpath;
			if (outpath[outpath.size() - 1] != '/') ss << '/';
			ss << "track.0.out";
		}
		else
		{
			ss << outpath;
		}

		size_t outnum = 1;
		std::ofstream outfile(ss.str(), std::ios_base::binary);

		sr::data::TrackReaderOptions opt;
		opt.take_all_particles = takeallparticles;
		opt.particle_filter = std::move(particles);
		opt.take_all_planets = takeallplanets;
		opt.planet_filter = std::move(planet_filter);

		sr::data::read_tracks(inpath, opt,
			[&](sr::data::HostPlanetSnapshot& pl, sr::data::HostParticleSnapshot& pa, double time)
			{
				sr::data::save_binary_track(outfile, pl, pa, time, false, false);

				if (splitbytes != 0 && outfile.tellp() > static_cast<int>(splitbytes))
				{
					ss = std::ostringstream();

					ss << outpath;
					if (outpath[outpath.size() - 1] != '/') ss << '/';
					ss << "track." << outnum++ << ".out";

					outfile = std::ofstream(ss.str(), std::ios_base::binary);
				}
			});
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
