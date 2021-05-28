#include "../src/data.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

static const char USAGE[] = R"(find-max-e
Usage:
    find-max-e [options] <input> <output>

Options:
    -h, --help           Show this screen.
    -t <t>, --time <t>   Set the max time to read the track until
)";

struct ParticleInfo
{
	double emax;
	double emax_t;
};

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "find-max-e");

	try
	{
		std::string inpath = args["<input>"].asString();
		std::string outpath = args["<output>"].asString();

		std::unordered_map<uint32_t, ParticleInfo> map;

		sr::data::TrackReaderOptions opt;
		opt.take_all_planets = false;
		opt.take_all_particles = true;

		if (args["--time"])
		{
			opt.max_time = std::stod(args["--time"].asString());
		}

		sr::data::read_tracks(inpath, opt,
			[&](sr::data::HostPlanetSnapshot& pl, sr::data::HostParticleSnapshot& pa, double time)
			{
				(void) pl;

				for (size_t i = 0; i < pa.n; i++)
				{
					if (map.count(pa.id[i]) == 0 || pa.r[i].y > map[pa.id[i]].emax)
					{
						map[pa.id[i]].emax = pa.r[i].y;
						map[pa.id[i]].emax_t = time;
					}
				}
			});

		std::ofstream outfile(outpath);
		outfile << "id,emax,emax_t" << std::endl;
		for (auto& pair : map)
		{
			outfile << pair.first << "," << pair.second.emax << "," << pair.second.emax_t << std::endl;
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
