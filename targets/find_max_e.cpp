#include "../data.h"
#include "../util.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include "../cxxopts.h"
#pragma GCC diagnostic pop

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

struct ParticleInfo
{
	double emax;
	double emax_t;
};

int main(int argc, char** argv)
{
	cxxopts::Options options("find-max-e", "Find the maximum e and maximum time of e of each particle from tracks");
	options.add_options()
		("i,input", "Input file or directory", cxxopts::value<std::string>())
		("o,output", "Output file", cxxopts::value<std::string>());

	options.parse_positional({ "input", "output" });

	try
	{
		auto result = options.parse(argc, argv);

		if (result.count("i") == 0)
		{
			throw cxxopts::OptionException("Required argument -i");
		}
		if (result.count("o") == 0)
		{
			throw cxxopts::OptionException("Required argument -o");
		}

		std::string inpath = result["i"].as<std::string>();
		std::string outpath = result["o"].as<std::string>();

		std::ofstream outfile(outpath);

		std::unordered_map<uint32_t, ParticleInfo> map;
		sr::data::read_tracks(inpath, true, std::vector<uint32_t>(), true,
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

		for (auto& pair : map)
		{
			outfile << pair.first << " " << pair.second.emax << " " << pair.second.emax_t << std::endl;
		}
	}
	catch (cxxopts::OptionException& e)
	{
		std::cout << e.what() << std::endl;
		std::cout << options.help() << std::endl;
		return -1;
	}

	return 0;
}
