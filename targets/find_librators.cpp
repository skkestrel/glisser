#include "../src/data.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

static const char USAGE[] = R"(find-librators
Usage:
    find-librators [options] --mmr <a> <input> <output>

Options:
    -h, --help           Show this screen.
    -t <t>, --time <t>   Set the max time
    --slow               Use slow mode
    --tolerance <deg>    Tolerance [default: 10]
    --center <deg>       Libration center [default: 180]
    --mmr <a>            mmr symbol (example: '3:1@4')
)";

struct ParticleInfo
{
	bool ok;
	bool okopp;
	bool alive;

	ParticleInfo()
	{
		ok = true;
		okopp = true;
	}
};

using mmr_t = std::tuple<int32_t, int32_t, uint32_t>;

double mean_anomaly(double e, double f)
{
	double E = std::acos((e + std::cos(f)) / (1 + e * std::cos(f)));
	E = std::copysign(E, f);
	return E - e * std::sin(E);
}

double principal_angle(double t)
{
        t = t - 2 * M_PI * std::round(t / (2 * M_PI));
        if (t < -M_PI)
            t += 2 * M_PI;
        if (t > M_PI)
            t -= 2 * M_PI;
	return t;
}

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "find-librators");

	try
	{
		std::string inpath = args["<input>"].asString();
		std::string outpath = args["<output>"].asString();
		std::string mmrstring = args["--mmr"].asString();

		mmr_t mmr;

		{
			std::stringstream ss(mmrstring);
			std::string token;

			std::getline(ss, token, ':');
			std::get<0>(mmr) = std::stoi(token);
			std::getline(ss, token, '@');
			std::get<1>(mmr) = std::stoi(token);
			std::getline(ss, token, '\0');
			std::get<2>(mmr) = std::stoi(token);
		}

		std::unordered_map<uint32_t, ParticleInfo> map;

		sr::data::TrackReaderOptions opt;
		opt.take_all_planets = true;
		opt.take_all_particles = true;

		if (args["--time"])
		{
			opt.max_time = std::stod(args["--time"].asString());
		}

    		if (args["--slow"].asBool())
		{
			throw std::runtime_error("not supported");
		}
		else
		{
			double tol = std::stod(args["--tolerance"].asString()) / 180 * M_PI;
			double ctr = std::stod(args["--center"].asString()) / 180 * M_PI;
			double opp = principal_angle(ctr + M_PI);

			sr::data::read_tracks(inpath, opt,
				[&](sr::data::HostPlanetSnapshot& pl, sr::data::HostParticleSnapshot& pa, double time)
				{
					for (auto& pair : map)
					{
						pair.second.alive = false;
					}

					int planet_index = -1;
					for (size_t i = 0; i < pl.n; i++)
					{
						if (std::get<2>(mmr) == pl.id[i])
						{
							planet_index = static_cast<int>(i);
						}
					}

					if (planet_index < 0)
					{
						std::ostringstream ss;
						ss << "planet not found at time " << time;
						std::runtime_error(ss.str());
					}

					double pl_e = pl.r[planet_index].y;
					double pl_f = pl.v[planet_index].z;
					double pl_O = pl.v[planet_index].x;
					double pl_o = pl.v[planet_index].y;
					double pl_M = mean_anomaly(pl_e, pl_f);

					for (size_t i = 0; i < pa.n; i++)
					{
						map[pa.id[i]].alive = true;
						if (!map[pa.id[i]].ok && !map[pa.id[i]].okopp)
						{
							continue;
						}

						double pa_e = pa.r[i].y;
						double pa_f = pa.v[i].z;
						double pa_O = pa.v[i].x;
						double pa_o = pa.v[i].y;
						double pa_M = mean_anomaly(pa_e, pa_f);

						double arg = std::get<0>(mmr) * (pa_O + pa_o + pa_M) - std::get<1>(mmr) * (pl_O + pl_o + pl_M) +
							(std::get<1>(mmr) - std::get<0>(mmr)) * (pa_O + pa_o);
						arg = principal_angle(arg);

						double dist = std::min((2 * M_PI) - std::abs(arg - opp), std::abs(arg - opp));
						double dist_opp = std::min((2 * M_PI) - std::abs(arg - ctr), std::abs(arg - ctr));

						if (dist < tol)
						{
							map[pa.id[i]].ok = false;
						}

						if (dist_opp < tol)
						{
							map[pa.id[i]].okopp = false;
						}
					}
				});
		}

		std::ofstream outfile(outpath);
		outfile << "id,lib,xlib" << std::endl;
		for (auto& pair : map)
		{
			outfile << pair.first << "," << (pair.second.ok && pair.second.alive) << "," <<
				(pair.second.okopp && !pair.second.ok && pair.second.alive) << std::endl;
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
