#include "../src/data.h"
#include "../src/util.h"
#include "../src/convert.h"
#include "../docopt/docopt.h"

#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

static const char USAGE[] = R"(cj-metrics
Usage:
    cj_metrics [options] <input> <output>

Options:
    -h, --help           Show this screen.
)";

struct ParticleInfo
{
	double init_cj, final_cj;
	double init_a, final_a;
	double max_cj_error;
};

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "cj-metrics");

	try
	{
		std::string inpath = args["<input>"].asString();
		std::string outpath = args["<output>"].asString();

		std::unordered_map<uint32_t, ParticleInfo> map;

		sr::data::TrackReaderOptions opt;
		opt.take_all_planets = true;
		opt.take_all_particles = true;

		bool first_run = true;
		double m_n = 2.0438420219676245e-3;
		double m_s = 4 * M_PI * M_PI;

		size_t planet_index = 1;

		sr::data::read_tracks(inpath, opt,
			[&](sr::data::HostPlanetSnapshot& pl, sr::data::HostParticleSnapshot& pa, double time)
			{
				double a_n = pl.r[planet_index].x; f64_3 r_n, v_n;

				sr::convert::from_elements(m_s + m_n, pl.r[planet_index].x, pl.r[planet_index].y, pl.r[planet_index].z, pl.v[planet_index].x, pl.v[planet_index].y, pl.v[planet_index].z, &r_n, &v_n);

				double mu = m_n / m_s;
				double factor = -m_n / (m_s + m_n);
				f64_3 r_s = r_n * factor;

				for (size_t i = 0; i < pa.n; i++)
				{
					f64_3 r, v;
					sr::convert::from_elements(m_s + m_n, pa.r[i].x, pa.r[i].y, pa.r[i].z, pa.v[i].x, pa.v[i].y, pa.v[i].z, &r, &v);

					double R = std::sqrt(r.lensq()) / a_n;
					double R1 = std::sqrt((r - r_s).lensq()) / a_n;
					double R2 = std::sqrt((r - r_n).lensq()) / a_n;


					double cj = pl.r[planet_index].x / pa.r[i].x + 2 * std::sqrt(pa.r[i].x / pl.r[planet_index].x * (1 - pa.r[i].y * pa.r[i].y)) * std::cos(pa.r[i].z) - 2 / R + 2 * (1 - mu) / R1 + 2 * mu / R2;

					if (first_run)
					{
						map[pa.id[i]].init_cj = cj;
						map[pa.id[i]].init_a = pa.r[i].x;
						map[pa.id[i]].max_cj_error = 0;
					}
					else
					{
						double err = std::abs(map[pa.id[i]].init_cj - cj) / map[pa.id[i]].init_cj;
						if (map[pa.id[i]].max_cj_error < err) map[pa.id[i]].max_cj_error = err;
					}

					map[pa.id[i]].final_a = pa.r[i].x;
					map[pa.id[i]].final_cj = cj;
				}
				first_run = false;
			});

		std::ofstream outfile(outpath);
		outfile << "id,cj_error_f,max_cj_error,da" << std::endl;
		for (auto& pair : map)
		{
			outfile << pair.first << "," << std::abs(pair.second.final_cj - pair.second.init_cj) / pair.second.init_cj << "," << pair.second.max_cj_error << "," << std::abs(pair.second.final_a - pair.second.init_a) << std::endl;
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
