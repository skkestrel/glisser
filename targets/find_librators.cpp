#include "../src/data.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

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

	std::array<double, 21> window;
	double prev_minmax;
	int was_min_or_max;

	double center_sum;
	int center_ct;
	double first_amp;
	bool adv_ok;


	std::array<double, 41> window2;
	double prev_minmax2;
	int was_min_or_max2;

	double center_sum2;
	int center_ct2;
	double first_amp2;
	bool adv_ok2;

	ParticleInfo()
	{
		adv_ok = true;
		center_sum = 0;
		center_ct = 0;
		first_amp = -1;
		was_min_or_max = -1;

		adv_ok2 = true;
		center_sum2 = 0;
		center_ct2 = 0;
		first_amp2 = -1;
		was_min_or_max2 = -1;

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

			size_t index = 0;

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
						if (!map[pa.id[i]].ok && !map[pa.id[i]].okopp && !map[pa.id[i]].adv_ok && !map[pa.id[i]].adv_ok2)
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

						auto& ref = map[pa.id[i]];

						if (index < ref.window.size())
						{
							ref.window[index] = arg;
						}
						else
						{
							for (size_t j = 0; j < ref.window.size() - 1; j++)
							{
								ref.window[j] = ref.window[j + 1];
							if (!map[pa.id[i]].ok && !map[pa.id[i]].okopp && !map[pa.id[i]].adv_ok)
							{
								continue;
							}
							}
							ref.window[ref.window.size() - 1] = arg;

							// normalize window
							for (size_t j = 0; j < ref.window.size() - 1; j++)
							{
								while (std::abs(ref.window[j] - ref.window[j + 1]) > M_PI)
								{
									ref.window[j + 1] -= (ref.window[j + 1] - ref.window[j - 1] > 0) ? 2 * M_PI : -2 * M_PI;
								}
							}

							// check min/max
							bool max_ok = true;
							bool min_ok = true;
							double cur_minmax = ref.window[ref.window.size() / 2];

							for (size_t j = 0; j < ref.window.size(); j++)
							{
								if (j == ref.window.size() / 2) continue;

								if (ref.window[j] > cur_minmax)
								{
									max_ok = false;
								}
								else
								{
									min_ok = false;
								}
							}

							bool first_peak = false;
							bool prev_peak_correct = false;
							double amp = 0;
							double mid = 0;

							if (ref.was_min_or_max == -1) first_peak = true;
							if (max_ok) {
								prev_peak_correct = ref.was_min_or_max == 0 || first_peak;
								amp = cur_minmax - ref.prev_minmax;
								mid = ref.prev_minmax + amp;
							}
							if (min_ok) {
								prev_peak_correct = ref.was_min_or_max == 1 || first_peak;
								amp = ref.prev_minmax - cur_minmax;
								mid = cur_minmax + amp;
							}

							if (max_ok || min_ok)
							{
								if (prev_peak_correct && !first_peak)
								{
									while (amp < 0)
									{
										amp += 2 * M_PI;
									}
									while (amp > 2 * M_PI)
									{
										amp -= 2 * M_PI;
									}
									amp /= 2;

									if (ref.first_amp < 0)
									{
										ref.first_amp = amp;
									}
									else
									{
										// verify that amplitude and center are about the same
										if (std::abs(ref.first_amp - amp) > 0.5)
										{
											ref.adv_ok = false;
										}
									}

									// TODO need a circular mean here
									while (mid > 2 * M_PI)
									{
										mid -= 2 * M_PI;
									}
									while (mid < 0)
									{
										mid += 2 * M_PI;
									}
									ref.center_sum += mid;
									ref.center_ct += 1;
								}
								else if (!first_peak)
								{
									ref.adv_ok = false;
								}
								ref.prev_minmax = cur_minmax;

								if (max_ok) {
									ref.was_min_or_max = 1;
								}
								if (min_ok) {
									ref.was_min_or_max = 0;
								}
							}
						}


						if (index < ref.window2.size())
						{
							ref.window2[index] = arg;
						}
						else
						{
							for (size_t j = 0; j < ref.window2.size() - 1; j++)
							{
								ref.window2[j] = ref.window2[j + 1];
							}
							ref.window2[ref.window2.size() - 1] = arg;

							// normalize window
							for (size_t j = 0; j < ref.window2.size() - 1; j++)
							{
								while (std::abs(ref.window2[j] - ref.window2[j + 1]) > M_PI)
								{
									ref.window2[j + 1] -= (ref.window2[j + 1] - ref.window2[j - 1] > 0) ? 2 * M_PI : -2 * M_PI;
								}
							}

							// check min/max
							bool max_ok = true;
							bool min_ok = true;
							double cur_minmax = ref.window2[ref.window2.size() / 2];

							for (size_t j = 0; j < ref.window2.size(); j++)
							{
								if (j == ref.window2.size() / 2) continue;

								if (ref.window2[j] > cur_minmax)
								{
									max_ok = false;
								}
								else
								{
									min_ok = false;
								}
							}

							bool first_peak = false;
							bool prev_peak_correct = false;
							double amp = 0;
							double mid = 0;

							if (ref.was_min_or_max2 == -1) first_peak = true;
							if (max_ok) {
								prev_peak_correct = ref.was_min_or_max2 == 0 || first_peak;
								amp = cur_minmax - ref.prev_minmax2;
								mid = ref.prev_minmax2 + amp;
							}
							if (min_ok) {
								prev_peak_correct = ref.was_min_or_max2 == 1 || first_peak;
								amp = ref.prev_minmax2 - cur_minmax;
								mid = cur_minmax + amp;
							}

							if (max_ok || min_ok)
							{
								if (prev_peak_correct && !first_peak)
								{
									while (amp < 0)
									{
										amp += 2 * M_PI;
									}
									while (amp > 2 * M_PI)
									{
										amp -= 2 * M_PI;
									}
									amp /= 2;

									if (ref.first_amp2 < 0)
									{
										ref.first_amp2 = amp;
									}
									else
									{
										// verify that amplitude and center are about the same
										if (std::abs(ref.first_amp2 - amp) > 0.5)
										{
											ref.adv_ok2 = false;
										}
									}

									// TODO need a circular mean here
									while (mid > 2 * M_PI)
									{
										mid -= 2 * M_PI;
									}
									while (mid < 0)
									{
										mid += 2 * M_PI;
									}
									ref.center_sum2 += mid;
									ref.center_ct2 += 1;
								}
								else if (!first_peak)
								{
									ref.adv_ok2 = false;
								}
								ref.prev_minmax2 = cur_minmax;

								if (max_ok) {
									ref.was_min_or_max2 = 1;
								}
								if (min_ok) {
									ref.was_min_or_max2 = 0;
								}
							}
						}

						if (dist < tol)
						{
							map[pa.id[i]].ok = false;
						}

						if (dist_opp < tol)
						{
							map[pa.id[i]].okopp = false;
						}

					}
					index += 1;
				});
		}

		std::ofstream outfile(outpath);
		outfile << "id,lib,xlib" << std::endl;
		for (auto& pair : map)
		{
			outfile << pair.first << "," << (pair.second.ok && pair.second.alive) << "," <<
				(pair.second.okopp && !pair.second.ok && pair.second.alive) << ",";
			
			if (pair.second.adv_ok) {
				outfile << "1," << pair.second.center_sum / pair.second.center_ct << "," << pair.second.first_amp << std::endl;
			}
			else if (pair.second.adv_ok2) {
				outfile << "1," << pair.second.center_sum2 / pair.second.center_ct2 << "," << pair.second.first_amp2 << std::endl;
			}
			else {
				outfile << "0,0,0" << std::endl;
			}
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
