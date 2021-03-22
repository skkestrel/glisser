/*************************************************************

GLISSE Integrator

*************************************************************/

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <cmath>
#include <iomanip>


#include <execinfo.h>
#include <csignal>

#include "../src/executor_facade.h"

#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

static const char USAGE[] = R"(sr(_cpu)
Usage:
    sr(_cpu) [options] [<config>]

Options:
    -h, --help         Show this screen.
)";

volatile sig_atomic_t end_loop = 0;

void term(int signum)
{
	if (signum == SIGINT)
	{
		end_loop = 1;
	}
}

int main(int argc, char** argv)
{
	std::ios_base::sync_with_stdio(false);

	if (sizeof(double) != 8)
	{
		std::cerr << "sizeof double is not 8 - this code will likely not work!" << std::endl;
	}

	if (sizeof(float) != 4)
	{
		std::cerr << "sizeof float is not 4 - this code will likely not work!" << std::endl;
	}

	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "sr");

	std::string configin = "config.in";
	if (args["<config>"]) configin = args["<config>"].asString();
	
	std::cout << "Reading from configuration file " << configin << std::endl;
	
	std::ifstream configfile(configin);

	sr::data::Configuration config_mut;

	try
	{
		read_configuration(configfile, &config_mut);
	}
	catch (std::exception& e)
	{
		std::cerr << "Could not read configuration." << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

	const sr::data::Configuration& config = config_mut;


	sr::util::make_dir(config.outfolder);

	if (!sr::util::is_dir_empty(config.outfolder))
	{
		std::cout << "Output folder is not empty! Do you want to continue?" << std::endl;
		std::cout << "Type \"Yes\" exactly as shown to continue: ";
	
		std::string s;
		std::getline(std::cin, s);

		if (s != "Yes") return -1;
	}

	{
		std::ofstream configstream(sr::util::joinpath(config.outfolder, "config.in"));
		write_configuration(configstream, config);
	}

	sr::util::make_dir(sr::util::joinpath(config.outfolder, "dumps"));
	sr::util::make_dir(sr::util::joinpath(config.outfolder, "tracks"));

	std::ofstream coutlog(sr::util::joinpath(config.outfolder, "stdout"));
	sr::util::teestream tout(std::cout, coutlog);

	tout << "Host uses little-endian floats? " << (sr::data::is_float_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian doubles? " << (sr::data::is_double_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian ints? " << (sr::data::is_int_little_endian() ? "yes" : "no") << std::endl;

	sr::data::HostData hd;

	sr::exec::ExecutorFacade ex(hd, config, tout);

	ex.t = config.t_0;

	if (load_data(hd.planets, hd.particles, config)) return -1;
	save_data(hd.planets.base, hd.particles, config, sr::util::joinpath(config.outfolder, "state.in"));

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	std::ofstream encounterlog(sr::util::joinpath(config.outfolder, "encounter.out"));
	ex.encounter_output = &encounterlog;

	std::ofstream timelog(sr::util::joinpath(config.outfolder, "time.out"));
	timelog << "start " << std::put_time(&tm, "%c %Z") << std::endl;

	ex.init();

	uint32_t counter = 0;
	uint32_t dump_num = 0;

	bool crashed = false;
	std::ofstream trackout;

	std::ofstream plhistout;

	signal(SIGTERM, term);
	signal(SIGINT, term);

	uint32_t track_num = 1;

	try
	{
		if (config.track_every != 0)
		{
			trackout = std::ofstream(sr::util::joinpath(config.outfolder, "tracks/track.0.out"), std::ios_base::binary);

			sr::data::HostParticleSnapshot snapshot_copy = ex.hd.particles.base;
			snapshot_copy.sort_by_id(0, snapshot_copy.n_alive);

			// TODO: make to_elements a user parameter.
			sr::data::save_binary_track(trackout, ex.hd.planets_snapshot, snapshot_copy, ex.t, config.write_ele_track, config.write_bary_track, config.write_single_track);
		}

		if (config.planet_hist_every != 0)
		{
			if (config.write_binary_hist)
			{
				plhistout = std::ofstream(sr::util::joinpath(config.outfolder, "plhist.out"), std::ios_base::binary);
				sr::data::save_binary_hist(plhistout, ex.hd.planets_snapshot, ex.t, config.write_single_hist);
			}
			else
			{
				plhistout = std::ofstream(sr::util::joinpath(config.outfolder, "plhist.out"));
				sr::data::save_txt_hist(plhistout, ex.hd.planets_snapshot, ex.t, config.write_single_hist);
			}
		}

		double cumulated_cputime, cumulated_gputime, cumulated_looptime;
		size_t cumulated_encounter;

		cumulated_cputime = cumulated_gputime = cumulated_looptime = 0;
		cumulated_encounter = 0;
		while (ex.t < config.t_f)
		{
			double cputime, gputime, looptime;
			size_t n_encounter;

			bool output_safe = ex.loop(&cputime, &gputime, &looptime, &n_encounter);
			cumulated_cputime += cputime;
			cumulated_gputime += gputime;
			cumulated_looptime += looptime;
			cumulated_encounter += n_encounter;
			// tout << cputime << " " << looptime << " " << n_encounter << std::endl;

			// if we are not in a safe time to output, skip everything
		    if (!output_safe)
			{
				continue;
			}

			counter++;
			ex.add_job([&timelog, &tout, &ex, &config, counter, cputime, gputime, looptime, n_encounter, &cumulated_cputime, &cumulated_gputime, &cumulated_looptime, &cumulated_encounter]()
				{	
					bool log_out = config.log_every != 0 && (counter % config.log_every == 0);
					if (!log_out) return;

					double elapsed = ex.time();
					double total = elapsed * (config.t_f - config.t_0) / (ex.t - config.t_0);

					if (log_out)
					{
						tout << std::setprecision(4) << std::endl;
						tout << "t=" << ex.t << " (" << std::setw(4) << elapsed / total * 100 << "%) | " << std::setw(6) << elapsed << "m elapsed | "
							 << std::setw(6) << total << "m total | " << std::setw(6) <<  total - elapsed << "m remain | ";
						tout << ex.hd.particles.n_alive() << " particles alive " << std::endl;

						tout << "GPU time (First)    : " << std::setw(6) << gputime                << " | CPU time (First):     " 
							 << std::setw(6) << cputime                 << " | LOOP time (First):     " 
							 << std::setw(6) << looptime                << " (ms)" << " | N encounters (First):     " 
							 << n_encounter << std::endl;

						tout << "GPU time (Cumulative): " << std::setw(6) << cumulated_gputime/1000 << " | CPU time (Cumulative): " 
							 << std::setw(6) << cumulated_cputime/1000  << " | LOOP time (Cumulative): " 
							 << std::setw(6) << cumulated_looptime/1000 << "  (s)" << " | N encounters (Cumulative): " 
							 << cumulated_encounter << std::endl;

						cumulated_cputime = cumulated_gputime = cumulated_looptime = 0;
						cumulated_encounter = 0;
					}
				});
			
			bool dump = config.dump_every != 0 && counter % config.dump_every == 0;
			bool track = config.track_every != 0 && counter % config.track_every == 0;
			bool plhist = config.planet_hist_every != 0 && counter % config.planet_hist_every == 0;

			if (dump || track || plhist)
			{
				// NOTE: this is not needed when we use a resync function that always downloads everything
				// ex.download_data();

				if (dump)
				{
					sr::data::Configuration out_config = config.output_config();
					// out_config.t_f = config.t_f - config.t_0 + ex.t;
					out_config.t_0 = ex.t;
					// out_config.writesplit = false;
					// out_config.writebinary = true;

					ex.add_job([&tout, &ex, &out_config, &config, &dump_num]()
						{
							tout << "Dumping to disk. t = " << ex.t << std::endl;
							std::ostringstream ss;
							ss << "dumps/config." << dump_num << ".out";

							std::ofstream configout(sr::util::joinpath(config.outfolder, ss.str()));
							write_configuration(configout, out_config);

							ss = std::ostringstream();
							ss << "dumps/state." << dump_num << ".out";
							save_data(ex.hd.planets_snapshot, ex.hd.particles, config, sr::util::joinpath(config.outfolder, ss.str()));

							dump_num++;
						});
				}

				if (track)
				{
					if (config.split_track_file > 0 && trackout.tellp() > config.split_track_file)
					{
						std::ostringstream ss;
						ss << "tracks/track." << track_num++ << ".out";
						trackout = std::ofstream(sr::util::joinpath(config.outfolder, ss.str()), std::ios_base::binary);
					}

					ex.add_job([&trackout, &ex, &config]()
						{
							sr::data::HostParticleSnapshot snapshot_copy = ex.hd.particles.base;
							snapshot_copy.sort_by_id(0, snapshot_copy.n_alive);
							sr::data::save_binary_track(trackout, ex.hd.planets_snapshot, snapshot_copy, ex.t, config.write_ele_track, config.write_bary_track, config.write_single_track);
						});
				}

				if (plhist)
				{
					ex.add_job([&plhistout, &ex, &config]()
						{
							if (config.write_binary_hist)
							{
								sr::data::save_binary_hist(plhistout, ex.hd.planets_snapshot, ex.t, config.write_single_hist);
							}
							else
							{
								sr::data::save_txt_hist(plhistout, ex.hd.planets_snapshot, ex.t, config.write_single_hist);
							}
						});
				}
			}

			if (end_loop)
			{
				tout << "Caught signal. What do you want to do?" << std::endl;
				tout << "dump <output_config> <output_file> | quit | continue" << std::endl;

				std::string s;
				while (std::getline(std::cin, s))
				{
					std::vector<std::string> tokens;
					std::stringstream ss(s);
					std::string token;
					while (ss >> token)
					{
						tokens.push_back(token);
					}

					if (tokens.size() < 1) continue;

					if (tokens[0] == "quit")
					{
						throw std::runtime_error("Quit requested by user");
					}
					else if (tokens[0] == "continue")
					{
						end_loop = 0;
						break;
					}
					else if (tokens[0] == "dump")
					{
						if (tokens.size() != 3)
						{
							tout << "?" << std::endl;
						}
						if (!dump && !track)
						{
							ex.download_data(0, ex.hd.particles.n());
						}

						// TODO dump on next 
						sr::data::Configuration out_config = config.output_config();
						out_config.t_0 = ex.t;

						tout << "Dumping to disk. t = " << ex.t << std::endl;
						std::ofstream configout(tokens[1]);
						write_configuration(configout, out_config);
						save_data(ex.hd.planets_snapshot, ex.hd.particles, config, tokens[2]);
					}
					else
					{
						tout << "?" << std::endl;
					}
				}
				// TODO add option to edit configuration
			}

		}

		ex.finish();
	}
	catch (const std::exception& e)
	{
		ex.finish();

		void* array[50];
		size_t size = backtrace(array, 50);
		backtrace_symbols_fd(array, static_cast<int>(size), 1);

		tout << "Exception caught: " << std::endl;
		tout << e.what() << std::endl;
		
		if (config.resolve_encounters)
		{
			tout << "End state may be corrupted. Resuming from a dump is recommended." << std::endl;
		}
		crashed = true;
	}

	ex.download_data(0, ex.hd.particles.n());

	// save last time into the history
	if (config.write_binary_hist)
	{
		sr::data::save_binary_hist(plhistout, ex.hd.planets_snapshot, ex.t, config.write_single_hist);
	}
	else
	{
		sr::data::save_txt_hist(plhistout, ex.hd.planets_snapshot, ex.t, config.write_single_hist);
	}

	tout << "Saving to disk." << std::endl;
	save_data(hd.planets_snapshot, hd.particles, config, sr::util::joinpath(config.outfolder, "state.out"));

	sr::data::Configuration out_config = config.output_config();
	// out_config.t_f = config.t_f - config.t_0 + ex.t;
	out_config.t_0 = ex.t;

	std::ofstream configout(sr::util::joinpath(config.outfolder, "config.out"));
	write_configuration(configout, out_config);

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return crashed ? -1 : 0;
}
