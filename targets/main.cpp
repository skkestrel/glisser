/*************************************************************

.---. .            .         .   .-,--.                        
\___  |- ,-. ,-. ,-| . . ,-. |-   `|__/ ,-. .  , ,-. ,-. . ,-. 
    \ |  ,-| |   | | | | `-. |    /  \  |-' | /  |-' |   | |-' 
`---' `' `-^ '   `-^ `-^ `-' `'   `-  ` `-' `'   `-' '   ' `-' 

*************************************************************/

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#include <execinfo.h>
#include <sys/stat.h>
#include <csignal>

#include "../executor_facade.h"
#include "../data.h"
#include "../wh.h"
#include "../convert.h"
#include "../util.h"

volatile sig_atomic_t end_loop = 0;

void term(int signum)
{
	(void) signum;
	end_loop = 1;
}

int main(int argv, char** argc)
{
	signal(SIGTERM, term);
	signal(SIGINT, term);

	if (sizeof(double) != 8)
	{
		std::cerr << "sizeof double is not 8 - this code will likely not work!" << std::endl;
	}

	if (sizeof(float) != 4)
	{
		std::cerr << "sizeof float is not 4 - this code will likely not work!" << std::endl;
	}

	std::string configin = "config.in";
	if (argv >= 2) configin = std::string(argc[1]);
	
	std::cout << "Reading from configuration file " << configin << std::endl;
	
	std::ifstream configfile(configin);

	Configuration config;
	if (read_configuration(configfile, &config)) return -1;
	

	{
		std::ofstream configstream(joinpath(config.outfolder, "config.in"));
		write_configuration(configstream, config);
	}

	mkdir(config.outfolder.c_str(), ACCESSPERMS);

	if (!is_dir_empty(config.outfolder))
	{
		std::cout << "Output folder is not empty!" << std::endl;
		return -1;
	}

	mkdir(joinpath(config.outfolder, "dump").c_str(), ACCESSPERMS);


	std::ofstream coutlog(joinpath(config.outfolder, "stdout"));
	teestream tout(std::cout, coutlog);

	tout << "Host uses little-endian floats? " << (is_float_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian doubles? " << (is_double_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian ints? " << (is_int_little_endian() ? "yes" : "no") << std::endl;

	HostData hd;
	ExecutorFacade ex(hd, tout);

	*ex.t = config.t_0;
	*ex.t_0 = config.t_0;
	*ex.dt = config.dt;
	*ex.t_f = config.t_f;
	*ex.tbsize = config.tbsize;
	*ex.resolve_encounters = config.resolve_encounters;
	*ex.ce_n1 = config.ce_n1;
	*ex.ce_n2 = config.ce_n2;

	if (load_data(hd, config)) return -1;

	std::ofstream trackout;
	if (config.trackbinary)
	{
		trackout = std::ofstream(joinpath(config.outfolder, "track.out"), std::ios_base::binary);
	}
	else
	{
		trackout = std::ofstream(joinpath(config.outfolder, "track.out"));
	}

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	std::ofstream encounterlog(joinpath(config.outfolder, "encounter.out"));
	*ex.encounter_output = &encounterlog;

	std::ofstream timelog(joinpath(config.outfolder, "time.out"));
	timelog << "start " << std::put_time(&tm, "%c %Z") << std::endl;

	ex.init();

	size_t counter = 0;
	size_t dump_num = 0;

	bool crashed = false;
	try
	{
		while (*ex.t < *ex.t_f)
		{
			ex.loop();
			counter++;

			ex.add_job([&timelog, &tout, &ex, &config, counter]()
				{
					bool output_energy = config.energy_every != 0 && (counter % config.energy_every == 0);
					bool log_out = config.print_every != 0 && (counter % config.print_every == 0);


					if (!log_out && !output_energy) return;

					double e_;
					f64_3 l_;
					calculate_planet_metrics(ex.hd.planets, &e_, &l_);
					double elapsed = ex.time();
					double total = elapsed * (*ex.t_f - *ex.t_0) / (*ex.t - *ex.t_0);

					if (output_energy)
					{
						timelog << "time " << elapsed << " " << ex.hd.particles.n_alive << " " << ex.hd.particles.n_encounter << std::endl;
						timelog << "ep " << e_ << std::endl;
						timelog << "lp " << l_ << std::endl;
					}

					if (log_out)
					{
						tout << "t=" << *ex.t << " (" << elapsed / total * 100 << "% " << elapsed << "m elapsed, "
							<< total << "m total " << total - elapsed << "m remain)" << std::endl;
						tout << "Error = " << (e_ - *ex.e_0) / *ex.e_0 * 100 << ", " <<
							ex.hd.particles.n_alive << " particles remaining, " << ex.hd.particles.n_encounter << " in encounter" << std::endl;
					}
				});
			
			bool dump = config.dump_every != 0 && counter % config.dump_every == 0;
			bool track = config.track_every != 0 && counter % config.track_every == 0;

			if (dump || track)
			{
				ex.download_data();

				if (dump)
				{
					ex.add_job([&tout, &ex, &config, &dump_num]()
						{
							tout << "Dumping to disk. t = " << *ex.t << std::endl;
							std::ostringstream ss;
							ss << "dump/config." << dump_num << ".out";

							config.t_f = *ex.t_f - *ex.t_0 + *ex.t;
							config.t_0 = *ex.t;
							std::ofstream configout(joinpath(config.outfolder, ss.str()));
							write_configuration(configout, config);

							ss = std::ostringstream();
							ss << "dump/state." << dump_num << ".out";
							save_data(ex.hd, config, joinpath(config.outfolder, ss.str()), true);
						});
				}

				if (track)
				{
					ex.add_job([&trackout, &ex, &config]()
						{
							if (config.trackbinary)
							{
								write_binary(trackout, *ex.t);
								write_binary(trackout, ex.hd.planets.n_alive - 1);

								for (size_t i = 1; i < ex.hd.planets.n_alive; i++)
								{
									int esign;
									double a, e, in, capom, om, f;
									to_elements(ex.hd.planets.m[i] + ex.hd.planets.m[0], ex.hd.planets.r[i], ex.hd.planets.v[i],
										&esign, &a, &e, &in, &capom, &om, &f);

									write_binary(trackout, static_cast<uint32_t>(ex.hd.planets.id[i]));
									write_binary(trackout, static_cast<float>(a));
									write_binary(trackout, static_cast<float>(e));
									write_binary(trackout, static_cast<float>(in));
								}

								write_binary(trackout, ex.hd.particles.n_alive);
								for (size_t i = 0; i < ex.hd.particles.n_alive; i++)
								{
									int esign;
									double a, e, in, capom, om, f;
									to_elements(ex.hd.planets.m[0], ex.hd.particles.r[i], ex.hd.particles.v[i],
										&esign, &a, &e, &in, &capom, &om, &f);

									write_binary(trackout, static_cast<uint32_t>(ex.hd.particles.id[i]));
									write_binary(trackout, static_cast<float>(a));
									write_binary(trackout, static_cast<float>(e));
									write_binary(trackout, static_cast<float>(in));
								}
							}
							else
							{
								trackout << std::setprecision(7);
								trackout << *ex.t << std::endl;
								trackout << ex.hd.planets.n_alive - 1 << std::endl;
								for (size_t i = 1; i < ex.hd.planets.n_alive; i++)
								{
									int esign;
									double a, e, in, capom, om, f;
									to_elements(ex.hd.planets.m[i] + ex.hd.planets.m[0], ex.hd.planets.r[i], ex.hd.planets.v[i],
										&esign, &a, &e, &in, &capom, &om, &f);
									trackout << ex.hd.planets.id[i] << " " << a << " " << e << " " << in << std::endl;
								}

								trackout << ex.hd.particles.n_alive << std::endl;
								for (size_t i = 0; i < ex.hd.particles.n_alive; i++)
								{
									int esign;
									double a, e, in, capom, om, f;
									to_elements(ex.hd.planets.m[0], ex.hd.particles.r[i], ex.hd.particles.v[i],
										&esign, &a, &e, &in, &capom, &om, &f);
									trackout << ex.hd.particles.id[i] << " " << a << " " << e << " " << in << std::endl;
								}
							}

						});
				}
			}


			if (end_loop)
			{
				tout << "Caught signal." << std::endl;
				throw std::exception();
			}
		}
	}
	catch (const std::exception& e)
	{
		void* array[50];
		size_t size = backtrace(array, 50);
		backtrace_symbols_fd(array, static_cast<int>(size), 2);

		tout << "Exception caught: " << std::endl;
		tout << e.what() << std::endl;
		tout << "Recovering data." << std::endl;
		crashed = true;
	}

	ex.finish();
	tout << "Saving to disk." << std::endl;
	save_data(hd, config, joinpath(config.outfolder, "state.out"));
	config.t_f = *ex.t_f - *ex.t_0 + *ex.t;
	config.t_0 = *ex.t;
	std::ofstream configout(joinpath(config.outfolder, "config.out"));
	write_configuration(configout, config);

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return crashed ? -1 : 0;
}
