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
#include <sys/stat.h>
#include <thread>
#include <iomanip>

#include "types.cuh"
#include "executor.cuh"
#include "data.h"
#include "wh.h"
#include "convert.h"
#include "util.h"

#include <execinfo.h>
#include <csignal>

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

	std::string configin = "config.in";
	if (argv >= 2) configin = std::string(argc[1]);
	
	std::cout << "Reading from configuration file " << configin << std::endl;
	
	std::ifstream configfile(configin);

	Configuration config;
	if (read_configuration(configfile, &config)) return -1;
	

	{
		std::ofstream configin(joinpath(config.outfolder, "config.in"));
		write_configuration(configin, config);
	}

	mkdir(config.outfolder.c_str(), ACCESSPERMS);
	mkdir(joinpath(config.outfolder, "dumps").c_str(), ACCESSPERMS);

	std::ofstream coutlog(joinpath(config.outfolder, "stdout"));
	teestream tout(std::cout, coutlog);

	tout << "Host uses little-endian floats? " << (is_float_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian doubles? " << (is_double_little_endian() ? "yes" : "no") << std::endl;
	tout << "Host uses little-endian ints? " << (is_int_little_endian() ? "yes" : "no") << std::endl;

	HostData hd;
	DeviceData dd;
	Executor ex(hd, dd, tout);

	ex.t = config.t;
	ex.t_0 = config.t;
	ex.dt = config.dt;
	ex.t_f = config.t_f;
	ex.tbsize = config.tbsize;
	ex.resolve_encounters = config.resolve_encounters;
	ex.ce_factor = config.ce_factor;
	ex.print_every = config.print_every;

	if (load_data(hd, config)) return -1;

	std::ofstream timelog(joinpath(config.outfolder, "time.out"));
	std::ofstream disclog(joinpath(config.outfolder, "discard.out"));

	std::ofstream track, trackb;
	if (config.enable_ascii_track) track = std::ofstream(joinpath(config.outfolder, "track.out"));
	if (config.enable_binary_track) trackb = std::ofstream(joinpath(config.outfolder, "trackb.out"), std::ios_base::binary);

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	ex.timing_output = &timelog;
	ex.discard_output = &disclog;

	timelog << "start " << std::put_time(&tm, "%c %Z") << std::endl;

	ex.init();

	size_t counter = 0;
	size_t dump_num = 0;

	bool crashed = false;
	try
	{
		while (ex.t < ex.t_f)
		{
			ex.loop();
			counter++;

			if (config.dump_every != 0 && counter % config.dump_every == 0)
			{
				ex.download_data();

				ex.add_job([&]()
					{
						tout << "Dumping to disk. t = " << ex.t << std::endl;
						std::ostringstream ss;
						ss << "dump/config." << dump_num << ".out" << std::endl;

						save_data(hd, config, true, dump_num++);

						config.t_f = ex.t_f - ex.t_0 + ex.t;
						config.t_0 = ex.t;
						std::ofstream configout(joinpath(config.outfolder, ss.str()));
						write_configuration(configout, config);
					});
			}

			if (config.track_every != 0 && counter % config.track_every == 0)
			{
				ex.download_data();

				ex.add_job([&]()
					{
						track << std::setprecision(7);
						track << ex.t << std::endl;
						track << hd.planets.n_alive - 1 << std::endl;

						write_binary(trackb, ex.t);
						write_binary(trackb, hd.planets.n_alive - 1);

						for (size_t i = 1; i < hd.planets.n_alive; i++)
						{
							int esign;
							double a, e, in, capom, om, f;
							to_elements(hd.planets.m[i] + hd.planets.m[0], hd.planets.r[i], hd.planets.v[i],
								&esign, &a, &e, &in, &capom, &om, &f);

							int id = static_cast<uint32_t>(hd.planets.id[i]);
							float af = static_cast<float>(a);
							float ef = static_cast<float>(e);
							float if_ = static_cast<float>(in);

							write_binary(trackb, id);
							write_binary(trackb, af);
							write_binary(trackb, ef);
							write_binary(trackb, if_);
							track << id << " " << af << " " << ef << " " << if_ << std::endl;
						}
						track << hd.particles.n_alive << std::endl;
						write_binary(trackb, static_cast<uint32_t>(hd.particles.n_alive));
						for (size_t i = 0; i < hd.particles.n_alive; i++)
						{
							int esign;
							double a, e, in, capom, om, f;
							to_elements(hd.planets.m[0], hd.particles.r[i], hd.particles.v[i],
								&esign, &a, &e, &in, &capom, &om, &f);

							int id = static_cast<uint32_t>(hd.particles.id[i]);
							float af = static_cast<float>(a);
							float ef = static_cast<float>(e);
							float if_ = static_cast<float>(in);
							write_binary(trackb, id);
							write_binary(trackb, af);
							write_binary(trackb, ef);
							write_binary(trackb, if_);
							track << id << " " << af << " " << ef << " " << if_ << std::endl;
						}
					});
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
		backtrace_symbols_fd(array, size, 2);

		tout << "Exception caught: " << std::endl;
		tout << e.what() << std::endl;
		tout << "Recovering data." << std::endl;
		crashed = true;
	}

	ex.finish();
	tout << "Saving to disk. t = " << ex.t << std::endl;
	save_data(hd, config);
	config.t_f = ex.t_f - ex.t_0 + ex.t;
	config.t_0 = ex.t;
	std::ofstream configout(joinpath(config.outfolder, "config.out"));
	write_configuration(configout, config);

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return crashed ? -1 : 0;
}
