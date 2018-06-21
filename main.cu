/*************************************************************

.---. .            .         .   .-,--.                        
\___  |- ,-. ,-. ,-| . . ,-. |-   `|__/ ,-. .  , ,-. ,-. . ,-. 
    \ |  ,-| |   | | | | `-. |    /  \  |-' | /  |-' |   | |-' 
`---' `' `-^ '   `-^ `-^ `-' `'   `-  ` `-' `'   `-' '   ' `-' 

*************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#include "types.cuh"
#include "executor.cuh"
#include "data.h"
#include "wh.h"
#include "convert.h"

#include <execinfo.h>

#include <csignal>

Executor* EXECUTOR;

volatile sig_atomic_t end_loop = 0;

template<typename T>
void write_binary(std::ostream& o, const T& t)
{
	o.write(reinterpret_cast<const char*>(&t), sizeof(t));
}

void term(int signum)
{
	(void) signum;
	end_loop = 1;
}

std::string getpath(const std::string& base, const std::string& file)
{
	return base + "/" + file;
}

int main(int argv, char** argc)
{
	std::string configin = "config.in";
	if (argv >= 2) configin = std::string(argc[1]);
	std::cout << "Reading from configuration file " << configin << std::endl;

	std::ifstream configfile(configin);

	signal(SIGTERM, term);
	signal(SIGINT, term);

	Configuration config;
	if (read_configuration(configfile, &config)) return -1;
	write_configuration(std::cout, config);

	HostData hd;
	DeviceData dd;
	Executor ex(hd, dd, std::cout);
	EXECUTOR = &ex;

	ex.t = config.t;
	ex.t_0 = config.t;
	ex.dt = config.dt;
	ex.t_f = config.t_f;
	ex.tbsize = config.tbsize;
	ex.resolve_encounters = config.resolve_encounters;
	ex.ce_factor = config.ce_factor;
	ex.print_every = config.print_every;

	if (load_data(hd, config.plin, config.icsin, config.tbsize, config.ce_factor, config.max_particle, config.readmomenta)) return -1;

	std::ofstream timelog(getpath(config.outfolder, "time.out"));
	std::ofstream disclog(getpath(config.outfolder, "discard.out"));

	std::ofstream periodic, periodicb;
	if (config.enable_ascii_track) periodic = std::ofstream(getpath(config.outfolder, "periodic.out"));
	if (config.enable_binary_track) periodicb = std::ofstream(getpath(config.outfolder, "periodicb.out"), std::ios_base::binary);

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	ex.timing_output = &timelog;
	ex.discard_output = &disclog;

	timelog << "start " << std::put_time(&tm, "%c %Z") << std::endl;

	ex.init();

	size_t counter = 0;

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
						std::cout << "Saving to disk. t = " << ex.t << std::endl;
						save_data(hd, getpath(config.outfolder, "pl.dump.out"), getpath(config.outfolder, "ics.dump.out"));
					});
			}

			if (config.periodic_every != 0 && counter % config.periodic_every == 0)
			{
				ex.download_data();

				ex.add_job([&]()
					{
						periodic << std::setprecision(7);
						periodic << ex.t << std::endl;
						periodic << hd.planets.n_alive - 1 << std::endl;

						write_binary(periodicb, ex.t);
						write_binary(periodicb, hd.planets.n_alive - 1);

						for (size_t i = 1; i < hd.planets.n_alive; i++)
						{
							int esign;
							double a, e, in, capom, om, f;
							to_elements(hd.planets.m[i] + hd.planets.m[0], hd.planets.r[i], hd.planets.v[i],
								&esign, &a, &e, &in, &capom, &om, &f);

							int id = static_cast<int>(hd.planets.id[i]);
							float af = static_cast<float>(a);
							float ef = static_cast<float>(e);
							float if_ = static_cast<float>(in);

							write_binary(periodicb, id);
							write_binary(periodicb, af);
							write_binary(periodicb, ef);
							write_binary(periodicb, if_);
							periodic << id << " " << af << " " << ef << " " << if_ << std::endl;
						}
						periodic << hd.particles.n_alive << std::endl;
						write_binary(periodicb, hd.particles.n_alive);
						for (size_t i = 0; i < hd.particles.n_alive; i++)
						{
							int esign;
							double a, e, in, capom, om, f;
							to_elements(hd.planets.m[0], hd.particles.r[i], hd.particles.v[i],
								&esign, &a, &e, &in, &capom, &om, &f);

							int id = static_cast<int>(hd.particles.id[i]);
							float af = static_cast<float>(a);
							float ef = static_cast<float>(e);
							float if_ = static_cast<float>(in);
							write_binary(periodicb, id);
							write_binary(periodicb, af);
							write_binary(periodicb, ef);
							write_binary(periodicb, if_);
							periodic << id << " " << af << " " << ef << " " << if_ << std::endl;
						}
					});
			}

			if (end_loop)
			{
				std::cout << "Caught signal." << std::endl;
				throw std::exception();
			}
		}
		ex.finish();

		std::cout << "Saving to disk." << std::endl;
		save_data(hd, getpath(config.outfolder, "pl.out"), getpath(config.outfolder, "ics.out"));
	}
	catch (std::exception e)
	{
		void* array[50];
		size_t size = backtrace(array, 50);
		backtrace_symbols_fd(array, size, 2);

		std::cout << "Recovering data." << std::endl;

		ex.finish();
		std::cout << "Saving to disk. t = " << ex.t << std::endl;

		save_data(hd, getpath(config.outfolder, "pl.crash.out"), getpath(config.outfolder, "ics.crash.out"));
	}

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return 0;
}
