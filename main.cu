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

int main(int argv, char** argc)
{
	if (argv < 4)
	{
		std::cerr << "Please specify time step and final time " << argc[0] << " <CURTIME> <TIMESTEP> <FINALTIME> [<MAXPARTICLES>]" << std::endl;
		return -1;
	}

	HostData hd;
	DeviceData dd;
	Executor ex(hd, dd, std::cout);

	ex.t = std::stod(argc[1]);
	ex.t_0 = ex.t;
	ex.dt = std::stod(argc[2]);
	ex.t_f = std::stod(argc[3]);
	ex.tbsize = 1024;
	ex.print_every = 10;

	size_t max_particle = 0;
	if (argv >= 5) max_particle = static_cast<size_t>(std::stoi(argc[4]));

	if (load_data(hd, "pl.in", "ics.in", ex.tbsize, max_particle, false)) return -1;

	std::ofstream timelog("time.out");
	std::ofstream disclog("disc.out");

	std::time_t t = std::time(nullptr);
	std::tm tm = *std::localtime(&t);

	ex.timing_output = &timelog;
	ex.discard_output = &disclog;

	timelog << "start " << std::put_time(&tm, "%c %Z") << std::endl;

	ex.init();
	std::cout << "Saving to disk." << std::endl;
	save_data(hd, "pl.part.out", "ics.part.out");

	size_t counter = 0;
	size_t save_every = 50;

	while (ex.t < ex.t_f)
	{
		ex.loop();

		if (counter++ % save_every == 0)
		{
			std::cout << "Saving to disk. t = " << ex.t << " minus one timeblock" << std::endl;
			save_data(hd, "pl.out", "ics.out");
		}
	}
	ex.finish();

	std::cout << "Saving to disk." << std::endl;
	save_data(hd, "pl.out", "ics.out");

	t = std::time(nullptr);
	tm = *std::localtime(&t);
	timelog << "end " << std::put_time(&tm, "%c %Z") << std::endl;

	return 0;
}
