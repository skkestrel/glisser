#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#include "../docopt/docopt.h"
#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"
#include "../src/interp.h"

static const char USAGE[] = R"(lookup-info
Usage:
    lookup-info [options] <file>

Options:
    -h, --help                         Show this screen.
)";

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "lookup-info");

	try
	{
		sr::data::Configuration config = sr::data::Configuration::create_dummy();

		sr::data::HostData hd;
		hd.planets = sr::data::HostPlanetPhaseSpace(99);

		sr::interp::Interpolator interp(config, hd.planets, args["<file>"].asString());

		size_t i = 0;
		while (true)
		{
			if (i++ % 10000 == 0) std::cout << "\r" << interp.t1 << "                 ";
			interp.next(hd.planets);
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << e.what() << std::endl;
		return -1;
	}
	std::cout << std::endl;

	return 0;
}
