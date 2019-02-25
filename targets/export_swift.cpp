#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"
#include "../docopt/docopt.h"

static const char USAGE[] = R"(export-swift
Usage:
    export-swift [options] <config> <output-folder>

Options:
    -h, --help         Show this screen.
)";

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "export-swift");

	std::string configin = args["<config>"].asString();
	
	std::cout << "Reading from configuration file " << configin << std::endl;
	
	std::ifstream configfile(configin);

	sr::data::Configuration config;

	try
	{
		read_configuration(configfile, &config);
	}
	catch (std::exception& e)
	{
		std::cerr << "Could not read configuration." << std::endl;
		std::cerr << e.what() << std::endl;
		return -1;
	}

	std::string out_folder = args["<output-folder>"].asString();
	sr::util::make_dir(out_folder);

	if (!sr::util::is_dir_empty(out_folder))
	{
		std::cout << "Output folder is not empty! Do you want to continue?" << std::endl;
		std::cout << "Type \"Yes\" exactly as shown to continue: ";
	
		std::string s;
		std::getline(std::cin, s);

		if (s != "Yes") return -1;
	}

	sr::data::HostData hd;

	if (load_data(hd.planets, hd.particles, config)) return -1;

	std::ofstream plin(sr::util::joinpath(out_folder, "pl.in"));
	std::ofstream tpin(sr::util::joinpath(out_folder, "tp.in"));
	std::ofstream paramin(sr::util::joinpath(out_folder, "param.in"));

	paramin << config.t_0 << " " << config.t_f << " " << config.dt << std::endl;
	paramin << "999999 9999999" << std::endl;
	paramin << "F T F F T F" << std::endl;
	paramin << "0.5 500 200 -1 T" << std::endl;
	paramin << "/dev/null" << std::endl;
	paramin << "unknown" << std::endl;

	plin << hd.planets.n_alive() << std::endl;
	for (size_t i = 0; i < hd.planets.n_alive(); i++)
	{
		plin << hd.planets.m()[i] << std::endl;

		if (i == 0 && (hd.planets.r()[i].lensq() != 0 || hd.planets.v()[i].lensq() != 0))
		{
			throw std::string("planet 0 must be the sun and have 0 vectors");
		}

		plin << hd.planets.r()[i] << std::endl;
		plin << hd.planets.v()[i] << std::endl;
	}

	tpin << hd.particles.n_alive() << std::endl;
	for (size_t i = 0; i < hd.particles.n_alive(); i++)
	{
		tpin << hd.particles.r()[i] << std::endl;
		tpin << hd.particles.v()[i] << std::endl;
		tpin << "0 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
		tpin << "0 0 0 0 0" << std::endl;
		tpin << "0 0 0 0 0" << std::endl;
		tpin << "0 0 0" << std::endl;
	}
}
