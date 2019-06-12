#include <sstream>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#include "../src/data.h"
#include "../src/wh.h"
#include "../src/convert.h"
#include "../src/util.h"
#include "../docopt/docopt.h"



static const char USAGE[] = R"(make-state

generates a state file from a planet file
angles all in degrees
example: make-state pl.in -a1,10 -e0.1,0.2

Usage:
    make-state [options] <planet-file> <output-file>

Options:
    -h, --help          Show this screen.
    -n <N>              Number of particles to generate [default: 1000]
    -a <range>          a range
    -e <range>          e range
    -q <range>          q range
    -i <range>          i range
    -O <range>          O range
    -o <range>          o range
    -M <range>          M range
    -G <g>              G [default: 1]
    -B, --barycentric   Generate in barycentric coords
)";

using namespace sr::data;
using namespace sr::convert;
const double EPS = 1e-13;

void split_and_load(std::array<double, 12>& range, std::string s, size_t n)
{
	if (std::count(s.begin(), s.end(), ',') == 0)
	{
		range[2 * n + 1] = range[2 * n] = std::stod(s);
	}
	else
	{
		auto ss = std::stringstream(s);

		std::string token;
		std::getline(ss, token, ',');
		range[2 * n] = std::stod(token);
		std::getline(ss, token, ',');
		range[2 * n + 1] = std::stod(token);
	}
}

int main(int argc, char** argv)
{
	std::map<std::string, docopt::value> args = docopt::docopt(USAGE, { argv + 1, argv + argc }, true, "convert-state");

	try
	{
		Configuration config;
		HostData hd;

		std::ifstream inplanet(args["<planet-file>"].asString());
		if (sr::data::load_planet_data(hd.planets, config, inplanet))
		{
			std::cerr << "Could not load planet data" << std::endl;
			return -1;
		}

		double mu = hd.planets.m()[0] * std::stod(args["-G"].asString());
		
		bool gen_bary = static_cast<bool>(args["--barycentric"]);
		if (gen_bary)
		{
			sr::convert::to_bary(hd);
		}

		const double min = 1e-6;

		hd.particles = sr::data::HostParticlePhaseSpace(std::stol(args["-n"].asString()));

		std::array<double, 12> range({ 1, 1, 0.1, 0.1, 0.00001, 0.00001, 0, 360, 0, 360, 0, 360 });

		bool use_q = false;

		if (args["-a"]) split_and_load(range, args["-a"].asString(), 0);
		if (args["-e"]) split_and_load(range, args["-e"].asString(), 1);
		if (args["-i"]) split_and_load(range, args["-i"].asString(), 2);
		if (args["-O"]) split_and_load(range, args["-O"].asString(), 3);
		if (args["-o"]) split_and_load(range, args["-o"].asString(), 4);
		if (args["-M"]) split_and_load(range, args["-M"].asString(), 5);

		if (args["-q"])
		{
			use_q = true;
			split_and_load(range, args["-q"].asString(), 1);
		}

		for (size_t i = 4; i < 12; i++)
		{
			range[i] = range[i] / 360. * 2 * M_PI;
		}

		if (range[0] < min) range[0] = min;
		if (range[1] < min) range[1] = min;
		if (range[2] < min) range[2] = min;

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> adis(range[0], range[1]);
		std::uniform_real_distribution<> edis(range[2], range[3]);
		std::uniform_real_distribution<> idis(range[4], range[5]);
		std::uniform_real_distribution<> Odis(range[6], range[7]);
		std::uniform_real_distribution<> odis(range[8], range[9]);
		std::uniform_real_distribution<> Mdis(range[10], range[11]);

		for (uint32_t i = 0; i < hd.particles.n(); i++)
		{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
			double a = range[0] == range[1] ? range[0] : adis(gen);
			double e = range[2] == range[3] ? range[2] : edis(gen);

			// q = a(1-e)
			// e = 1 - q/a
			if (use_q) e = 1 - e / a;

			double inc = range[4] == range[5] ? range[4] : idis(gen);
			double O = range[6] == range[7] ? range[6] : Odis(gen);
			double o = range[8] == range[9] ? range[8] : odis(gen);
			double M = range[10] == range[11] ? range[10] : Mdis(gen);
#pragma GCC diagnostic pop

			// std::cout << a << " " << e << " " << inc << " " << O << " " << o << " " << M << std::endl;


			double sindE, cosdE;
			double ecosE = e;
			double esinE = 0;
			double dE = M + ecosE * std::sin(M);  /* input guess */

			uint32_t it;
			if (sr::wh::kepeq(M,esinE,ecosE,&dE,&sindE,&cosdE, &it)) throw std::runtime_error("?");

			double cosanom = (cosdE - e)/(1.0 - e*cosdE);
			double sinanom = std::sqrt(1.0 - e*e) * sindE/(1.0 - e*cosdE);
			double anom = std::atan2(sinanom,cosanom);
		   
			sr::convert::from_elements(mu,a,e,inc,O,o, anom, &hd.particles.r()[i], &hd.particles.v()[i]);
			hd.particles.id()[i] = i;
			hd.particles.deathflags()[i] = 0;
			hd.particles.deathtime()[i] = 0;

			sr::convert::to_elements(mu, hd.particles.r()[i], hd.particles.v()[i], nullptr, &a, &e, &inc, &O, &o, &M);
			std::cout << a << " " << e << " " << inc << " " << O << " " << o << " " << M << std::endl;
		}

		if (gen_bary)
		{
			sr::convert::to_helio(hd);
		}

		config.hybridout = args["<output-file>"].asString();
		sr::data::save_data(hd.planets.base, hd.particles, config, config.hybridout);
	}
	catch (std::runtime_error& e)
	{
		std::cout << "Error occurred" << std::endl;
		std::cout << e.what() << std::endl;
		return -1;
	}

	return 0;
}
