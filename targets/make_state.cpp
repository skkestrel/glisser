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
Usage:
    make-state [options] <planet-file> <output-file>

Options:
    -h, --help          Show this screen.
    -n <N>              Number of particles to generate [default: 1000]
    -a <range>          a range [default: 0,10]
    -e <range>          e range [default: 0,0.1]
    -i <range>          i range [default: 0,50]
    -O <range>          O range [default: 0,360]
    -o <range>          o range [default: 0,360]
    -M <range>          M range [default: 0,360]
    -G <g>              G [default: 1]
    -B, --barycentric   Generate in barycentric coords
)";

using namespace sr::data;
using namespace sr::convert;
const double EPS = 1e-13;

void split_and_load(std::array<double, 12>& range, std::stringstream& ss, size_t n)
{
	std::string token;
	std::getline(ss, token, ',');
	range[2 * n] = std::stod(token);
	std::getline(ss, token, ',');
	range[2 * n + 1] = std::stod(token);
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

		hd.particles = sr::data::HostParticlePhaseSpace(std::stol(args["-n"].asString()), false);

		std::array<double, 12> range;
		std::stringstream ss(args["-a"].asString());
		split_and_load(range, ss, 0);

		ss = std::stringstream(args["-e"].asString());
		split_and_load(range, ss, 1);

		ss = std::stringstream(args["-i"].asString());
		split_and_load(range, ss, 2);

		ss = std::stringstream(args["-O"].asString());
		split_and_load(range, ss, 3);

		ss = std::stringstream(args["-o"].asString());
		split_and_load(range, ss, 4);

		ss = std::stringstream(args["-M"].asString());
		split_and_load(range, ss, 5);

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

		for (size_t i = 0; i < hd.particles.n(); i++)
		{
			double a = adis(gen);
			double e = edis(gen);
			double inc = idis(gen);
			double O = Odis(gen);
			double o = odis(gen);
			double M = Mdis(gen);

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
			hd.particles.id()[i] = static_cast<uint32_t>(i);
			hd.particles.deathflags()[i] = 0;
			hd.particles.deathtime()[i] = 0;

			// sr::convert::to_elements(mu, hd.particles.r()[i], hd.particles.v()[i], nullptr, &a, &e, &inc, &O, &o, &M);
			// std::cout << a << " " << e << " " << inc << " " << O << " " << o << " " << M << std::endl;
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
