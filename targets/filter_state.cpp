#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include "../cxxopts.h"
#pragma GCC diagnostic pop

#include "../data.h"
#include "../wh.h"
#include "../convert.h"
#include "../util.h"

struct Criterion
{
	std::string variable;
	int comparison; // -1: var lt konst, 0: var eq konst, 1: var gt const
	double konst;
};

const double EPS = 1e-13;
int main(int argc, char** argv)
{
	cxxopts::Options options("filter-state", "Find particles in a state file that match criteria");
	options.add_options()
		("i,input", "Input file", cxxopts::value<std::string>())
		("b,binary", "Read binary input")
		("B,barycentric", "Calculate barycentric elements")
		("positional", "Criteria", cxxopts::value<std::vector<std::string>>());

	options.parse_positional({ "input", "positional" });
	auto result = options.parse(argc, argv);

	sr::data::Configuration config;
	sr::data::HostData hd;

	config.hybridin = result["i"].as<std::string>();
	config.readbinary = result.count("b") > 0;
	config.readmomenta = false;

	load_data(hd, config);

	if (result.count("B") > 0)
	{
		sr::convert::to_bary(hd);
	}
	else
	{
		sr::convert::to_helio(hd);
	}

	std::vector<Criterion> criteria;

	for (const std::string& str : result["positional"].as<std::vector<std::string>>())
	{
		Criterion crit;
		size_t foundgt = str.find('>');
		size_t foundeq = str.find('=');
		size_t foundlt = str.find('<');

		size_t split;

		if (foundgt != std::string::npos)
		{
			crit.comparison = 1;
			split = foundgt;
		}
		else if (foundeq != std::string::npos)
		{
			crit.comparison = 0;
			split = foundeq;
		}
		else if (foundlt != std::string::npos)
		{
			crit.comparison = -1;
			split = foundlt;
		}
		else
		{
			std::cout << "?" << std::endl;
			return -1;
		}

		crit.variable = str.substr(0, split);
		crit.konst = std::stod(str.substr(split + 1, std::string::npos));

		criteria.push_back(crit);
	}	


	std::vector<size_t> candidates;
	for (size_t i = 0; i < hd.particles.n; i++)
	{
		int esign;
		double A, E, I, CAPOM, OM, F;
		sr::convert::to_elements(hd.planets.m[0], hd.particles.r[i], hd.particles.v[i], &esign, &A, &E, &I, &CAPOM, &OM, &F);

		bool ok = true;
		for (const auto& crit : criteria)
		{
			double val;
			if (crit.variable == "a") val = A;
			else if (crit.variable == "e") val = E;
			else if (crit.variable == "i") val = I;
			else if (crit.variable == "O") val = CAPOM;
			else if (crit.variable == "o") val = OM;
			else if (crit.variable == "f") val = F;
			else if (crit.variable == "deathtime") val = hd.particles.deathtime[i];
			else if (crit.variable == "killer")
			{
				if (hd.particles.deathflags[i] == 0) val = -1;
				else val = static_cast<double>(hd.particles.deathflags[i] >> 8);
			}
			else
			{
				std::cout << "Unknown value " << crit.variable << std::endl;
				return -1;
			}

			if (crit.comparison == 1)
			{
				ok = ok && (val > crit.konst);
			}
			else if (crit.comparison == -1)
			{
				ok = ok && (val < crit.konst);
			}
			else if (crit.comparison == 0)
			{
				ok = ok && std::abs(val - crit.konst) < EPS;
			}

			if (!ok) break;
		}
		
		if (ok)
		{
			candidates.push_back(i);
		}
	}

	std::cout << "ID      |   A       E       I       Om      om      F" << std::endl;
	std::cout << std::setfill(' ');

	int num = 0;
	for (auto& i : candidates)
	{
		int esign;
		double A, E, I, CAPOM, OM, F;
		sr::convert::to_elements(hd.planets.m[0], hd.particles.r[i], hd.particles.v[i], &esign, &A, &E, &I, &CAPOM, &OM, &F);

		std::cout << std::setw(7) << hd.particles.id[i] << " | ";
		std::cout << std::setfill(' ') << std::setw(8) << std::fixed << std::setprecision(3)
			<< std::setw(8) << A << std::setw(8) << E << std::setw(8) << I
			<< std::setw(8) << CAPOM << std::setw(8) << OM << std::setw(8) << F << std::endl;
		num++;

		if (num % 10 == 0 && num != 0)
		{
			std::cout << candidates.size() - num << " left (" << candidates.size() << " total) (press enter)";
			std::string s;
			std::getline(std::cin, s);
		}
	}
}
