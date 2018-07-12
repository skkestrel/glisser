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
		("I,initial", "Initial state input file", cxxopts::value<std::string>())
		("u,union", "Take union of criteria instead of intersection")
		("b,binary", "Read binary input")
		("B,barycentric", "Calculate barycentric elements")
		("criteria", "Criteria (eg. a>25)", cxxopts::value<std::vector<std::string>>());

	options.parse_positional({ "input", "criteria" });

	try
	{
		auto result = options.parse(argc, argv);

		bool use_union = result.count("u") > 0;

		if (result.count("i") == 0)
		{
			throw cxxopts::OptionException("Required option -i");
		}

		std::vector<Criterion> criteria;

		if (result.count("criteria") > 0)
		{
			for (const std::string& str : result["criteria"].as<std::vector<std::string>>())
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
					std::ostringstream ss;
					ss << "Could not parse criterion " << str;
					throw cxxopts::OptionException(ss.str());
				}

				crit.variable = str.substr(0, split);
				crit.konst = std::stod(str.substr(split + 1, std::string::npos));

				criteria.push_back(crit);
			}	
		}

		sr::data::Configuration config = sr::data::Configuration::create_dummy();
		sr::data::HostData hd;
		sr::data::HostData hd_init;

		config.hybridin = result["i"].as<std::string>();
		config.readbinary = result.count("b") > 0;
		config.readmomenta = false;

		load_data(hd, config);
		hd.particles.sort_by_id(0, hd.particles.n);
		if (result.count("B") > 0)
		{
			sr::convert::to_bary(hd);
		}
		else
		{
			sr::convert::to_helio(hd);
		}


		bool has_init = result.count("I") > 0;
		if (has_init)
		{
			config.hybridin = result["I"].as<std::string>();
			load_data(hd_init, config);

			hd_init.particles.sort_by_id(0, hd_init.particles.n);

			if (hd_init.particles.id != hd.particles.id)
			{
				throw cxxopts::OptionException("Initial state is not congruent with input state");
			}

			if (result.count("B") > 0)
			{
				sr::convert::to_bary(hd_init);
			}
			else
			{
				sr::convert::to_helio(hd_init);
			}
		}

		std::vector<size_t> candidates;
		for (size_t i = 0; i < hd.particles.n; i++)
		{
			int esign, esign_I;
			double A, E, I, CAPOM, OM, F;
			double A_I, E_I, I_I, CAPOM_I, OM_I, F_I;

			sr::convert::to_elements(hd.planets.m[0], hd.particles.r[i], hd.particles.v[i], &esign, &A, &E, &I, &CAPOM, &OM, &F);

			if (has_init)
			{
				sr::convert::to_elements(hd_init.planets.m[0], hd_init.particles.r[i], hd_init.particles.v[i], &esign_I, &A_I, &E_I, &I_I, &CAPOM_I, &OM_I, &F_I);
			}

			bool ok = !use_union;

			for (const auto& crit : criteria)
			{
				double val;
				if (crit.variable == "a") val = A;
				else if (crit.variable == "e") val = E;
				else if (crit.variable == "i") val = I;
				else if (crit.variable == "O") val = CAPOM;
				else if (crit.variable == "o") val = OM;
				else if (crit.variable == "f") val = F;
				else if (crit.variable == "id") val = static_cast<double>(hd.particles.id[i]);
				else if (crit.variable == "a_i") val = A_I;
				else if (crit.variable == "e_i") val = E_I;
				else if (crit.variable == "i_i") val = I_I;
				else if (crit.variable == "O_i") val = CAPOM_I;
				else if (crit.variable == "o_i") val = OM_I;
				else if (crit.variable == "f_i") val = F_I;
				else if (crit.variable == "deathtime") val = hd.particles.deathtime[i];
				else if (crit.variable == "killer")
				{
					if (hd.particles.deathflags[i] == 0) val = -1;
					else val = static_cast<double>(hd.particles.deathflags[i] >> 8);
				}
				else
				{
					std::ostringstream ss;
					ss << "Unknown value " << crit.variable;
					throw cxxopts::OptionException(ss.str());
				}

				bool newbool;
				if (crit.comparison == 1)
				{
					newbool = (val > crit.konst);
				}
				else if (crit.comparison == -1)
				{
					newbool = (val < crit.konst);
				}
				else if (crit.comparison == 0)
				{
					newbool = std::abs(val - crit.konst) < EPS;
				}

				if (use_union)
				{
					ok = ok || newbool;
					if (ok) break;
				}
				else
				{
					ok = ok && newbool;
					if (!ok) break;
				}
			}

			if (ok)
			{
				candidates.push_back(i);
			}
		}

		std::cout << "ID      |   a        e        i        Om       om       f" << std::endl;
		std::cout << std::setfill(' ');

		int num = 0;
		for (auto& i : candidates)
		{
			int esign;
			double A, E, I, CAPOM, OM, F;
			sr::convert::to_elements(hd.planets.m[0], hd.particles.r[i], hd.particles.v[i], &esign, &A, &E, &I, &CAPOM, &OM, &F);

			std::cout << std::setw(7) << hd.particles.id[i] << " | ";
			std::cout << std::setfill(' ') << std::setw(9) << std::fixed << std::setprecision(4)
				<< std::setw(9) << A << std::setw(9) << E << std::setw(9) << I
				<< std::setw(9) << CAPOM << std::setw(9) << OM << std::setw(9) << F << std::endl;
			num++;

			if (num % 10 == 0 && num != 0)
			{
				std::cout << candidates.size() - num << " left (" << candidates.size() << " total) (press return or q) ";
				int ch = '\0';
				while (ch != '\n' && ch != 'q')
				{
					ch = std::getchar();
				}

				if (ch == 'q') return 0;
			}
		}
	}
	catch (cxxopts::OptionException& e)
	{
		std::cout << e.what() << std::endl;
		std::cout << options.help() << std::endl;
		return -1;
	}

	return 0;
}
