#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <thread>
#include <iomanip>

#include "../data.h"
#include "../wh.h"
#include "../convert.h"
#include "../util.h"

int main(int argv, char** argc)
{
	Configuration config;
	HostData hd;

	bool ishelio;
	bool binary = false;
	bool momentum = false;

	for (size_t i = 1; i < argv; i++)
	{
		std::string arg(argc[i - 1]);

		if (arg == "momentum")
		{
			momentum = true;
		}
		else if (arg == "velocity")
		{
			momentum = false;
		}
		else if (arg == "binary")
		{
			binary = true;
		}
		else if (arg == "ascii")
		{
			binary = false;
		}
		else if (arg == "read")
		{
			config.hybridin = std::string(argc[i]);
			config.readbinary = binary;
			config.readmomenta = momentum;

			load_data(hd, config);

			if (hd.planets.r[0].lensq() == 0 && hd.planets.r[0].lensq() == 0) ishelio = true;
			else
			{
				to_bary(hd);
			}

			i++;
		}
		else if (arg == "read-split")
		{
			config.plin = std::string(argc[i]);
			config.icsin = std::string(argc[i + 1]);

			config.readmomenta = momentum;
			config.readbinary = false;
			config.readsplit = true;

			load_data(hd, config);

			if (hd.planets.r[0].lensq() == 0 && hd.planets.r[0].lensq() == 0) ishelio = true;
			else
			{
				to_bary(hd);
			}

			i += 2;
		}
		else if (arg == "write")
		{
			config.hybridout = std::string(argc[i]);
			config.writebinary = binary;
			config.writemomenta = momentum;

			hd.planets_snapshot = HostPlanetSnapshot(hd.planets);
			save_data(hd, config, config.hybridout);

			i++;
		}
		else if (arg == "to-bary")
		{
			to_bary(hd);
		}
		else if (arg == "to-helio")
		{
			to_helio(hd);
		}
		else if (arg == "to-elements")
		{
			double a, e, i, capom, om, f;
			int esign;

			double totalmass = 0;
			for (size_t j = 0; j < hd.planets.n; j++)
			{
				totalmass += hd.planets.m[j];
			}

			for (size_t j = 1; j < hd.planets.n; j++)
			{
				if (ishelio)
				{
					to_elements(hd.planets.m[j] + hd.planets.m[0], hd.planets.r[j], hd.planets.v[j], &esign, &a, &e, &i, &capom, &om, &f);
				}
				else
				{
					to_elements(totalmass, hd.planets.r[j], hd.planets.v[j], &esign, &a, &e, &i, &capom, &om, &f);
				}

				if (esign == 0)
				{
					std::cout << "Parabolic orbit detected!" << std::endl;
				}

				hd.planets.r[j].x = a;
				hd.planets.r[j].y = e;
				hd.planets.r[j].z = i;
				hd.planets.v[j].x = capom;
				hd.planets.v[j].y = om;
				hd.planets.v[j].z = f;
			}

			for (size_t j = 0; j < hd.particles.n; j++)
			{
				if (ishelio)
				{
					to_elements(hd.planets.m[0], hd.particles.r[j], hd.particles.v[j], &esign, &a, &e, &i, &capom, &om, &f);
				}
				else
				{
					to_elements(totalmass, hd.particles.r[j], hd.particles.v[j], &esign, &a, &e, &i, &capom, &om, &f);
				}

				if (esign == 0)
				{
					std::cout << "Parabolic orbit detected!" << std::endl;
				}

				hd.particles.r[j].x = a;
				hd.particles.r[j].y = e;
				hd.particles.r[j].z = i;
				hd.particles.v[j].x = capom;
				hd.particles.v[j].y = om;
				hd.particles.v[j].z = f;
			}
		}
		else
		{
			std::cout << "?" << std::endl;
			return -1;
		}
	}
}
