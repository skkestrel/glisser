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

	for (size_t i = 1; i < argv; i++)
	{
		std::string arg(argc[i - 1]);

		if (arg == "read")
		{
			config.hybridin = std::string(argc[i]);
			load_data(hd, config);
			hd.planets_snapshot = HostPlanetSnapshot(hd.planets);

			if (hd.planets.r[0].lensq() == 0 && hd.planets.r[0].lensq() == 0) ishelio = true;
			else
			{
				to_bary(hd);
			}

			i++;
		}

		if (arg == "write")
		{
			config.hybridout = std::string(argc[i]);
			save_data(hd, config, config.hybridout);

			i++;
		}

		if (arg == "read-split")
		{
			config.plin = std::string(argc[i]);
			config.icsin = std::string(argc[i + 1]);
			config.readsplit = true;
			load_data(hd, config);
			hd.planets_snapshot = HostPlanetSnapshot(hd.planets);

			if (hd.planets.r[0].lensq() == 0 && hd.planets.r[0].lensq() == 0) ishelio = true;
			else
			{
				to_bary(hd);
			}

			i += 2;
		}

		if (arg == "to-bary")
		{
			to_bary(hd);
		}

		if (arg == "to-helio")
		{
			to_helio(hd);
		}

		if (arg == "to-elements")
		{
			double a, e, i, capom, om, f;
			int esign;

			double totalmass = 0;
			for (size_t j = 0; j < hd.planets_snapshot.n; j++)
			{
				totalmass += hd.planets_snapshot.m[j];
			}

			for (size_t j = 1; j < hd.planets_snapshot.n; j++)
			{
				if (ishelio)
				{
					to_elements(hd.planets_snapshot.m[j] + hd.planets_snapshot.m[0], hd.planets_snapshot.r[j], hd.planets_snapshot.v[j], &esign, &a, &e, &i, &capom, &om, &f);
				}
				else
				{
					to_elements(totalmass, hd.planets_snapshot.r[j], hd.planets_snapshot.v[j], &esign, &a, &e, &i, &capom, &om, &f);
				}

				if (esign == 0)
				{
					std::cout << "Parabolic orbit detected!" << std::endl;
				}

				hd.planets_snapshot.r[j].x = a;
				hd.planets_snapshot.r[j].y = e;
				hd.planets_snapshot.r[j].z = i;
				hd.planets_snapshot.v[j].x = capom;
				hd.planets_snapshot.v[j].y = om;
				hd.planets_snapshot.v[j].z = f;
			}

			for (size_t j = 0; j < hd.particles.n; j++)
			{
				if (ishelio)
				{
					to_elements(hd.planets_snapshot.m[0], hd.particles.r[j], hd.particles.v[j], &esign, &a, &e, &i, &capom, &om, &f);
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
	}
}
