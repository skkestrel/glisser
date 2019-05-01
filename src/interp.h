#pragma once
#include "types.h"
#include "util.h"
#include "data.h"
#include "convert.h"
#include <unordered_map>
#include <fstream>

namespace sr
{
namespace interp
{
	class Interpolator
	{
		public:
		Interpolator();
		Interpolator(const sr::data::Configuration& config, sr::data::HostPlanetPhaseSpace& pl, std::string file);
		void next(sr::data::HostPlanetPhaseSpace& pl);
		void fill(sr::data::HostPlanetPhaseSpace& pl, size_t nstep, double t, double dt);

		private:
		std::ifstream input;
		Vf64_3 aei0, aei1;
		Vf64_3 oom0, oom1;
		Vf64_3 daei, doom;
		float64_t t0, t1;
		std::vector<double> mmfreq;
		std::unordered_map<uint32_t, size_t> idmap;

		bool resolve_encounters;
		size_t fast_factor;
	};

	class EOSError
	{
	};
}
}
