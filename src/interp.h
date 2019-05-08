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

		// relative_t should be t - t0
		void fill(sr::data::HostPlanetPhaseSpace& pl, size_t nstep, double relative_t, double dt);

		// This is the effective dt, calculated by taking the nearest timestep to the user-defined timestep
		// this ensures that itme chunk boundaries always occur on the lookup boundary
		float64_t eff_dt;

		// number of effective timesteps
		size_t n_ts;

		// the current timestep number in the current lookup interval
		size_t cur_ts;

		// t - t0; this is needed because when t0 is very large, error in t - t0 accumulates extremely quickly
		float64_t rel_t;

		Vf64_3 aei0, aei1, aei_m1;
		Vf64_3 oom0, oom1, oom_m1;
		float64_t t0, t1, t_m1;

		private:
		std::ifstream input;
		Vf64_3 daei, doom;

		double user_dt;

		std::vector<double> mmfreq;
		std::unordered_map<uint32_t, size_t> idmap;
	};

	class EOSError : public std::runtime_error
	{
		public:
			inline EOSError() : std::runtime_error("EOS reached for lookup file.") { }
	};
}
}
