#pragma once
#include "data.h"
#include <unistd.h>

namespace sr
{
namespace swift
{
	class SwiftEncounterIntegrator
	{
	public:
		SwiftEncounterIntegrator(const sr::data::Configuration& config, double t, size_t prev_tbsize, size_t cur_tbsize);
		void begin_integrate(const sr::data::HostPlanetPhaseSpace& pl, const sr::data::HostParticlePhaseSpace& pa);
		void end_integrate(sr::data::HostParticlePhaseSpace& pa);

		void write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, double time, std::string dest) const;
		void write_param_in(std::string dest) const;
		void write_tp_in(const sr::data::HostParticlePhaseSpace& pa, size_t chunk_begin, size_t chunk_end, std::string dest) const;
		void write_pl_in(const sr::data::HostPlanetPhaseSpace& pl, std::string dest) const;

		std::vector<std::vector<int>> istat;
		std::vector<std::vector<int>> irstat;
		
		const int NSTAT = 3 + 11 - 1;
		const int NSTATR = NSTAT;

		// TODO this shoudl read from interp file
		void set_planetary_history();
	private:
		double t;
		size_t prev_tbsize;
		size_t cur_tbsize;
		const sr::data::Configuration& config;


		std::vector<pid_t> _children;
	};
}
}
