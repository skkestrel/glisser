#pragma once
#include "data.h"
#include "interp.h"
#include <unistd.h>

namespace sr
{
namespace swift
{
	struct ChildProcess
	{
		::pid_t pid;
		int piper;
		std::string tpout;

		inline ChildProcess(::pid_t _pid, int r, std::string _tpout) : pid(_pid), piper(r), tpout(_tpout) { }
	};

	class SwiftEncounterIntegrator
	{
	public:
		SwiftEncounterIntegrator(const sr::data::Configuration& config, double t, double dt, size_t prev_tbsize, size_t cur_tbsize);

		void begin_integrate(const sr::data::HostPlanetPhaseSpace& pl, const sr::data::HostParticlePhaseSpace& pa, const sr::interp::Interpolator& interp, bool old);
		void end_integrate(sr::data::HostParticlePhaseSpace& pa);

		void write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, std::string dest) const;
		void write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, const sr::interp::Interpolator& interp, std::string dest, bool old) const;
		void write_param_in(std::string dest) const;
		void write_tp_in(const sr::data::HostParticlePhaseSpace& pa, size_t chunk_begin, size_t chunk_end, std::string dest) const;
		void write_pl_in(const sr::data::HostPlanetPhaseSpace& pl, std::string dest) const;

		std::vector<std::vector<int>> istat;
		std::vector<std::vector<double>> rstat;
	private:
		double t, dt;
		size_t prev_tbsize;
		size_t cur_tbsize;
		const sr::data::Configuration& config;


		std::vector<ChildProcess> _children;
	};
}
}
