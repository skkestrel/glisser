#pragma once
#include "data.h"
#include "interp.h"
#include <unistd.h>
#include <unordered_map>

namespace sr
{
namespace swift
{
	struct ChildProcess
	{
		::pid_t pid;
		int piper;
		std::string tpout;
		size_t chunk_begin, chunk_end;

		inline ChildProcess(::pid_t _pid, int r, std::string _tpout, size_t _chunk_begin, size_t _chunk_end) : pid(_pid), piper(r), tpout(_tpout), chunk_begin(_chunk_begin), chunk_end(_chunk_end) { }
	};

	class SwiftEncounterIntegrator
	{
	public:
		SwiftEncounterIntegrator();
		SwiftEncounterIntegrator(const sr::data::Configuration& config, size_t npa);

		void begin_integrate(
				const sr::data::HostPlanetPhaseSpace& pl,
				const sr::data::HostParticlePhaseSpace& pa,
				const sr::interp::Interpolator& interp,
				bool old,
				double t,
				double rel_t,
				double dt,
				size_t prev_tbsize,
				size_t cur_tbsize
		);

		void end_integrate(sr::data::HostParticlePhaseSpace& pa);

		void write_planetary_history(const sr::data::HostPlanetPhaseSpace& pl, const sr::interp::Interpolator& interp, std::string dest);
		void write_param_in(std::string dest) const;
		void write_tp_in(const sr::data::HostParticlePhaseSpace& pa, size_t chunk_begin, size_t chunk_end, std::string dest) const;
		void write_pl_in(const sr::data::HostPlanetPhaseSpace& pl, std::string dest) const;
		void write_stat(std::string dest) const;

		// maps particle id to index in istat and rstat
		std::unordered_map<uint32_t, size_t> statmap;

		std::vector<std::vector<int32_t>> istat;
		std::vector<std::vector<double>> rstat;
	private:
		bool old_log;
		double t, rel_t, dt;
		size_t prev_tbsize;
		size_t cur_tbsize;

		size_t swift_statlen;
		std::string swift_path;
		std::string outfolder;
		size_t swift_part_min;
		size_t num_swift;
		double outer_radius;

		Vu32 planet_id_list;

		std::vector<ChildProcess> _children;
	};
}
}
