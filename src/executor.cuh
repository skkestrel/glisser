#include "data.cuh"
#include "data.h"
#include "wh.cuh"
#include "swift.h"
#include "interp.h"
#include <ctime>
#include <chrono>
#include <functional>
#include <ostream>

namespace sr
{
namespace exec
{
	using namespace sr::wh;
	using namespace sr::util;
	using namespace sr::data;

	struct Executor
	{
		HostData& hd;
		DeviceData& dd;
		sr::wh::WHCudaIntegrator integrator;
		sr::swift::SwiftEncounterIntegrator swift;
		sr::interp::Interpolator interpolator;

		cudaStream_t main_stream, dth_stream, htd_stream, sort_stream;
		cudaEvent_t start_event, cpu_finish_event, gpu_finish_event;

		float64_t t;

		float64_t e_0;

		DeviceParticlePhaseSpace rollback_state;

		size_t prev_tbsize;
		size_t cur_tbsize;
		
		float64_t prev_dt;
		float64_t cur_dt;

		std::ostream& output;
		std::ostream* encounter_output;
		std::ofstream discard_output;
		std::ofstream temp_log;

		size_t resync_counter;
		bool starting_lookup_interval;

		const Configuration& config;

		std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

		std::vector<std::function<void()>> work;

		Executor(const Executor&) = delete;
		Executor(HostData& hd, DeviceData& dd, const Configuration& config, std::ostream& out);

		void init();
		void write_encounter(size_t begin, size_t end, double prev_t);
		void upload_data(size_t begin, size_t length);
		void upload_planet_log();
		void download_data(size_t begin, size_t length);

		double time() const;
		bool loop(double* cputime, double* gputime, double* totaltime, size_t* nencounter);
		size_t handle_encounters(bool do_work);
		void add_job(const std::function<void()>& job);
		void resync();
		size_t resync2();
		void finish();
		void swap_logs();
		void update_planets();

		void alloc_packed();

		double t_initswift, t_backup, t_enc, t_writeswift, t_swift, t_readswift, t_delayswift, t_io, t_planet, t_planetup, t_encup, t_sort, t_dl, t_rollback, t_enc2, t_resync;

		int8_t* cpu_packed_mem;
		int8_t* gpu_packed_mem;

		size_t packed_mem_size;

		size_t presort_index;
	};
}
}
