#include "data.cuh"
#include "data.h"
#include "wh.cuh"
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

	struct ExecutorData
	{
		std::vector<f64_3> r, v;
		std::vector<uint32_t> id, deathtime_index;
		std::vector<uint16_t> deathflags;

		ExecutorData();
		ExecutorData(size_t size);
	};

	struct Executor
	{
		HostData& hd;
		DeviceData& dd;
		sr::wh::WHCudaIntegrator integrator;
		sr::interp::Interpolator interpolator;

		cudaStream_t main_stream, dth_stream, htd_stream, par_stream;
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

		size_t resync_counter;
		bool ending_lookup_interval;

		const Configuration& config;

		std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

		std::vector<std::function<void()>> work;

		Executor(const Executor&) = delete;
		Executor(HostData& hd, DeviceData& dd, const Configuration& config, std::ostream& out);

		void init();
		void upload_data(size_t begin, size_t length);
		void upload_planet_log();
		void download_data(bool ignore_errors = false);


		double time() const;
		void loop(double* cputime, double* gputime);
		void handle_encounters(bool do_work);
		void add_job(const std::function<void()>& job);
		void resync();
		void finish();
		void swap_logs();
		void update_planets();
		void assert_true(bool cond, std::string msg);
	};
}
}
