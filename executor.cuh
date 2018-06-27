#include "data.cuh"
#include "data.h"
#include "integrator.cuh"
#include "wh.cuh"
#include <ctime>
#include <chrono>
#include <functional>
#include <ostream>

struct ExecutorData
{
	std::vector<f64_3> r, v;
	std::vector<uint32_t> id, deathtime_index;
	std::vector<uint16_t> deathflags;

	std::vector<uint8_t> encounter_planet_id;

	std::unique_ptr<std::vector<size_t>> gather_indices;

	ExecutorData();
	ExecutorData(size_t size);
};

struct Executor
{
	const Configuration& config;

	HostData& hd;
	DeviceData& dd;
	WHCudaIntegrator integrator;

	ExecutorData ed;

	cudaStream_t main_stream, dth_stream, htd_stream, par_stream;
	cudaEvent_t start_event, cpu_finish_event, gpu_finish_event;

	float64_t t;
	float64_t e_0;

	std::ostream& output;
	std::ostream* encounter_output;

	std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

	std::vector<std::function<void()>> work;

	Executor(const Executor&) = delete;
	Executor(HostData& hd, DeviceData& dd, const Configuration& config, std::ostream& out);

	void init();
	void upload_data();
	void upload_planet_log();
	void download_data();

	double time() const;
	void loop();
	void add_job(const std::function<void()>& job);
	void resync();
	void finish();
	void step_and_upload_planets();
};
