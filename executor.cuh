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

	ExecutorData();
	ExecutorData(size_t size);
};

struct Executor
{
	HostData& hd;
	DeviceData& dd;
	WHCudaIntegrator integrator;

	ExecutorData ed;

	cudaStream_t main_stream, dth_stream, htd_stream, par_stream;

	float64_t t_0, t, dt, t_f;
	float64_t e_0;

	size_t tbsize, ce_n1, ce_n2;

	bool resolve_encounters;

	std::ostream& output;
	std::ostream* encounter_output;

	std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

	std::vector<std::function<void()>> work;

	Executor(const Executor&) = delete;
	Executor(HostData& hd, DeviceData& dd, std::ostream& out);

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
