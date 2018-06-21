#include "data.cuh"
#include "data.h"
#include <ctime>
#include <chrono>
#include <iostream> 
#include <functional>

struct Executor
{
	HostData& hd;
	DeviceData& dd;

	cudaStream_t main_stream, dth_stream, htd_stream, par_stream;

	float64_t t_0, t, dt, t_f;
	float64_t e_0;

	size_t print_every, print_counter;
	size_t tbsize, ce_factor;

	bool resolve_encounters;

	std::ostream& output;
	std::ostream* timing_output;
	std::ostream* discard_output;

	std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

	std::vector<std::function<void()>> work;

	Executor(HostData& hd, DeviceData& dd, std::ostream& out);

	void init();
	void upload_data();
	void upload_planet_log();
	void download_data();

	double time() const;
	void loop();
	void add_job(std::function<void()> job);
	void resync();
	void animate();
	void finish();
	void step_and_upload_planets();
};
