#include "data.cuh"
#include "data.h"
#include <ctime>
#include <chrono>
#include <iostream> 

struct Executor
{
	HostData& hd;
	DeviceData& dd;

	cudaStream_t main_stream, dth_stream, htd_stream, par_stream;

	float64_t t_0, t, dt, t_f;
	float64_t e_0;

	size_t print_every, print_counter;
	size_t tbsize;

	std::ostream& output;
	std::ostream* timing_output;
	std::ostream* discard_output;

	std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

	Executor(HostData& hd, DeviceData& dd, std::ostream& out);

	void init();
	void upload_data();
	void upload_planet_log();
	void download_data();

	double time() const;
	void loop();
	void resync();
	void animate();
	void finish();
};
