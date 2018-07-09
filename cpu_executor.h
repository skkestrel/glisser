#include "integrator.h"
#include "wh.h"
#include <ctime>
#include <chrono>
#include <functional>
#include <ostream>

struct CPUExecutor
{
	HostData& hd;
	WHIntegrator integrator;

	float64_t t;
	float64_t e_0;

	std::ostream& output;
	std::ostream* encounter_output;

	const Configuration& config;

	std::chrono::time_point<std::chrono::high_resolution_clock> starttime;

	std::vector<std::function<void()>> work;

	CPUExecutor(const CPUExecutor&) = delete;
	CPUExecutor(HostData& hd, const Configuration& config, std::ostream& out);

	void init();

	double time() const;
	void loop(double* cputime);
	void add_job(const std::function<void()>& job);

	void resync();
	void finish();
	void step_planets();
};
