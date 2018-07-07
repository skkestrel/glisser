#pragma once
#include "data.h"
#include "wh.h"
#include <functional>
#include <memory>

struct Executor;
struct DeviceData;

struct ExecutorFacade
{
	HostData& hd;
	std::unique_ptr<DeviceData> dd;
	std::unique_ptr<Executor> impl;

	float64_t& t;
	float64_t& e_0;

	std::ostream*& encounter_output;


	ExecutorFacade(HostData& hd, const Configuration& config, std::ostream& out);
	~ExecutorFacade();

	void init();
	void download_data(bool ignore_errors = false);
	double time() const;
	void loop(double* cputimeout, double* gputimeout);
	void add_job(const std::function<void()>& job);
	void finish();
};
