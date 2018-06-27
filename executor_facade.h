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
	std::ostream& output;

	float64_t* t_0, *t, *dt, *t_f;
	float64_t* e_0;

	size_t* tbsize, *ce_n1, *ce_n2;
	bool* resolve_encounters;
	std::ostream** encounter_output;

	std::unique_ptr<Executor> impl;
	std::unique_ptr<DeviceData> dd;

	ExecutorFacade(HostData& hd, std::ostream& out);
	~ExecutorFacade();

	void init();
	void download_data();
	double time() const;
	void loop();
	void add_job(const std::function<void()>& job);
	void finish();
};
