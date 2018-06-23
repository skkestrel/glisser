#pragma once
#include "data.h"
#include "wh.h"
#include <functional>

struct Executor;
struct DeviceData;

struct ExecutorFacade
{
	HostData& hd;
	std::ostream& output;

	float64_t* t_0, *t, *dt, *t_f;
	float64_t* e_0;

	size_t* tbsize, *ce_factor;
	bool* resolve_encounters;
	WHIntermediate* wh_alloc;
	std::ostream** encounter_output;


	Executor* impl;
	DeviceData* dd;

	ExecutorFacade(HostData& hd, std::ostream& out);
	~ExecutorFacade();

	void init();
	void download_data();
	double time() const;
	void loop();
	void add_job(const std::function<void()>& job);
	void finish();
};
