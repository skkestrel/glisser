#include "executor_facade.h"
#include "executor.cuh"

ExecutorFacade::ExecutorFacade(HostData& hd, std::ostream& out) : hd(hd), output(out)
{
	dd = new DeviceData();
	impl = new Executor(hd, *dd, out);

	wh_alloc = &impl->wh_alloc;
	t_0 = &impl->t_0;
	t = &impl->t;
	dt = &impl->dt;
	t_f = &impl->t_f;
	e_0 = &impl->e_0;
	tbsize = &impl->tbsize;
	ce_factor = &impl->ce_factor;
	resolve_encounters = &impl->resolve_encounters;
	encounter_output = &impl->encounter_output;
}

ExecutorFacade::~ExecutorFacade()
{
	delete impl;
}

void ExecutorFacade::init()
{
	impl->init();
}

void ExecutorFacade::download_data()
{
	impl->download_data();
}

double ExecutorFacade::time() const
{
	return impl->time();
}

void ExecutorFacade::loop()
{
	impl->loop();
}

void ExecutorFacade::finish()
{
	impl->finish();
}

void ExecutorFacade::add_job(const std::function<void()>& job)
{
	impl->add_job(job);
}
