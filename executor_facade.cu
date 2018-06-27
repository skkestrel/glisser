#include "executor_facade.h"
#include "util.h"
#include "data.cuh"
#include "executor.cuh"

ExecutorFacade::ExecutorFacade(HostData& hd, std::ostream& out) : hd(hd), output(out)
{
	dd = std::make_unique<DeviceData>();
	impl = std::make_unique<Executor>(hd, *dd.get(), out);

	t_0 = &impl->t_0;
	t = &impl->t;
	dt = &impl->dt;
	t_f = &impl->t_f;
	e_0 = &impl->e_0;
	tbsize = &impl->tbsize;
	ce_n1 = &impl->ce_n1;
	ce_n2 = &impl->ce_n2;
	resolve_encounters = &impl->resolve_encounters;
	encounter_output = &impl->encounter_output;
}

ExecutorFacade::~ExecutorFacade()
{
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
