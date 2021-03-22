#include "executor_facade.h"
#include "util.h"
#include "data.cuh"
#include "executor.cuh"

namespace sr
{
namespace exec
{
	using namespace sr::data;

	ExecutorFacade::ExecutorFacade(HostData& _hd, const Configuration& config, std::ostream& out) :
		hd(_hd),
		dd(std::make_unique<DeviceData>()),
		impl(std::make_unique<Executor>(_hd, *dd.get(), config, out)),
		t(impl->t),
		e_0(impl->e_0),
		encounter_output(impl->encounter_output)
	{
	}

	ExecutorFacade::~ExecutorFacade()
	{
	}

	void ExecutorFacade::init()
	{
		impl->init();
	}

	void ExecutorFacade::download_data(size_t begin, size_t length)
	{
		impl->download_data(begin, length);
	}

	double ExecutorFacade::time() const
	{
		return impl->time();
	}

	bool ExecutorFacade::loop(double* cputimeout, double* gputimeout, double* totaltimeout, size_t* nencounter)
	{
		return impl->loop(cputimeout, gputimeout, totaltimeout, nencounter);
	}

	void ExecutorFacade::finish()
	{
		impl->finish();
	}

	void ExecutorFacade::add_job(const std::function<void()>& job)
	{
		impl->add_job(job);
	}
}
}
