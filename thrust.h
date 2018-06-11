#ifdef GCC
	#include <vector>
	namespace thrust
	{
		template<typename T>
		using host_vector = std::vector<T>;
		template<typename T>
		using device_vector = std::vector<T>;
	}

	#define __host__
	#define __device__
	#define cudaStreamCreate(...)
	#define cudaStreamSynchronize(...)

	using cudaStream_t = int;
#else
	#include <thrust/system/cuda/execution_policy.h>
	#include <thrust/for_each.h>
	#include <thrust/iterator/zip_iterator.h>
	#include <thrust/host_vector.h>
	#include <thrust/device_vector.h>
#endif

