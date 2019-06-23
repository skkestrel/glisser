#pragma once
#include <vector>
#include <cstdint>
#include "thrust.cuh"
#include "types.cuh"

template<typename T, typename Alloc>
inline cudaError_t memcpy_dth(std::vector<T, Alloc>& dest, const thrust::device_vector<T>& src, cudaStream_t stream, size_t destbegin = 0, size_t srcbegin = 0, size_t len = static_cast<uint32_t>(-1))
{
	if (len == static_cast<uint32_t>(-1))
	{
		len = src.size();
	}
	if (dest.size() < destbegin + len)
	{
		throw std::exception();
	}

	return cudaMemcpyAsync(dest.data() + destbegin, src.data().get() + srcbegin, len * sizeof(T), cudaMemcpyDeviceToHost, stream);
}

template<typename T>
inline cudaError_t memcpy_dtd(thrust::device_vector<T>& dest, const thrust::device_vector<T>& src, cudaStream_t stream, size_t destbegin = 0, size_t srcbegin = 0, size_t len = static_cast<uint32_t>(-1))
{
	if (len == static_cast<uint32_t>(-1))
	{
		len = src.size();
	}
	if (dest.size() < destbegin + len)
	{
		throw std::exception();
	}

	return cudaMemcpyAsync(dest.data().get() + destbegin, src.data().get() + srcbegin, len * sizeof(T), cudaMemcpyDeviceToDevice, stream);
}

template<typename T, typename Alloc>
inline cudaError_t memcpy_htd(thrust::device_vector<T>& dest, const std::vector<T, Alloc>& src, cudaStream_t stream, size_t destbegin = 0, size_t srcbegin = 0, size_t len = static_cast<uint32_t>(-1))
{
	if (len == static_cast<uint32_t>(-1))
	{
		len = src.size();
	}
	if (dest.size() < destbegin + len)
	{
		throw std::exception();
	}

	return cudaMemcpyAsync(dest.data().get() + destbegin, src.data() + srcbegin, len * sizeof(T), cudaMemcpyHostToDevice, stream);
}

