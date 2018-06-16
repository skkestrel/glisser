#pragma once
#include "thrust.cuh"
#include "types.h"

using Hvf64 = thrust::host_vector<float64_t>;
using Hvu32 = thrust::host_vector<uint32_t>;
using Hvu16 = thrust::host_vector<uint16_t>;
using Hvu8 = thrust::host_vector<uint8_t>;
using Hvf32 = thrust::host_vector<float32_t>;

using Dvf64 = thrust::device_vector<float64_t>;
using Dvu32 = thrust::device_vector<uint32_t>;
using Dvu16 = thrust::device_vector<uint16_t>;
using Dvu8 = thrust::device_vector<uint8_t>;
using Dvf32 = thrust::device_vector<float32_t>;

using Hvf64_3 = thrust::host_vector<f64_3>;
using Dvf64_3 = thrust::device_vector<f64_3>;
