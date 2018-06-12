#pragma once
#include <cstdint>
#include "thrust.h"

#define M_PI 3.14159265358979323846
#define M_2PI M_PI * 2



using float64_t = double;
using float32_t = float;

using Hvf64 = thrust::host_vector<float64_t>;
using Hvu32 = thrust::host_vector<uint32_t>;
using Hvu8 = thrust::host_vector<uint8_t>;
using Hvf32 = thrust::host_vector<float32_t>;

using Dvf64 = thrust::device_vector<float64_t>;
using Dvu32 = thrust::device_vector<uint32_t>;
using Dvu8 = thrust::device_vector<uint8_t>;
using Dvf32 = thrust::device_vector<float32_t>;

template<typename T>
struct v_3
{
	T x, y, z;

	__host__ __device__
	inline v_3() { }
	__host__ __device__
	explicit inline v_3(T t) : x(t), y(t), z(t) { }
	__host__ __device__
	inline v_3(T x, T y, T z) : x(x), y(y), z(z) { }

	__host__ __device__
	inline v_3 cross(const v_3<T>& b) const
	{
		return v_3<T>(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
	}

	__host__ __device__
	inline T lensq() const
	{
		return x * x + y * y + z * z;
	}
	
	__host__ __device__
	inline const v_3<T>& operator*=(T a)
	{
		x *= a;
		y *= a;
		z *= a;
		return *this;
	}

	__host__ __device__
	inline const v_3<T>& operator/=(T a)
	{
		x /= a;
		y /= a;
		z /= a;
		return *this;
	}

	__host__ __device__
	inline const v_3<T>& operator/=(const v_3<T>& a)
	{
		x /= a.x;
		y /= a.y;
		z /= a.z;
		return *this;
	}

	__host__ __device__
	inline const v_3<T>& operator*=(const v_3<T>& a)
	{
		x *= a.x;
		y *= a.y;
		z *= a.z;
		return *this;
	}

	__host__ __device__
	inline const v_3<T>& operator-=(const v_3<T>& a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;
		return *this;
	}

	__host__ __device__
	inline const v_3<T>& operator+=(const v_3<T>& a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	}

	__host__ __device__
	inline v_3<T> operator+(const v_3<T>& b) const
	{
		return v_3<T>(x + b.x, y + b.y, z + b.z);
	}

	__host__ __device__
	inline v_3<T> operator-(const v_3<T>& b) const
	{
		return v_3<T>(x - b.x, y - b.y, z - b.z);
	}

	__host__ __device__
	inline v_3<T> operator-() const
	{
		return v_3<T>(-x, -y, -z);
	}

	__host__ __device__
	inline v_3<T> operator/(const v_3<T>& b) const
	{
		return v_3<T>(x / b.x, y / b.y, z / b.z);
	}

	__host__ __device__
	inline v_3<T> operator*(const v_3<T>& b) const
	{
		return v_3<T>(x * b.x, y * b.y, z * b.z);
	}

	__host__ __device__
	inline v_3<T> operator/(T b) const
	{
		return v_3<T>(x / b, y / b, z / b);
	}

	__host__ __device__
	inline v_3<T> operator*(T b) const
	{
		return v_3<T>(x * b, y * b, z * b);
	}

	__host__ __device__
	inline v_3<T> operator*(v_3<T>& b) const
	{
		return v_3<T>(x * b, y * b, z * b);
	}
};

using f64_3 = v_3<float64_t>;

using Hvf64_3 = thrust::host_vector<f64_3>;
using Dvf64_3 = thrust::device_vector<f64_3>;
