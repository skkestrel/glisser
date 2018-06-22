#pragma once
#include <vector>
#include <cstdint>
#include <ostream>

#define M_PI 3.14159265358979323846
#define M_2PI M_PI * 2

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

using float64_t = double;
using float32_t = float;

using Vf64 = std::vector<float64_t>;
using Vu32 = std::vector<uint32_t>;
using Vu16 = std::vector<uint16_t>;
using Vu8 = std::vector<uint8_t>;
using Vf32 = std::vector<float32_t>;

template<typename T>
struct v_3
{
	T x, y, z;

	__host__ __device__
	inline v_3(const v_3<T>& other) : x(other.x), y(other.y), z(other.z) { }
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

template<typename T>
inline std::ostream& operator<<(std::ostream& stream, const v_3<T>& v)
{
	return stream << v.x << " " << v.y << " " << v.z;
}

using f64_3 = v_3<float64_t>;
using Vf64_3 = std::vector<f64_3>;
