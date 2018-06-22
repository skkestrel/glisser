#pragma once
#include "types.h"
#include <string>
#include <utility>
#include <cstdlib>
#include <istream>
#include <ostream>
using size_t = std::size_t;

struct HostParticlePhaseSpace
{
	size_t n, n_alive, n_encounter;

	Vf64_3 r, v;
	Vf64_3 a;

	Vu16 deathflags;
	Vf32 deathtime; Vu32 id; 
	inline HostParticlePhaseSpace() { }
	inline HostParticlePhaseSpace(size_t n) : n(n), n_alive(n), n_encounter(0), r(n), v(n), a(n), deathflags(n), deathtime(n), id(n) { }

	void stable_partition_alive(size_t begin = 0, size_t length = static_cast<size_t>(-1));
};

struct HostPlanetPhaseSpace
{
	size_t n, n_alive;
	Vf64 m, eta;
	Vf64_3 r, v, rj, vj;
	Vf64_3 a;

	Vf64_3 r_log, v_log;
	Vf64_3 r_log_old, v_log_old;
	Vf64_3 r_log_slow, v_log_slow;
	Vf64_3 r_log_slow_old, v_log_slow_old;
	Vf64_3 h0_log, h0_log_old;
	Vf64_3 h0_log_slow, h0_log_slow_old;

	f64_3 bary_r, bary_v;

	Vu32 id;

	inline HostPlanetPhaseSpace() { }
	inline HostPlanetPhaseSpace(size_t n, size_t tb_size, size_t ce_factor):
		n(n), n_alive(n), m(n), eta(n), r(n), v(n), rj(n), vj(n),
		a(n), id(n)
	{
		r_log = v_log = r_log_old = v_log_old = Vf64_3((n - 1) * tb_size * ce_factor);
		r_log_slow = v_log_slow = r_log_slow_old = v_log_slow_old = Vf64_3((n - 1) * tb_size);

		h0_log = h0_log_old = Vf64_3(tb_size * ce_factor);
		h0_log_slow = h0_log_slow_old = Vf64_3(tb_size);
	}
};

struct HostPlanetSnapshot
{
	size_t n, n_alive;
	Vf64 m;
	Vf64_3 r, v;
	Vu32 id;

	inline HostPlanetSnapshot() { }
	inline HostPlanetSnapshot(const HostPlanetPhaseSpace& o) :
		n(o.n), n_alive(o.n_alive), m(o.m), r(o.r), v(o.v), id(o.id) { }
};

struct HostData
{
	HostParticlePhaseSpace particles;
	HostPlanetPhaseSpace planets;
	HostPlanetSnapshot planets_snapshot;

	inline HostData() { }
};

struct Configuration
{
	double t, t_0, t_f, dt;
	size_t tbsize, ce_factor, print_every, dump_every, track_every, max_particle;
	bool resolve_encounters, readmomenta, trackbinary, readhybrid, writehybrid, dumpbinary, writehybridbinary, readhybridbinary;

	std::string icsin, plin, hybridin;
	std::string outfolder;

	Configuration();
};


template<typename T>
inline T reverse_2byte(T in)
{
	static_assert(sizeof(T) == 2, "");
	T ret;
	char *inb = reinterpret_cast<char*>(&in);
	char *retb = reinterpret_cast<char*>(&ret);

	// swap the bytes into a temporary buffer
	retb[0] = inb[1];
	retb[1] = inb[0];

	return ret;
}

template<typename T>
inline T reverse_4byte(T in)
{
	static_assert(sizeof(T) == 4, "");
	T ret;
	char *inb = reinterpret_cast<char*>(&in);
	char *retb = reinterpret_cast<char*>(&ret);

	// swap the bytes into a temporary buffer
	retb[0] = inb[3];
	retb[1] = inb[2];
	retb[2] = inb[1];
	retb[3] = inb[0];

	return ret;
}

template<typename T>
inline T reverse_8byte(T in)
{
	static_assert(sizeof(T) == 8, "");
	T ret;
	char *inb = reinterpret_cast<char*>(&in);
	char *retb = reinterpret_cast<char*>(&ret);

	// swap the bytes into a temporary buffer
	retb[0] = inb[7];
	retb[1] = inb[6];
	retb[2] = inb[5];
	retb[3] = inb[4];
	retb[4] = inb[3];
	retb[5] = inb[2];
	retb[6] = inb[1];
	retb[7] = inb[0];

	return ret;
}

template<typename T>
inline T reverse_bytes(T in) { (void) in; return T::unimplemented; }

template<>
inline int16_t reverse_bytes<int16_t>(int16_t in) { return reverse_2byte(in); }
template<>
inline uint16_t reverse_bytes<uint16_t>(uint16_t in) { return reverse_2byte(in); }
template<>
inline int32_t reverse_bytes<int32_t>(int32_t in) { return reverse_4byte(in); }
template<>
inline uint32_t reverse_bytes<uint32_t>(uint32_t in) { return reverse_4byte(in); }
template<>
inline int64_t reverse_bytes<int64_t>(int64_t in) { return reverse_8byte(in); }
template<>
inline uint64_t reverse_bytes<uint64_t>(uint64_t in) { return reverse_8byte(in); }
template<>
inline double reverse_bytes<double>(double in) { return reverse_8byte(in); }
template<>
inline float reverse_bytes<float>(float in) { return reverse_4byte(in); }

inline bool is_int_little_endian()
{
	uint32_t i = 0x01;
	return reinterpret_cast<uint8_t*>(&i)[0] == 0x01;
}

inline bool is_float_little_endian()
{
	float i = 2;
	return reinterpret_cast<uint8_t*>(&i)[3] == 0x40;
}

inline bool is_double_little_endian()
{
	double i = 2;
	return reinterpret_cast<uint8_t*>(&i)[7] == 0x40;
}

template<typename T>
inline bool is_little_endian() { return T::unimplemented; }

template<>
inline bool is_little_endian<int16_t>() { return is_int_little_endian(); }
template<>
inline bool is_little_endian<uint16_t>() { return is_int_little_endian(); }
template<>
inline bool is_little_endian<int32_t>() { return is_int_little_endian(); }
template<>
inline bool is_little_endian<uint32_t>() { return is_int_little_endian(); }
template<>
inline bool is_little_endian<float>() { return is_float_little_endian(); }
template<>
inline bool is_little_endian<int64_t>() { return is_int_little_endian(); }
template<>
inline bool is_little_endian<uint64_t>() { return is_int_little_endian(); }
template<>
inline bool is_little_endian<double>() { return is_double_little_endian(); }

template<typename T>
inline T to_little_endian(T in)
{
	if (is_little_endian<T>())
	{
		return in;
	}
	else
	{
		return reverse_bytes(in);
	}
}

template<typename T>
inline void write_binary(std::ostream& o, const T& t)
{
	T c = to_little_endian(t);
	o.write(reinterpret_cast<const char*>(&c), sizeof(c));
}

template<typename T>
inline void read_binary(std::istream& i, T& t)
{
	i.read(reinterpret_cast<char*>(&t), sizeof(T));
	t = to_little_endian(t);
}

bool load_data(HostData& hd, const Configuration& config);
void save_data(const HostData& hd, const Configuration& config, bool dump=false, size_t dumpnum=0);
bool read_configuration(std::istream& in, Configuration* out);
void write_configuration(std::ostream& in, const Configuration& config);
std::string joinpath(const std::string& base, const std::string& file);
