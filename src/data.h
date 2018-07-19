#pragma once
#include <string>
#include <memory>
#include <utility>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <functional>
#include "util.h"
#include "types.h"

using size_t = std::size_t;

namespace sr
{
namespace data
{

	/// Deathflags guide:
	/// High byte:
	///     Planet ID that killed the particle
	/// Low byte:
	/// 0x80 Particle fell into sun
	/// 0x08 Particle absorbed by planet
	/// 0x04 Kepler didn't converge
	/// 0x02 Out of bounds
	/// 0x01 Particle is in close encounter
	
	struct HostParticleSnapshot
	{
		size_t n, n_alive;

		Vf64_3 r, v;
		Vu32 id;

		void gather(const std::vector<size_t>& indices, size_t begin, size_t length);
		std::unique_ptr<std::vector<size_t>> sort_by_id(size_t begin, size_t length);

		inline HostParticleSnapshot() { }
		inline HostParticleSnapshot(size_t n_) : n(n_), n_alive(n_), r(n_), v(n_), id(n_) { }
	};

	struct HostPlanetSnapshot
	{
		size_t n, n_alive;

		Vf64_3 r, v;
		Vu32 id;
		Vf64 m;

		inline HostPlanetSnapshot() { }
		inline HostPlanetSnapshot(size_t n_) : n(n_), n_alive(n_), r(n_), v(n_), id(n_), m(n_) { }
	};


	struct HostParticlePhaseSpace
	{
		HostParticleSnapshot base;

		inline Vf64_3& r() { return base.r; }
		inline const Vf64_3& r() const { return base.r; }

		inline Vf64_3& v() { return base.v; }
		inline const Vf64_3& v() const { return base.v; }

		inline Vu32& id() { return base.id; }
		inline const Vu32& id() const { return base.id; }

		inline size_t& n() { return base.n; }
		inline const size_t& n() const { return base.n; }

		inline size_t& n_alive() { return base.n_alive; }
		inline const size_t& n_alive() const { return base.n_alive; }

		inline size_t& n_encounter() { return _n_encounter; }
		inline const size_t& n_encounter() const { return _n_encounter; }

		inline Vu16& deathflags() { return _deathflags; }
		inline const Vu16& deathflags() const { return _deathflags; }

		inline Vf32& deathtime() { return _deathtime; }
		inline const Vf32& deathtime() const { return _deathtime; }

		inline Vu32& deathtime_index() { return _deathtime_index; }
		inline const Vu32& deathtime_index() const { return _deathtime_index; }

		inline HostParticlePhaseSpace() { }
		inline HostParticlePhaseSpace(size_t siz, bool cpu_only) : base(siz), _n_encounter(0), _deathflags(siz), _deathtime(siz)
		{ 
			if (cpu_only)
			{
				_deathtime_index = Vu32(siz);
			}
		}

		std::unique_ptr<std::vector<size_t>> stable_partition_alive(size_t begin, size_t length);
		std::unique_ptr<std::vector<size_t>> stable_partition_unflagged(size_t begin, size_t length);
		void gather(const std::vector<size_t>& indices, size_t begin, size_t length);
		std::unique_ptr<std::vector<size_t>> sort_by_id(size_t begin, size_t length);

		inline static uint8_t encounter_planet(uint16_t deathflags)
		{
			return static_cast<uint8_t>((deathflags & 0xFF00) >> 8);
		}

	private:
		size_t _n_encounter;

		Vu16 _deathflags;
		Vf32 _deathtime;

		Vu32 _deathtime_index;
	};

	struct HostPlanetPhaseSpace
	{
		HostPlanetSnapshot base;

		inline Vf64_3& r() { return base.r; }
		inline const Vf64_3& r() const { return base.r; }

		inline Vf64_3& v() { return base.v; }
		inline const Vf64_3& v() const { return base.v; }

		inline Vu32& id() { return base.id; }
		inline const Vu32& id() const { return base.id; }

		inline Vf64& m() { return base.m; }
		inline const Vf64& m() const { return base.m; }

		inline size_t& n() { return base.n; }
		inline const size_t& n() const { return base.n; }

		inline size_t& n_alive() { return base.n_alive; }
		inline const size_t& n_alive() const { return base.n_alive; }

		inline size_t& n_alive_old() { return _n_alive_old; }
		inline const size_t& n_alive_old() const { return _n_alive_old; }

		inline sr::util::LogQuartet<Vf64_3>& r_log() { return _r_log; }
		inline const sr::util::LogQuartet<Vf64_3>& r_log() const { return _r_log; }

		inline sr::util::LogQuartet<Vf64_3>& v_log() { return _v_log; }
		inline const sr::util::LogQuartet<Vf64_3>& v_log() const { return _v_log; }

		inline HostPlanetPhaseSpace() { }
		inline HostPlanetPhaseSpace(size_t siz, size_t tb_size, size_t ce_factor) :
			base(siz), _n_alive_old(siz)
		{
			_r_log = sr::util::LogQuartet<Vf64_3>(((n() - 1) * tb_size), ce_factor);
			_v_log = sr::util::LogQuartet<Vf64_3>(((n() - 1) * tb_size), ce_factor);
		}

		inline void swap_logs()
		{
			std::swap(n_alive(), n_alive_old());

			r_log().swap_logs();
			v_log().swap_logs();
		}

		template<bool old>
		inline size_t log_index_at(size_t timestep, size_t planet_id) const
		{
			return timestep * ((old ? n_alive() : _n_alive_old) - 1) + planet_id - 1;
		}

	private:
		size_t _n_alive_old;

		sr::util::LogQuartet<Vf64_3> _r_log;
		sr::util::LogQuartet<Vf64_3> _v_log;
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
		double t_0, t_f, dt, big_g;
		uint32_t num_thread;
		uint32_t tbsize, print_every, dump_every, track_every, energy_every, max_particle;
		double wh_ce_r1, wh_ce_r2;
		uint32_t wh_ce_n1, wh_ce_n2;
		uint32_t split_track_file;

		bool use_gpu;

		double cull_radius;

		bool resolve_encounters, readmomenta, writemomenta, trackbinary, readsplit, writesplit, dumpbinary, writebinary, readbinary;

		std::string icsin, plin, hybridin, hybridout;
		std::string outfolder;

		inline uint32_t fast_timestep_factor() const
		{
			return resolve_encounters ? (wh_ce_n1 * wh_ce_n2) : 1;
		}

		Configuration();

		Configuration output_config() const;
		inline static Configuration create_dummy()
		{
			Configuration config;
			config.resolve_encounters = false;
			config.tbsize = 0;
			config.wh_ce_n1 = 0;
			config.wh_ce_n2 = 0;

			return config;
		}
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
	inline T read_binary(std::istream& i)
	{
		T t;
		i.read(reinterpret_cast<char*>(&t), sizeof(T));
		return to_little_endian(t);
	}

	template<typename T>
	inline void read_binary(std::istream& i, T& t)
	{
		i.read(reinterpret_cast<char*>(&t), sizeof(T));
		t = to_little_endian(t);
	}

	template<typename T>
	void gather(std::vector<T>& values, const std::vector<size_t>& indices, size_t begin, size_t length)
	{
		std::vector<T> copy(values.begin() + begin, values.begin() + begin + length);
		for (size_t i = begin; i < begin + length; i++)
		{
			values[i] = copy[indices[i - begin]];
		}
	}

	size_t stable_partition_alive_indices(const std::vector<uint16_t>& flags, size_t begin, size_t length, std::unique_ptr<std::vector<size_t>>* indices);

	bool load_data_hybrid_binary(HostData& hd, const Configuration& config, std::istream& in);
	bool load_data_hybrid(HostData& hd, const Configuration& config, std::istream& in);
	bool load_data_split(HostData& hd, const Configuration& config, std::istream& plin, std::istream& icsin);

	bool load_data(HostData& hd, const Configuration& config);

	void save_data_hybrid_binary(const HostData& hd, const Configuration& config, std::ostream& out);
	void save_data_hybrid(const HostData& hd, const Configuration& config, std::ostream& out);
	void save_data_split(const HostData& hd, const Configuration& config, std::ostream& plout, std::ostream& icsout);

	void save_data(const HostData& hd, const Configuration& config, const std::string& outfile, bool dump=false);

	void read_configuration(std::istream& in, Configuration* out);
	void write_configuration(std::ostream& in, const Configuration& config);

	void save_binary_track(std::ostream& trackout, const HostPlanetSnapshot& pl, const HostParticleSnapshot& pa, double time, bool to_elements);
	void load_binary_track(std::istream& trackin, HostPlanetSnapshot& pl, HostParticleSnapshot& pa, double& time, bool skipplanets, bool skipparticles);
	void process_track(std::istream& input, bool takeallparticles, const std::vector<uint32_t>& particle_filter, bool removeplanets,
		const std::function<void(HostPlanetSnapshot&, HostParticleSnapshot&, double)>& callback);

	void read_tracks(const std::string& path, bool takeallparticles, const std::vector<uint32_t>& particle_filter, bool removeplanets,
		const std::function<void(HostPlanetSnapshot&, HostParticleSnapshot&, double)>& callback);

	const size_t TRACK_PARTICLE_STRIDE = 28;
	const size_t TRACK_PLANET_STRIDE = 28;
}
}
