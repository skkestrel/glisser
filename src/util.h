#pragma once
#include <streambuf>
#include <algorithm>
#include <ostream>
#include <memory>

#if __cplusplus < 201404L
namespace std
{
	template<typename T, typename... Args>
	inline std::unique_ptr<T> make_unique(Args&&... args)
	{
	    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
	}
}
#endif

namespace std
{
	inline unsigned stou(std::string const & str, size_t * idx = 0, int base = 10)
	{
		unsigned long result = std::stoul(str, idx, base);
		if (result > std::numeric_limits<unsigned>::max())
		{
			throw std::out_of_range("stou");
		}
		return static_cast<unsigned>(result);
	}
}

namespace sr
{
namespace util
{
	enum class PathType
	{
		File,
		Directory,
		DoesNotExist,
		Other
	};

	PathType get_path_type(const std::string& path);
	bool does_file_exist(const std::string& filename);
	bool is_dir_empty(const std::string& dirname);
	void make_dir(const std::string& path);

	class teebuf : public std::streambuf
	{
		public:
			// Construct a streambuf which tees output to both input
			// streambufs.
			inline teebuf(std::streambuf* _sb1, std::streambuf* _sb2)
				: sb1(_sb1), sb2(_sb2) { }
		private:
			// This tee buffer has no buffer. So every character "overflows"
			// and can be put directly into the teed buffers.
			inline virtual int overflow(int c)
			{
				if (c == EOF)
				{
					return !EOF;
				}
				else
				{
					int r1 = sb1->sputc(static_cast<char>(c));
					int r2 = sb2->sputc(static_cast<char>(c));
					return r1 == EOF || r2 == EOF ? EOF : c;
				}
			}

			// Sync both teed buffers.
			inline virtual int sync()
			{
				int const r1 = sb1->pubsync();
				int const r2 = sb2->pubsync();
				return r1 == 0 && r2 == 0 ? 0 : -1;
			}   
		private:
			std::streambuf* sb1;
			std::streambuf* sb2;
	};

	class teestream : public std::ostream
	{
	public:
	    // Construct an ostream which tees output to the supplied
	    // ostreams.
	    teestream(std::ostream& o1, std::ostream& o2);
	private:
	    teebuf tbuf;
	};

	inline teestream::teestream(std::ostream& o1, std::ostream& o2) : std::ostream(&tbuf), tbuf(o1.rdbuf(), o2.rdbuf()) { }

	namespace detail
	{
		template<typename Self, typename Return, bool slow, bool old>
		struct Selector
		{
			static Return get(Self x)
			{
				return Self::not_implemented;
			}
		};

		template<typename Self, typename Return>
		struct Selector<Self, Return, false, false>
		{
			static Return get(Self x) { return x->log; }
		};

		template<typename Self, typename Return>
		struct Selector<Self, Return, false, true>
		{
			static Return get(Self x) { return x->old; }
		};

		template<typename Self, typename Return>
		struct Selector<Self, Return, true, false>
		{
			static Return get(Self x) { return x->slow; }
		};

		template<typename Self, typename Return>
		struct Selector<Self, Return, true, true>
		{
			static Return get(Self x) { return x->slow_old; }
		};
	}

	template<typename Vector>
	struct LogQuartet
	{
		size_t len;
		size_t len_old;

		Vector log;
		Vector old;
		Vector slow;
		Vector slow_old;

		inline LogQuartet() { }

		inline LogQuartet(size_t slow_size, size_t speed_factor)
		{
			slow = slow_old = Vector(slow_size);

			if (speed_factor > 1)
			{
				log = old = Vector(slow_size * speed_factor);
			}
		}

		inline void swap_logs()
		{
			std::swap(log, old);
			std::swap(slow, slow_old);
			std::swap(len, len_old);
		}

		template<bool slow, bool old>
		inline Vector& get()
		{
			return detail::Selector<LogQuartet<Vector>*, Vector&, slow, old>::get(this);
		}

		template<bool slow, bool old>
		inline const Vector& get() const
		{
			return detail::Selector<const LogQuartet<Vector>*, const Vector&, slow, old>::get(this);
		}
	};

	inline std::string joinpath(const std::string& base, const std::string& app)
	{
		return base + "/" + app;
	}
}
}
