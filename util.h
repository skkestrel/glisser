#pragma once
#include <streambuf>
#include <dirent.h>
#include <unistd.h>
#include <algorithm>

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

inline bool does_file_exist(const std::string& filename)
{
	return access(filename.c_str(), F_OK) != -1;
}


inline bool is_dir_empty(const std::string& dirname)
{
	int n = 0;
	DIR *dir = opendir(dirname.c_str());
	if (!dir) return 1;
	while (readdir(dir))
	{
		if (++n > 2) break;
	}
	closedir(dir);
	return n <= 2;
}

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
	Vector log;
	Vector old;
	Vector slow;
	Vector slow_old;

	inline LogQuartet(size_t slow_size, size_t speed_factor)
	{
		slow = slow_old = Vector(slow_size);

		if (speed_factor > 1)
		{
			log = old = Vector(slow_size * speed_factor);
		}
	}

	inline void swap_old()
	{
		std::swap(log, old);
		std::swap(slow, slow_old);
	}

	template<typename... Args>
	LogQuartet(Args... args);

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

template<typename Vector>
template<typename... Args>
inline LogQuartet<Vector>::LogQuartet(Args... args)
	: log(args...), old(args...), slow(args...), slow_old(args...)
{ }
