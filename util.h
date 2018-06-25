#include <streambuf>
#include <dirent.h>
#include <unistd.h>

inline bool does_file_exist(const std::string& filename)
{
	return access(filename.c_str(), F_OK) != -1;
}


inline bool is_dir_empty(const std::string& dirname)
{
	int n = 0;
	struct dirent *d;
	DIR *dir = opendir(dirname.c_str());
	if (!dir) return 1;
	while ((d = readdir(dir)))
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
