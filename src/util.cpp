#include "util.h"

#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

namespace sr
{
namespace util
{
	PathType get_path_type(const std::string& path)
	{
		struct stat s;
		if (stat(path.c_str(), &s) == 0)
		{
			if (s.st_mode & S_IFDIR)
			{
				return PathType::Directory;
			}
			else if (s.st_mode & S_IFREG)
			{
				return PathType::File;
			}
			else
			{
				return PathType::Other;
			}
		}
		else
		{
			return PathType::DoesNotExist;
		}
	}
	bool does_file_exist(const std::string& filename)
	{
		return get_path_type(filename) == PathType::File;
	}

	bool is_dir_empty(const std::string& dirname)
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

	void make_dir(const std::string& path)
	{
		mkdir(path.c_str(), ACCESSPERMS);
	}
}
}
