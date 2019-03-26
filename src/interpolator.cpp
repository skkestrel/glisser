#include <data.h>
#include <types.h>

namespace sr
{
	Interpolator::Interpolator(std::HostData& hd, std::string file)
	{
		input = std::ostream(file, std::ios_base::binary);

		hd.planets.m()[0] = sr::data::read_binary<float64_t>(input);
		sr::data::skip_binary(input, 32 - 8);

		t0 = sr::data::read_binary<float64_t>(input);

		hd.planets.n() = sr::data::read_binary<uint32_t>(input) + 1;
		hd.planets.n()
	}
}
