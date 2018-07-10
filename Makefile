.PHONY: all convert
WFLAGS = -Wcast-align \
-Wshadow \
-Wcast-qual -Wconversion \
-Wdisabled-optimization \
-Wfloat-equal -Wformat=2 \
-Wformat-nonliteral -Wformat-security  \
-Wformat-y2k \
-Wimport  -Winit-self \
-Winvalid-pch   \
-Wmissing-field-initializers -Wmissing-format-attribute   \
-Wmissing-include-dirs -Wmissing-noreturn \
-Wpointer-arith \
-Wstack-protector \
-Wstrict-aliasing=2 -Wswitch-default \
-Wswitch-enum \
-Wunreachable-code -Wunused \
-Wunused-parameter \
-Wvariadic-macros \
-Wwrite-strings

default:
	mkdir -p bin
	nvcc targets/main.cpp *.cpp *.cu -lineinfo -g -maxrregcount 64 -arch=sm_35 --std=c++11 --compiler-options "-Wall -Wextra ${WFLAGS} -fstack-protector" -o bin/sr -O3


#-fsanitize=address
cpu:
	mkdir -p bin
	g++ targets/main.cpp *.cpp -g -DNO_CUDA --std=c++11 -Wall -Wextra -Wpedantic ${WFLAGS} -o bin/sr_cpu -O3

convert:
	mkdir -p bin
	g++ targets/convert.cpp *.cpp -g --std=c++11 -Wall -Wextra -Wpedantic ${WFLAGS} -o bin/convert -fsanitize=address -O3

prune-track:
	mkdir -p bin
	g++ targets/prune_track.cpp *.cpp -g --std=c++11 -Wall -Wextra -Wpedantic ${WFLAGS} -o bin/prune-track -fsanitize=address -O3
