.PHONY: all clean

SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj
DOCOPT_DIR = docopt
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES)) $(DOCOPT_DIR)/docopt.o

WFLAGS = -Wcast-align -Wshadow -Wcast-qual -Wconversion -Wdisabled-optimization \
-Wfloat-equal -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k \
-Wimport  -Winit-self -Winvalid-pch -Wmissing-field-initializers -Wmissing-format-attribute   \
-Wmissing-include-dirs -Wmissing-noreturn -Wpointer-arith -Wstack-protector \
-Wstrict-aliasing=2 -Wswitch-default -Wswitch-enum -Wunreachable-code -Wunused \
-Wunused-parameter -Wvariadic-macros -Wwrite-strings

CPPFLAGS = ${WFLAGS} -g --std=c++11 -Wall -Wextra -Wpedantic ${WFLAGS} -fsanitize=address -O3 -DNO_CUDA


LDFLAGS = -lasan

clean:
	rm $(OBJ_DIR)/*
	rm $(BIN_DIR)/*

default:
	mkdir -p $(BIN_DIR)
	nvcc targets/main.cpp *.cpp *.cu -lineinfo -g -maxrregcount 64 -arch=sm_35 --std=c++11 --compiler-options "-Wall -Wextra ${WFLAGS} -fstack-protector" -o $(BIN_DIR)/sr -O3

$(OBJ_DIR)/docopt.o:
	@mkdir -p obj
	@g++ $(DOCOPT_DIR)/docopt.cpp $(CPPFLAGS) -c -o $(OBJ_DIR)/docopt.o
	$(info CPP $(DOCOPT_DIR)/docopt.cpp)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p obj
	@g++ $(CPPFLAGS) -c -o $@ $<
	$(info CPP $<)

cpu: $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	@g++ targets/main.cpp $(OBJ_FILES) $(CPPFLAGS) $(LDFLAGS) -o $(BIN_DIR)/sr_cpu
	$(info LD $(BIN_DIR)/sr_cpu)

convert-state: $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	@g++ targets/convert_state.cpp $(OBJ_FILES) $(LDFLAGS) $(CPPFLAGS) -o $(BIN_DIR)/convert-state
	$(info LD $(BIN_DIR)/convert-state)

prune-track: $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	@g++ targets/prune_track.cpp $(OBJ_FILES) $(LDFLAGS) $(CPPFLAGS) -o $(BIN_DIR)/prune-track
	$(info LD $(BIN_DIR)/prune-track)

filter-state: $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	@g++ targets/filter_state.cpp $(OBJ_FILES) $(LDFLAGS) $(CPPFLAGS) -o $(BIN_DIR)/filter-state
	$(info LD $(BIN_DIR)/filter-state)

find-max-e: $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	@g++ targets/find_max_e.cpp $(OBJ_FILES) $(LDFLAGS) $(CPPFLAGS) -o $(BIN_DIR)/find-max-e
	$(info LD $(BIN_DIR)/find-max-e)
