.PHONY: all clean

CUDA_ARCH=sm_35

SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj
TARGETS_DIR = targets
DOCOPT_DIR = docopt
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES)) $(OBJ_DIR)/docopt/docopt.o
DEP_FILES = $(OBJ_FILES:%.o=%.d)

WFLAGS = -Wcast-align -Wshadow -Wcast-qual -Wconversion -Wdisabled-optimization \
-Wfloat-equal -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k \
-Wimport  -Winit-self -Winvalid-pch -Wmissing-field-initializers -Wmissing-format-attribute   \
-Wmissing-include-dirs -Wmissing-noreturn -Wpointer-arith -Wstack-protector \
-Wstrict-aliasing=2 -Wswitch-default -Wswitch-enum -Wunreachable-code -Wunused \
-Wunused-parameter -Wvariadic-macros -Wwrite-strings -Wunused-result

CPPFLAGS = ${WFLAGS} -g --std=c++11 -Wall -Wextra -Wpedantic ${WFLAGS} -DNO_CUDA -O3 # -fsanitize=address

LDFLAGS = -pthread # -lasan

glisser:
	@mkdir -p $(BIN_DIR)
	@nvcc $(TARGETS_DIR)/main.cpp $(DOCOPT_DIR)/docopt.cpp $(SRC_FILES) $(SRC_DIR)/*.cu -lineinfo -g -arch=$(CUDA_ARCH) -maxrregcount 64 --std=c++11 -D_GLIBC_USE_C99 --compiler-options "-Wall -Wextra ${WFLAGS} -fstack-protector" -o $(BIN_DIR)/glisser -O3 -lnvToolsExt

clean:
	rm -r $(OBJ_DIR)/* 
	rm -r $(BIN_DIR)/*

$(OBJ_DIR)/docopt/docopt.o: $(DOCOPT_DIR)/docopt.cpp
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/docopt
	$(info CPP $(OBJ_DIR)/docopt/docopt.o)
	@g++ $(DOCOPT_DIR)/docopt.cpp $(CPPFLAGS) -c -o $(OBJ_DIR)/docopt/docopt.o

-include $(DEP_FILES)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p obj
	$(info CPP $@)
	@g++ $(CPPFLAGS) -MMD -c -o $@ $<

define make-target =
@mkdir -p $(BIN_DIR)
@mkdir -p $(OBJ_DIR)/targets
$(info CPP $(OBJ_DIR)/targets/$1.o)
@g++ $(TARGETS_DIR)/$1.cpp $(CPPFLAGS) -c -o $(OBJ_DIR)/targets/$1.o
$(info LD $(BIN_DIR)/$2)
@g++ $(OBJ_DIR)/targets/$1.o $(OBJ_FILES) $(LDFLAGS) -o $(BIN_DIR)/$2
endef

convert-state: $(OBJ_FILES)
	$(call make-target,convert_state,convert-state)

filter-state: $(OBJ_FILES)
	$(call make-target,filter_state,filter-state)

prune-track: $(OBJ_FILES)
	$(call make-target,prune_track,prune-track)

find-max-e: $(OBJ_FILES)
	$(call make-target,find_max_e,find-max-e)

make-state: $(OBJ_FILES)
	$(call make-target,make_state,make-state)

find-librators: $(OBJ_FILES)
	$(call make-target,find_librators,find-librators)

track-info: $(OBJ_FILES)
	$(call make-target,track_info,track-info)

export-track: $(OBJ_FILES)
	$(call make-target,export_track,export-track)

export-swift: $(OBJ_FILES)
	$(call make-target,export_swift,export-swift)

import-swift: $(OBJ_FILES)
	$(call make-target,import_swift,import-swift)

lookup-info: $(OBJ_FILES)
	$(call make-target,lookup_info,lookup-info)

cj-metrics: $(OBJ_FILES)
	$(call make-target,cj_metrics,cj-metrics)
