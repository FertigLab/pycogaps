CXX := g++
CXXFLAGS := -fPIC -O2 -Wall -std=c++11
PY_INCLUDES := $(shell python3 -m pybind11 --includes)
INCLUDES := -Isrc/Rpackage/src/include $(PY_INCLUDES)
COGAPS_SOURCES := \
	$(wildcard src/Rpackage/src/*.cpp) \
	$(wildcard src/Rpackage/src/atomic/*.cpp) \
	$(wildcard src/Rpackage/src/data_structures/*.cpp) \
	$(wildcard src/Rpackage/src/file_parser/*.cpp) \
	$(wildcard src/Rpackage/src/gibbs_sampler/*.cpp) \
	$(wildcard src/Rpackage/src/math/*.cpp) \
	$(wildcard src/Rpackage/src/utils/*.cpp)

COGAPS_SOURCES := $(filter-out src/Rpackage/src/Cogaps.cpp, $(COGAPS_SOURCES))
COGAPS_SOURCES := $(filter-out src/Rpackage/src/RcppExports.cpp, $(COGAPS_SOURCES))
COGAPS_SOURCES := $(filter-out src/Rpackage/src/test-runner.cpp, $(COGAPS_SOURCES))

OBJECTS := \
	$(COGAPS_SOURCES:src/Rpackage/src/%.cpp=build/libcogaps/%.o) \
	build/CogapsPy.o

all : CogapsPy.so

CogapsPy.so : $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@ -o $@ $^

build/%.o : src/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/libcogaps/%.o : src/Rpackage/src/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/libcogaps/atomic/%.o : src/Rpackage/src/atomic/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/atomic
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/libcogaps/data_structures/%.o : src/Rpackage/src/data_structures/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/data_structures
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/libcogaps/file_parser/%.o : src/Rpackage/src/file_parser/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/file_parser
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/libcogaps/gibbs_sampler/%.o : src/Rpackage/src/gibbs_sampler/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/gibbs_sampler
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/libcogaps/math/%.o : src/Rpackage/src/math/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/math
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/libcogaps/utils/%.o : src/Rpackage/src/utils/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/utils
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean : 
	rm -rf build
