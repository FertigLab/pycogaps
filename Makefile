CXX := g++
CXXFLAGS := -fPIC
LIBS := -lpython2.7 -lboost_python

PY_SOURCES := \
	src/HelloWorld.cpp

COGAPS_SOURCES := \
	src/Rpackage/src/AtomicDomain.cpp \
	src/Rpackage/src/GapsRunner.cpp \
	src/Rpackage/src/GapsStatistics.cpp \
	src/Rpackage/src/GibbsSampler.cpp \
	src/Rpackage/src/ProposalQueue.cpp \
	src/Rpackage/src/data_structures/Matrix.cpp \
	src/Rpackage/src/data_structures/Vector.cpp \
	src/Rpackage/src/file_parser/CsvParser.cpp \
	src/Rpackage/src/file_parser/FileParser.cpp \
	src/Rpackage/src/file_parser/MtxParser.cpp \
	src/Rpackage/src/file_parser/TsvParser.cpp \
	src/Rpackage/src/math/Algorithms.cpp \
	src/Rpackage/src/math/Random.cpp

OBJECTS := \
	$(PY_SOURCES:src/%.cpp=build/%.o) \
	$(COGAPS_SOURCES:src/Rpackage/src/%.cpp=build/libcogaps/%.o)

all :
	@echo $(OBJECTS)

CogapsPy.so : $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@ -o $@ $^ $(LIBS)

build/%.o : src/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -c $< -o $@

build/libcogaps/%.o : src/Rpackage/src/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/libcogaps/data_structures/%.o : src/Rpackage/src/data_structures/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/data_structures
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/libcogaps/file_parser/%.o : src/Rpackage/src/file_parser/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/file_parser
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/libcogaps/math/%.o : src/Rpackage/src/math/%.cpp
	@mkdir -p build
	@mkdir -p build/libcogaps
	@mkdir -p build/libcogaps/math
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean : 
	rm -rf build