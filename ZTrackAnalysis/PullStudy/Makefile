CXX=g++
CXXFLAGS=-O3 -g -Wall `root-config --cflags` -I$(ATLAS_PATH)/include -I$(ROOT_UTIL_PATH) -I$(PYTHIA8_DIR)/include
LDFLAGS=-Wl,-rpath `root-config --libdir --glibs` -L$(ATLAS_PATH)/lib -L$(ROOT_UTIL_PATH) -L$(PYTHIA8_DIR)/lib -lGlobalParams -lUtilities -lAtlasUtils -lAtlasStyle -ldl

reqdirs= bin outputs logs errors Plots

directories:
	mkdir -p ${reqdirs}

all: directories gen analyze

gen: src/gen.cxx
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) $(PYTHIA8_DIR)/lib/libpythia8.a -o bin/gen

analyze: src/analyze.cxx
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) -o bin/analyze

clean:
	rm -rf bin logs errors
