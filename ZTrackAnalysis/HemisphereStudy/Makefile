CXX=clang++
CXXFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(ATLAS_PATH)/include -I$(ROOT_UTIL_PATH) -I$(FASTJET_PATH)/include
LDFLAGS=-Wl,-rpath,`root-config --libdir` `root-config --glibs` -L$(ATLAS_PATH)/lib -L$(ROOT_UTIL_PATH) -L$(FASTJET_PATH)/lib -lGlobalParams -lUtilities -lAtlas -ldl -lfastjet

all: 

RunFastJet: src/RunFastJet.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/RunFastJet.cc -o bin/RunFastJet

clean :
	rm -rf lib/*
