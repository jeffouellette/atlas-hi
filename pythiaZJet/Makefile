CXX=clang++
CXXFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(ATLAS_PATH)/include -I$(ROOT_UTIL_PATH) -I$(PYTHIA8_DIR)/include
LDFLAGS=-Wl,-rpath,`root-config --libdir` `root-config --glibs` -L$(ATLAS_PATH)/lib -L$(ROOT_UTIL_PATH) -L$(PYTHIA8_DIR)/lib -lGlobalParams -lUtilities -lAtlas -ldl

all: gammajetgen_LHC analyze_gammajet_LHC gammajetgen_RHIC analyze_gammajet_RHIC zgen_LHC analyze_z_LHC

gammajetgen_LHC: src/gammajetgen_LHC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/gammajetgen_LHC.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/gammajetgen_LHC

gammajetgen_RHIC: src/gammajetgen_RHIC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/gammajetgen_RHIC.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/gammajetgen_RHIC

zgen_LHC: src/zgen_LHC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/zgen_LHC.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/zgen_LHC

mbgen_LHC: src/mbgen_LHC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/mbgen_LHC.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/mbgen_LHC

analyze_gammajet_LHC: src/analyze_gammajet_LHC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/analyze_gammajet_LHC.cc -o lib/analyze_gammajet_LHC

analyze_gammajet_RHIC: src/analyze_gammajet_RHIC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/analyze_gammajet_RHIC.cc -o lib/analyze_gammajet_RHIC

analyze_z_LHC: src/analyze_z_LHC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/analyze_z_LHC.cc -o lib/analyze_z_LHC

clean :
	rm -rf lib/*
