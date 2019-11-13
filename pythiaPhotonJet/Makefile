CXX=clang++
CXXFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(ATLAS_PATH)/include -I$(ROOT_UTIL_PATH) -I$(PYTHIA8_DIR)/include
LDFLAGS=-Wl,-rpath,`root-config --libdir` `root-config --glibs` -L$(ATLAS_PATH)/lib -L$(ROOT_UTIL_PATH) -L$(PYTHIA8_DIR)/lib -lGlobalParams -lUtilities -lAtlas -ldl

all: gammajetgen analyze_gammajet wgen analyze_w gammajetgen_RHIC analyze_gammajet_shift wgen_RHIC analyze_w_shift

gammajetgen: src/gammajetgen.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/gammajetgen.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/gammajetgen

gammajetgen_RHIC: src/gammajetgen_RHIC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/gammajetgen_RHIC.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/gammajetgen_RHIC

wgen: src/wgen.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/wgen.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/wgen

wgen_RHIC: src/wgen_RHIC.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/wgen_RHIC.cc $(PYTHIA8_DIR)/lib/libpythia8.a -o lib/wgen_RHIC

analyze_gammajet: src/analyze_gammajet.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/analyze_gammajet.cc -o lib/analyze_gammajet

analyze_gammajet_shift: src/analyze_gammajet_shift.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/analyze_gammajet_shift.cc -o lib/analyze_gammajet_shift

analyze_w: src/analyze_w.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/analyze_w.cc -o lib/analyze_w

analyze_w_shift: src/analyze_w_shift.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/analyze_w_shift.cc -o lib/analyze_w_shift

clean :
	rm -rf lib/*
