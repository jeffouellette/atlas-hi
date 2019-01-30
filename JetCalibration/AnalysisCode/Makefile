CXX=clang++
CXXFLAGS=-Ofast -fPIC -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I./include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L./lib -lGlobalParams -lUtilities -lTrigger -lAtlas -ldl

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

libs = CalibUtils TreeVariables

analyses = EMTopoComparison ElectronContaminationStudy EnergyScaleChecks FCalDistribution Intercalibration JetAnalysis JetInsituCorrectionCheck MCValidator PhotonAnalysis PhotonJets RtrkComparison ZeeJets ZmumuJets ZMassCalc

hists = $(analyses:%=%Hist) ZJetsHist ZGammaJetCrossCheckHist PhotonRtrkComparisonHist

all : $(libs) $(analyses:%=lib/lib%.so) $(hists:%=lib/lib%.so) 

# Utilities can be directly compiled straight to libraries since they do not depend on other libraries compiled here
lib/lib%.so: obj/%.o
	$(CXX) $(LDFLAGS) -shared $< -o $@

CalibUtils:
	$(CC) -shared src/CalibUtils.cxx -o lib/libCalibUtils.so

RooZfit:
	$(CC) -shared -lRooFit -lRooFitCore src/RooZfit.cxx -o lib/libRooZfit.so

TreeVariables:
	$(CC) -shared src/TreeVariables.cxx -o lib/libTreeVariables.so

# Libraries need to be compiled into object files and then migrated to a shared library
lib/lib%.so: obj/%.o
	$(CXX) $(LDFLAGS) -lCalibUtils -lTreeVariables -shared $< -o $@

obj/%.o :	src/%.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

% : obj/%.o
	$(CXX) $(LDFLAGS) -lCalibUtils -lTreeVariables -shared $< -o lib/lib$@.so

clean :
	rm -rf ./lib/*.so*
