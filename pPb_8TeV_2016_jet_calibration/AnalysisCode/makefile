#ANALY_PATH=$(CALIB_PATH)/AnalysisCode

CXXFLAGS=-Ofast -fPIC -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I./include
#LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L./lib -lGlobalParams -lInitialization -lTrigger -lAtlas -ldl
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L./lib -lGlobalParams -lTrigger -lAtlas -ldl

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

#hists = EMTopoComparisonHist ElectronContaminationStudyHist EnergyScaleChecksHist FCalDistributionHist GammaJetsHist JetInsituCorrectionCheckHist RtrkComparisonHist ZeeJetsHist ZGammaJetCrossCheckHist ZmumuJetsHist ZMassCalcHist

analyses = EMTopoComparison ElectronContaminationStudy EnergyScaleChecks FCalDistribution GammaJets JetInsituCorrectionCheck RtrkComparison ZeeJets ZGammaJetCrossCheck ZmumuJets ZMassCalc

all : Utils TreeVariables $(analyses:%=lib/lib%.so) $(analyses:%=lib/lib%Hist.so)

# Utilities can be directly compiled straight to libraries since they do not depend on other libraries compiled here
lib/lib%.so: obj/%.o Utils TreeVariables
	$(CXX) $(LDFLAGS) -lUtils -lTreeVariables -shared $< -o $@

Utils:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared src/Utils.C -o lib/libUtils.so

TreeVariables:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared src/TreeVariables.C -o lib/libTreeVariables.so

RooZfit:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -lRooFit -lRooFitCore src/RooZfit.C -o lib/libRooZfit.so

# Libraries need to be compiled into object files and then migrated to a shared library
lib/lib%.so: obj/%.o Utils TreeVariables
	$(CXX) $(LDFLAGS) -lUtils -lTreeVariables -shared $< -o $@

obj/%.o :	src/%.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

% : obj/%.o
	$(CXX) $(LDFLAGS) -lUtils -lTreeVariables -shared $< -o lib/lib$@.so

clean :
	rm -rf ./lib/*.so*
