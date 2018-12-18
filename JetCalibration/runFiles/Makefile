CXXFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I$(CALIB_PATH)/AnalysisCode/include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L$(CALIB_PATH)/AnalysisCode/lib -lAtlas -lGlobalParams -lInitialization -lTrigger -ldl

analyses = runEMTopoComparison runElectronContaminationStudy runEnergyScaleChecks runFCalDistribution runJetAnalysis runJetInsituCorrectionCheck runMCValidator runPhotonAnalysis runPhotonJets runRtrkComparison runSubtractionStudy runZeeJets runZmumuJets runZMassCalc

hists = $(analyses:%=%Hist) runZGammaJetCrossCheckHist runZJetsHist runPhotonRtrkComparisonHist

all : $(hists) $(analyses)

run%: obj/run%.o
	$(CXX) $(LDFLAGS) -l$(@:run%=%) obj/$@.o -o bin/$@

obj/run%.o:
	$(CXX) $(CXXFLAGS) -c $(@:obj/run%.o=src/run%.cc) -o $@

clean :
	rm -rf ./bin/*
