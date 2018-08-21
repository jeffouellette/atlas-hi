CFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I$(CALIB_PATH)/AnalysisCode/include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L$(CALIB_PATH)/AnalysisCode/lib -lAtlas -lGlobalParams -lInitialization -lTrigger -ldl


all : runEMTopoComparison runEMTopoComparisonHist runElectronContaminationStudy runElectronContaminationStudyHist runEnergyScaleChecks runEnergyScaleChecksHist runJetInsituCorrectionCheck runJetInsituCorrectionCheckHist runRtrkComparison runRtrkComparisonHist runZGammaJetCrossCheck runZGammaJetCrossCheckHist

runEMTopoComparison :
	$(CXX) -o `pwd`/lib/runEMTopoComparison `pwd`/src/runEMTopoComparison.cc $(CFLAGS) $(LDFLAGS) -lEMTopoComparison

runEMTopoComparisonHist :
	$(CXX) -o `pwd`/lib/runEMTopoComparisonHist `pwd`/src/runEMTopoComparisonHist.cc $(CFLAGS) $(LDFLAGS) -lEMTopoComparisonHist 

runElectronContaminationStudy :
	$(CXX) -o `pwd`/lib/runElectronContaminationStudy `pwd`/src/runElectronContaminationStudy.cc $(CFLAGS) $(LDFLAGS) -lElectronContaminationStudy

runElectronContaminationStudyHist :
	$(CXX) -o `pwd`/lib/runElectronContaminationStudyHist `pwd`/src/runElectronContaminationStudyHist.cc $(CFLAGS) $(LDFLAGS) -lElectronContaminationStudyHist

runEnergyScaleChecks :
	$(CXX) -o `pwd`/lib/runEnergyScaleChecks `pwd`/src/runEnergyScaleChecks.cc $(CFLAGS) $(LDFLAGS) -lEnergyScaleChecks

runEnergyScaleChecksHist :
	$(CXX) -o `pwd`/lib/runEnergyScaleChecksHist `pwd`/src/runEnergyScaleChecksHist.cc $(CFLAGS) $(LDFLAGS) -lEnergyScaleChecksHist

runJetInsituCorrectionCheck :
	$(CXX) -o `pwd`/lib/runJetInsituCorrectionCheck `pwd`/src/runJetInsituCorrectionCheck.cc $(CFLAGS) $(LDFLAGS) -lJetInsituCorrectionCheck

runJetInsituCorrectionCheckHist :
	$(CXX) -o `pwd`/lib/runJetInsituCorrectionCheckHist `pwd`/src/runJetInsituCorrectionCheckHist.cc $(CFLAGS) $(LDFLAGS) -lJetInsituCorrectionCheckHist

runRtrkComparison :
	$(CXX) -o `pwd`/lib/runRtrkComparison `pwd`/src/runRtrkComparison.cc $(CFLAGS) $(LDFLAGS) -lRtrkComparison

runRtrkComparisonHist :
	$(CXX) -o `pwd`/lib/runRtrkComparisonHist `pwd`/src/runRtrkComparisonHist.cc $(CFLAGS) $(LDFLAGS) -lRtrkComparisonHist

runZGammaJetCrossCheck :
	$(CXX) -o `pwd`/lib/runZGammaJetCrossCheck `pwd`/src/runZGammaJetCrossCheck.cc $(CFLAGS) $(LDFLAGS) -lZGammaJetCrossCheck

runZGammaJetCrossCheckHist :
	$(CXX) -o `pwd`/lib/runZGammaJetCrossCheckHist `pwd`/src/runZGammaJetCrossCheckHist.cc $(CFLAGS) $(LDFLAGS) -lZGammaJetCrossCheckHist

clean :
	rm -rf ./lib/*
