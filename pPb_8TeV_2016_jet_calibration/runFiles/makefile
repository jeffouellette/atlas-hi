CFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I$(CALIB_PATH)/AnalysisCode/include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L$(CALIB_PATH)/AnalysisCode/lib -lAtlas -lGlobalParams -lInitialization -lTrigger -ldl


all : runEMTopoComparison runEMTopoComparisonHist runElectronContaminationStudy runElectronContaminationStudyHist runEnergyScaleChecks runEnergyScaleChecksHist runFCalDistribution runFCalDistributionHist runJetInsituCorrectionCheck runJetInsituCorrectionCheckHist runRtrkComparison runRtrkComparisonHist runZGammaJetCrossCheck runZGammaJetCrossCheckHist runZMassCalc runZMassCalcHist

runEMTopoComparison :
	$(CXX) -o lib/runEMTopoComparison src/runEMTopoComparison.cc $(CFLAGS) $(LDFLAGS) -lEMTopoComparison

runEMTopoComparisonHist :
	$(CXX) -o lib/runEMTopoComparisonHist src/runEMTopoComparisonHist.cc $(CFLAGS) $(LDFLAGS) -lEMTopoComparisonHist 

runElectronContaminationStudy :
	$(CXX) -o lib/runElectronContaminationStudy src/runElectronContaminationStudy.cc $(CFLAGS) $(LDFLAGS) -lElectronContaminationStudy

runElectronContaminationStudyHist :
	$(CXX) -o lib/runElectronContaminationStudyHist src/runElectronContaminationStudyHist.cc $(CFLAGS) $(LDFLAGS) -lElectronContaminationStudyHist

runEnergyScaleChecks :
	$(CXX) -o lib/runEnergyScaleChecks src/runEnergyScaleChecks.cc $(CFLAGS) $(LDFLAGS) -lEnergyScaleChecks

runEnergyScaleChecksHist :
	$(CXX) -o lib/runEnergyScaleChecksHist src/runEnergyScaleChecksHist.cc $(CFLAGS) $(LDFLAGS) -lEnergyScaleChecksHist

runFCalDistribution :
	$(CXX) -o lib/runFCalDistribution src/runFCalDistribution.cc $(CFLAGS) $(LDFLAGS) -lFCalDistribution

runFCalDistributionHist :
	$(CXX) -o lib/runFCalDistributionHist src/runFCalDistributionHist.cc $(CFLAGS) $(LDFLAGS) -lFCalDistributionHist

runJetInsituCorrectionCheck :
	$(CXX) -o lib/runJetInsituCorrectionCheck src/runJetInsituCorrectionCheck.cc $(CFLAGS) $(LDFLAGS) -lJetInsituCorrectionCheck

runJetInsituCorrectionCheckHist :
	$(CXX) -o lib/runJetInsituCorrectionCheckHist src/runJetInsituCorrectionCheckHist.cc $(CFLAGS) $(LDFLAGS) -lJetInsituCorrectionCheckHist

runRtrkComparison :
	$(CXX) -o lib/runRtrkComparison src/runRtrkComparison.cc $(CFLAGS) $(LDFLAGS) -lRtrkComparison

runRtrkComparisonHist :
	$(CXX) -o lib/runRtrkComparisonHist src/runRtrkComparisonHist.cc $(CFLAGS) $(LDFLAGS) -lRtrkComparisonHist

runZGammaJetCrossCheck :
	$(CXX) -o lib/runZGammaJetCrossCheck src/runZGammaJetCrossCheck.cc $(CFLAGS) $(LDFLAGS) -lZGammaJetCrossCheck

runZGammaJetCrossCheckHist :
	$(CXX) -o lib/runZGammaJetCrossCheckHist src/runZGammaJetCrossCheckHist.cc $(CFLAGS) $(LDFLAGS) -lZGammaJetCrossCheckHist

runZMassCalc :
	$(CXX) -o lib/runZMassCalc src/runZMassCalc.cc $(CFLAGS) $(LDFLAGS) -lZMassCalc

runZMassCalcHist :
	$(CXX) -o lib/runZMassCalcHist src/runZMassCalcHist.cc $(CFLAGS) $(LDFLAGS) -lZMassCalcHist

clean :
	rm -rf ./lib/*
