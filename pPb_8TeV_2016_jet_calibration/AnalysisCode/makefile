CXXFLAGS=-shared -Ofast -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I$(CALIB_PATH)/AnalysisCode/include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L$(CALIB_PATH)/AnalysisCode/lib -lGlobalParams -lInitialization -lTrigger -lAtlas -ldl

all : Params TreeVariables RooZfit EMTopoComparison EMTopoComparisonHist ElectronContaminationStudy ElectronContaminationStudyHist EnergyScaleChecks EnergyScaleChecksHist JetInsituCorrectionCheck JetInsituCorrectionCheckHist RtrkComparison RtrkComparisonHist ZGammaJetCrossCheck ZGammaJetCrossCheckHist

Params : src/Params.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libParams.so $(CALIB_PATH)/AnalysisCode/src/Params.C

TreeVariables : src/TreeVariables.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libTreeVariables.so $(CALIB_PATH)/AnalysisCode/src/TreeVariables.C

RooZfit : src/RooZfit.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lRooFit -lRooFitCore -o $(CALIB_PATH)/AnalysisCode/lib/libRooZfit.so $(CALIB_PATH)/AnalysisCode/src/RooZfit.C

EMTopoComparison : src/EMTopoComparison.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libEMTopoComparison.so $(CALIB_PATH)/AnalysisCode/src/EMTopoComparison.C

EMTopoComparisonHist : src/EMTopoComparisonHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams  -o $(CALIB_PATH)/AnalysisCode/lib/libEMTopoComparisonHist.so $(CALIB_PATH)/AnalysisCode/src/EMTopoComparisonHist.C

ElectronContaminationStudy : src/ElectronContaminationStudy.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libElectronContaminationStudy.so $(CALIB_PATH)/AnalysisCode/src/ElectronContaminationStudy.C

ElectronContaminationStudyHist : src/ElectronContaminationStudyHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -o $(CALIB_PATH)/AnalysisCode/lib/libElectronContaminationStudyHist.so $(CALIB_PATH)/AnalysisCode/src/ElectronContaminationStudyHist.C

EnergyScaleChecks : src/EnergyScaleChecks.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libEnergyScaleChecks.so $(CALIB_PATH)/AnalysisCode/src/EnergyScaleChecks.C

EnergyScaleChecksHist : src/EnergyScaleChecksHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -o $(CALIB_PATH)/AnalysisCode/lib/libEnergyScaleChecksHist.so $(CALIB_PATH)/AnalysisCode/src/EnergyScaleChecksHist.C

JetInsituCorrectionCheck : src/JetInsituCorrectionCheck.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libJetInsituCorrectionCheck.so $(CALIB_PATH)/AnalysisCode/src/JetInsituCorrectionCheck.C

JetInsituCorrectionCheckHist : src/JetInsituCorrectionCheckHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -o $(CALIB_PATH)/AnalysisCode/lib/libJetInsituCorrectionCheckHist.so $(CALIB_PATH)/AnalysisCode/src/JetInsituCorrectionCheckHist.C

RtrkComparison : src/RtrkComparison.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libRtrkComparison.so $(CALIB_PATH)/AnalysisCode/src/RtrkComparison.C

RtrkComparisonHist : src/RtrkComparisonHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -o $(CALIB_PATH)/AnalysisCode/lib/libRtrkComparisonHist.so $(CALIB_PATH)/AnalysisCode/src/RtrkComparisonHist.C

ZGammaJetCrossCheck : src/ZGammaJetCrossCheck.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libZGammaJetCrossCheck.so $(CALIB_PATH)/AnalysisCode/src/ZGammaJetCrossCheck.C

ZGammaJetCrossCheckHist : src/ZGammaJetCrossCheckHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -o $(CALIB_PATH)/AnalysisCode/lib/libZGammaJetCrossCheckHist.so $(CALIB_PATH)/AnalysisCode/src/ZGammaJetCrossCheckHist.C

clean :
	rm -rf ./lib/*.so
