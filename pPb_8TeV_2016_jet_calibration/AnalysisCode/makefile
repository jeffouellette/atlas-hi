CFLAGS=-shared -Ofast -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I$(CALIB_PATH)/AnalysisCode/include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L$(CALIB_PATH)/AnalysisCode/lib -lGlobalParams -lInitialization -lTrigger -lAtlas

all : TreeVariables RooZfit EMTopoComparison EMTopoComparisonHist ElectronContaminationStudy ElectronContaminationStudyHist EnergyScaleChecks EnergyScaleChecksHist JetInsituCorrectionCheck JetInsituCorrectionCheckHist RtrkComparison RtrkComparisonHist ZGammaJetCrossCheck ZGammaJetCrossCheckHist

TreeVariables :
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libTreeVariables.so $(CALIB_PATH)/AnalysisCode/src/TreeVariables.C

RooZfit :
	$(CXX) $(CFLAGS) $(LDFLAGS) -lRooFit -lRooFitCore -o $(CALIB_PATH)/AnalysisCode/lib/libRooZfit.so $(CALIB_PATH)/AnalysisCode/src/RooZfit.C

EMTopoComparison :
	$(CXX) $(CFLAGS) $(LDFLAGS) -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libEMTopoComparison.so $(CALIB_PATH)/AnalysisCode/src/EMTopoComparison.C

EMTopoComparisonHist :
	$(CXX) $(CFLAGS) $(LDFLAGS)  -o $(CALIB_PATH)/AnalysisCode/lib/libEMTopoComparisonHist.so $(CALIB_PATH)/AnalysisCode/src/EMTopoComparisonHist.C

ElectronContaminationStudy :
	$(CXX) $(CFLAGS) $(LDFLAGS) -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libElectronContaminationStudy.so $(CALIB_PATH)/AnalysisCode/src/ElectronContaminationStudy.C

ElectronContaminationStudyHist :
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libElectronContaminationStudyHist.so $(CALIB_PATH)/AnalysisCode/src/ElectronContaminationStudyHist.C

EnergyScaleChecks :
	$(CXX) $(CFLAGS) $(LDFLAGS) -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libEnergyScaleChecks.so $(CALIB_PATH)/AnalysisCode/src/EnergyScaleChecks.C

EnergyScaleChecksHist :
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libEnergyScaleChecksHist.so $(CALIB_PATH)/AnalysisCode/src/EnergyScaleChecksHist.C

JetInsituCorrectionCheck :
	$(CXX) $(CFLAGS) $(LDFLAGS) -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libJetInsituCorrectionCheck.so $(CALIB_PATH)/AnalysisCode/src/JetInsituCorrectionCheck.C

JetInsituCorrectionCheckHist :
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libJetInsituCorrectionCheckHist.so $(CALIB_PATH)/AnalysisCode/src/JetInsituCorrectionCheckHist.C

RtrkComparison :
	$(CXX) $(CFLAGS) $(LDFLAGS) -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libRtrkComparison.so $(CALIB_PATH)/AnalysisCode/src/RtrkComparison.C

RtrkComparisonHist :
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libRtrkComparisonHist.so $(CALIB_PATH)/AnalysisCode/src/RtrkComparisonHist.C

ZGammaJetCrossCheck :
	$(CXX) $(CFLAGS) $(LDFLAGS) -lTreeVariables -o $(CALIB_PATH)/AnalysisCode/lib/libZGammaJetCrossCheck.so $(CALIB_PATH)/AnalysisCode/src/ZGammaJetCrossCheck.C

ZGammaJetCrossCheckHist :
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(CALIB_PATH)/AnalysisCode/lib/libZGammaJetCrossCheckHist.so $(CALIB_PATH)/AnalysisCode/src/ZGammaJetCrossCheckHist.C

clean :
	rm -rf ./lib/*.so
