CXXFLAGS=-shared -Ofast -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I$(DIJET_PATH)/AnalysisCode/include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L$(DIJET_PATH)/AnalysisCode/lib -lGlobalParams -lInitialization -lTrigger -lAtlas -ldl

all : Params TreeVariables Util DijetAnalysis DijetAnalysisHist InclusiveJetAnalysis InclusiveJetAnalysisHist InclusiveRpPbAnalysis InclusiveRpPbAnalysisHist TriggerEfficiencyAnalysis TriggerEfficiencyAnalysisHist TriggerPtAnalysis TriggerPtAnalysisCombinedHist TriggerPtAnalysisHist

Params : src/Params.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(DIJET_PATH)/AnalysisCode/lib/libParams.so $(DIJET_PATH)/AnalysisCode/src/Params.C

TreeVariables : src/TreeVariables.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(DIJET_PATH)/AnalysisCode/lib/libTreeVariables.so $(DIJET_PATH)/AnalysisCode/src/TreeVariables.C

Util : src/Util.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -o $(DIJET_PATH)/AnalysisCode/lib/libUtil.so $(DIJET_PATH)/AnalysisCode/src/Util.C

DijetAnalysis : src/DijetAnalysis.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -lTreeVariables -o $(DIJET_PATH)/AnalysisCode/lib/libDijetAnalysis.so $(DIJET_PATH)/AnalysisCode/src/DijetAnalysis.C

DijetAnalysisHist : src/DijetAnalysisHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtils -o $(DIJET_PATH)/AnalysisCode/lib/libDijetAnalysisHist.so $(DIJET_PATH)/AnalysisCode/src/DijetAnalysisHist.C

InclusiveJetAnalysis : src/InclusiveJetAnalysis.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -lTreeVariables -o $(DIJET_PATH)/AnalysisCode/lib/libInclusiveJetAnalysis.so $(DIJET_PATH)/AnalysisCode/src/InclusiveJetAnalysis.C

InclusiveJetAnalysisHist : src/InclusiveJetAnalysisHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -o $(DIJET_PATH)/AnalysisCode/lib/libInclusiveJetAnalysisHist.so $(DIJET_PATH)/AnalysisCode/src/InclusiveJetAnalysisHist.C

InclusiveRpPbAnalysis : src/InclusiveRpPbAnalysis.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -lTreeVariables -o $(DIJET_PATH)/AnalysisCode/lib/libInclusiveRpPbAnalysis.so $(DIJET_PATH)/AnalysisCode/src/InclusiveRpPbAnalysis.C

InclusiveRpPbAnalysisHist : src/InclusiveRpPbAnalysisHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -o $(DIJET_PATH)/AnalysisCode/lib/libInclusiveRpPbAnalysisHist.so $(DIJET_PATH)/AnalysisCode/src/InclusiveRpPbAnalysisHist.C

TriggerEfficiencyAnalysis : src/TriggerEfficiencyAnalysis.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -lTreeVariables -o $(DIJET_PATH)/AnalysisCode/lib/libTriggerEfficiencyAnalysis.so $(DIJET_PATH)/AnalysisCode/src/TriggerEfficiencyAnalysis.C

TriggerEfficiencyAnalysisHist : src/TriggerEfficiencyAnalysisHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -o $(DIJET_PATH)/AnalysisCode/lib/libTriggerEfficiencyAnalysisHist.so $(DIJET_PATH)/AnalysisCode/src/TriggerEfficiencyAnalysisHist.C

TriggerPtAnalysis : src/RtrkComparison.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -lTreeVariables -o $(DIJET_PATH)/AnalysisCode/lib/libTriggerPtAnalysis.so $(DIJET_PATH)/AnalysisCode/src/RtrkComparison.C

TriggerPtAnalysisHist : src/RtrkComparisonHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -o $(DIJET_PATH)/AnalysisCode/lib/libTriggerPtAnalysisHist.so $(DIJET_PATH)/AnalysisCode/src/RtrkComparisonHist.C

TriggerPtAnalysisCombinedHist : src/TriggerPtAnalysisCombinedHist.C
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lParams -lUtil -lTreeVariables -o $(DIJET_PATH)/AnalysisCode/lib/libTriggerPtAnalysisCombinedHist.so $(DIJET_PATH)/AnalysisCode/src/TriggerPtAnalysisCombinedHist.C

clean :
	rm -rf ./lib/*.so
