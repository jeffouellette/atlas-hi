CXXFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I$(JETTRACKS_PATH)/AnalysisCode/include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L$(JETTRACKS_PATH)/AnalysisCode/lib -ldl

analyses = runJetTrackEtaPhi runRtrkTreeMaker

all : $(analyses)

run%: obj/run%.o
	$(CXX) $(LDFLAGS) -l$(@:run%=%) obj/$@.o -o bin/$@

obj/run%.o:
	$(CXX) $(CXXFLAGS) -c $(@:obj/run%.o=src/run%.cc) -o $@

clean :
	rm -rf ./bin/*
