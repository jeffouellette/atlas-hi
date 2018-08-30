CXX=clang++
CXXFLAGS=-shared -g -Wall `root-config --cflags` -I$(ATLAS_PATH)/include
LDFLAGS=`root-config --glibs` -L$(ATLAS_PATH)/lib

all : GlobalParams Trigger Initialization

GlobalParams :
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(ATLAS_PATH)/lib/libGlobalParams.so $(ATLAS_PATH)/src/GlobalParams.C

Trigger :
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(ATLAS_PATH)/lib/libTrigger.so $(ATLAS_PATH)/src/Trigger.C

Initialization :
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lTrigger -lGlobalParams -o $(ATLAS_PATH)/lib/libInitialization.so $(ATLAS_PATH)/src/Initialization.C

clean :
	rm -rf ./lib/*.so*
