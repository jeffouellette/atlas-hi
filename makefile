CXX=clang++
CXXFLAGS=-shared -g -Wall `root-config --cflags` -I$(ATLAS_PATH)/include
LDFLAGS=`root-config --glibs` -L$(ATLAS_PATH)/lib

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : GlobalParams Utils Trigger Initialization TreeVariables

GlobalParams :
	$(CC) -o $(ATLAS_PATH)/lib/libGlobalParams.so $(ATLAS_PATH)/src/GlobalParams.C

Utils :
	$(CC) -lGlobalParams -o $(ATLAS_PATH)/lib/libUtils.so $(ATLAS_PATH)/src/Utils.C

Trigger :
	$(CC) -o $(ATLAS_PATH)/lib/libTrigger.so $(ATLAS_PATH)/src/Trigger.C

Initialization :
	$(CC) -lTrigger -lGlobalParams -o $(ATLAS_PATH)/lib/libInitialization.so $(ATLAS_PATH)/src/Initialization.C

TreeVariables :
	$(CC) -o $(ATLAS_PATH)/lib/libTreeVariables.so $(ATLAS_PATH)/src/TreeVariables.C

clean :
	rm -rf lib/*.so*
