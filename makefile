CXXFLAGS=-shared -g -Wall `root-config --cflags` -I$(ATLAS_PATH)/include -fPIC
LDFLAGS=`root-config --glibs` -L$(ATLAS_PATH)/lib

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : GlobalParams Utilities Trigger Initialization

GlobalParams :
	$(CC) -o $(ATLAS_PATH)/lib/libGlobalParams.so $(ATLAS_PATH)/src/GlobalParams.cxx

Utilities :
	$(CC) -lGlobalParams -o $(ATLAS_PATH)/lib/libUtilities.so $(ATLAS_PATH)/src/Utilities.cxx

Trigger :
	$(CC) -o $(ATLAS_PATH)/lib/libTrigger.so $(ATLAS_PATH)/src/Trigger.cxx

Initialization :
	$(CC) -lTrigger -lGlobalParams -o $(ATLAS_PATH)/lib/libInitialization.so $(ATLAS_PATH)/src/Initialization.cxx

clean :
	rm -rf lib/*.so*
