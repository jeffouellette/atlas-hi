CXXFLAGS=-O3 -shared -g -Wall `root-config --cflags` -I${ATLAS_PATH}/include -I${ROOT_UTILS_PATH}/include -fPIC
LDFLAGS=`root-config --glibs` -L${ATLAS_PATH}/lib -L${ROOT_UTILS_PATH}/lib -lUtilities

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : GlobalParams AtlasUtils AtlasStyle

GlobalParams :
	$(CC) -o lib/libGlobalParams.so src/GlobalParams.cxx

#Utilities :
#	$(CC) -lGlobalParams -o ${ATLAS_PATH}/lib/libUtilities.so ${ATLAS_PATH}/src/Utilities.cxx

#Trigger :
#	$(CC) -o lib/libTrigger.so src/Trigger.cxx

#Initialization :
#	$(CC) -lTrigger -lGlobalParams -o lib/libInitialization.so src/Initialization.cxx

AtlasUtils : src/AtlasUtils.cxx AtlasStyle
	$(CC) -lAtlasStyle -o lib/libAtlasUtils.so src/AtlasUtils.cxx

AtlasStyle : src/AtlasStyle.cxx
	$(CC) -o lib/libAtlasStyle.so src/AtlasStyle.cxx


clean :
	rm -rf lib/*.so*
