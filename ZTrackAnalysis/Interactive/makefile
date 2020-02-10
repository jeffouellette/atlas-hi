CXXFLAGS=-O3 -g -Wall `root-config --cflags` -I${ATLAS_PATH}/include -I${HOME}/RootUtils -I${ATLAS_PATH}/ZTrackAnalysis/Interactive
LDFLAGS=`root-config --glibs` -L${ATLAS_PATH}/lib -L${HOME}/RootUtils -lGlobalParams -lUtilities -lAtlasUtils -lAtlasStyle

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

all : ProcessFile

ProcessFile : ProcessFile.cxx *.h
	$(CC) -o ${ATLAS_PATH}/ZTrackAnalysis/Interactive/ProcessFile ${ATLAS_PATH}/ZTrackAnalysis/Interactive/ProcessFile.cxx

clean :
	rm ${ATLAS_PATH}/ZTrackAnalysis/Interactive/ProcessFile
