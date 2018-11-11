CXXFLAGS=-Ofast -fPIC -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I./include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L./lib -lUtils -lGlobalParams -lTrigger -lTreeVariables -lAtlas -ldl

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

analyses = PhotonAnalysis

all : $(analyses:%=lib/lib%.so) $(analyses:%=lib/lib%Hist.so)

# Libraries need to be compiled into object files and then migrated to a shared library
lib/lib%.so: obj/%.o
	$(CXX) $(LDFLAGS) -lUtils -shared $< -o $@

obj/%.o :	src/%.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

% : obj/%.o
	$(CXX) $(LDFLAGS) -lUtils -shared $< -o lib/lib$@.so

clean :
	rm -rf ./lib/*.so*
