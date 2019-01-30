CXXFLAGS=-Ofast -fPIC -g -Wall `root-config --cflags` -I$(HOME)/RootUtils -I$(ATLAS_PATH)/include -I./include
LDFLAGS=`root-config --glibs` -L$(HOME)/RootUtils -L$(ATLAS_PATH)/lib -L./lib -lUtilities -lGlobalParams -lTrigger -lAtlas -ldl

CC=$(CXX) $(CXXFLAGS) $(LDFLAGS)

analyses = JetTrackEtaPhi

all : $(analyses:%=lib/lib%.so) $(analyses:%=lib/lib%Hist.so)

TreeVariables: src/TreeVariables.cxx
	$(CC) -shared -o lib/libTreeVariables.so src/TreeVariables.cxx

# Libraries need to be compiled into object files and then migrated to a shared library
lib/lib%.so: obj/%.o
	$(CXX) $(LDFLAGS) -lUtilities -lTreeVariables -shared $< -o $@

obj/%.o :	src/%.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

% : obj/%.o
	$(CXX) $(LDFLAGS) -lUtilities -lTreeVariables -shared $< -o lib/lib$@.so

clean :
	rm -rf ./lib/*.so*
