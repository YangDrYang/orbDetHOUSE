CC=g++
CFLAGS = -w -I.. -I./../eigen -I./../orbitDetermination -I./../orbitDetermination/sofa -I./../orbitDetermination/3rdparty  -I./../yaml-cpp/include -O3 -fmax-errors=5

# source is orbitDetermination folder
SRCDIR = ./../orbitDetermination
OBJDIR = ./bin
SRC = constants.o coordTrans.o satRefSys.o jplEph.o gravity.o forceModels.o auxillaryData.o

objects = $(addprefix $(OBJDIR)/, $(SRC))

3rdparty = $(wildcard ../orbitDetermination/3rdparty/*.cpp)
sofa = $(wildcard ../orbitDetermination/sofa/*.cpp)

# link object files
orbit: constants.o coordTrans.o satRefSys.o jplEph.o gravity.o forceModels.o auxillaryData.o
	$(CC) $(CFLAGS) $(3rdparty) $(sofa) $(wildcard *.o) orbit.cpp -o orbit.exe ../house.a -g

validateModel: $(objects)
	$(CC) $(CFLAGS) $(3rdparty) $(sofa) $(objects) validateModel.cpp -o validateModel.exe ../house.a -g


# yaml library flag must go at end (for some reason)
test: _testYaml.cpp
	$(CC) -I./../eigen  _testYaml.cpp -o _testYaml.exe -lyaml-cpp 

# compiles to object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJDIR)/*.o *.exe