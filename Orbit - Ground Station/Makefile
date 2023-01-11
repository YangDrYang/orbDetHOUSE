CC=g++
CFLAGS = -w -I.. -I./../eigen -I./../orbitDetermination -I./../orbitDetermination/sofa -I./../orbitDetermination/3rdparty  -I./../yaml-cpp/include -O3 -fmax-errors=5

3rdparty = $(wildcard ../orbitDetermination/3rdparty/*.cpp)
sofa = $(wildcard ../orbitDetermination/sofa/*.cpp)

orbit: constants.o coordTrans.o satRefSys.o jplEph.o gravity.o forceModels.o auxillaryData.o
	$(CC) $(CFLAGS) $(3rdparty) $(sofa) $(wildcard *.o) orbit.cpp -o orbit.exe ../house.a -g

validateModel: constants.o coordTrans.o satRefSys.o jplEph.o gravity.o forceModels.o auxillaryData.o
	$(CC) $(CFLAGS) $(3rdparty) $(sofa) $(wildcard *.o) validateModel.cpp -o validateModel.exe ../house.a -g

constants.o: ../orbitDetermination/constants.cpp
	$(CC) $(CFLAGS) -c ../orbitDetermination/constants.cpp

coordTrans.o: ../orbitDetermination/coordTrans.cpp
	$(CC) $(CFLAGS) -c ../orbitDetermination/coordTrans.cpp

satRefSys.o: ../orbitDetermination/satRefSys.cpp
	$(CC) $(CFLAGS) -c ../orbitDetermination/satRefSys.cpp

jplEph.o: ../orbitDetermination/jplEph.cpp
	$(CC) $(CFLAGS) -c ../orbitDetermination/jplEph.cpp

gravity.o: ../orbitDetermination/gravity.cpp
	$(CC) $(CFLAGS)  -c ../orbitDetermination/gravity.cpp

forceModels.o: ../orbitDetermination/forceModels.cpp
	$(CC) $(CFLAGS) -c ../orbitDetermination/forceModels.cpp

auxillaryData.o: ../orbitDetermination/auxillaryData.cpp
	$(CC) $(CFLAGS) -c ../orbitDetermination/auxillaryData.cpp

clean:
	rm -f *.o *.exe