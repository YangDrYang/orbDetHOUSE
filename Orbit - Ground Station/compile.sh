echo "Compiling orbit..."

g++ -o orbit.exe orbit.cpp ./../orbitDetermination/auxillaryData.cpp ./../orbitDetermination/satRefSys.cpp \
    ./../orbitDetermination/satRefSys.hpp ./../orbitDetermination/enum.h ./../orbitDetermination/coordTrans.cpp\
    ./../orbitDetermination/config.hpp ./../orbitDetermination/config.cpp \
    -I.. -I./../eigen -I./../orbitDetermination -I./../orbitDetermination/sofa ../house.a \
    -O3 -g -Wall -pedantic -fmax-errors=5

