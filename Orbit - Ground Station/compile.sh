echo "Compiling orbit..."

# g++ -o orbit.exe orbit.cpp ./../orbitDetermination/auxillaryData.cpp ./../orbitDetermination/satRefSys.cpp \
#     ./../orbitDetermination/satRefSys.hpp ./../orbitDetermination/enum.h ./../orbitDetermination/coordTrans.cpp \
#     ./../orbitDetermination/config.hpp ./../orbitDetermination/config.cpp ./../orbitDetermination/sofa/*.cpp \
#     ./../orbitDetermination/gravity.* ./../orbitDetermination/gTime.* \
#     ./../orbitDetermination/forceModels.* ./../orbitDetermination/jplEph.* ./../orbitDetermination/3rdparty/*.cpp \
#     -I.. -I./../eigen -I./../orbitDetermination -I./../orbitDetermination/sofa -I./../orbitDetermination/3rdparty ../house.a  \
#     -O3 -g -fmax-errors=5

g++ -o orbit.exe orbit.cpp ./../orbitDetermination/*.cpp ./../orbitDetermination/3rdparty/*.cpp ./../orbitDetermination/sofa/*.cpp \
    -I.. -I./../eigen -I./../orbitDetermination -I./../orbitDetermination/sofa -I./../orbitDetermination/3rdparty ../house.a  \
    -O3 -g -fmax-errors=5