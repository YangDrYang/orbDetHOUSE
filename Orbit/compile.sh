echo "Compiling orbit..."

g++ -o orbit.exe orbit.cpp -I.. -I./../eigen ../house.a \
             -O3 -g -Wall -pedantic

