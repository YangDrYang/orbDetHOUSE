CC=g++
CFLAGS = -I.. -I./../eigen -I./../orbitDetermination ../house.a \ -O3 -g -Wall -pedantic -fmax-errors=5


orbit: orbit.o
    $(CC) $(CFLAGS) orbit.exe 


orbit.o: $(CC) $(CFLAGS) -c orbit.cpp
