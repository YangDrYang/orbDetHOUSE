# compiler
CC = g++

CPPFLAGS = -g -Wall -I./eigen -I. -B

SUPPORTING_FILES = dyn.cpp filter_aux.cpp house.cpp ode.cpp pearsonator.cpp timer.cpp ukf.cpp 

all: projectile rigid_body

projectile: ./projectile/projectile.cpp $(SUPPORTING_FILES)
	$(CC) $(CPPFLAGS) -o projectile ./projectile/projectile.cpp $(SUPPORTING_FILES)

rigid_body: ./rigid_body/rigid_body.cpp $(SUPPORTING_FILES)
	$(CC) $(CPPFLAGS) -o rigid_body ./rigid_body/rigid_body.cpp $(SUPPORTING_FILES)


