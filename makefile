# compiler
CC = g++

OBJDIR = bin
LOCAL_DIR = usr/local
FILTER_DIR = filters
ORBDET_DIR = orbitDetermination


# ----- Don't modify below this point -----
# for compiling house.a
CPPFLAGS1 = --std=c++11 -Wall -pedantic -g -O3 
INCLUDE1 = -I/$(LOCAL_DIR)/include/eigen3/ -I/$(LOCAL_DIR)/boost_1_81_0 -I./$(FILTER_DIR)

# for compiling orbit determination
CPPFLAGS2 =  --std=c++11 -Wall -pedantic -g -O3 
INCLUDE2 = -w -I. -I/$(LOCAL_DIR)/include/eigen3/ -I/$(LOCAL_DIR)/boost_1_81_0 -I./$(FILTER_DIR) -I./$(ORBDET_DIR) -I./$(ORBDET_DIR)/sofa \
-I./$(ORBDET_DIR)/3rdparty -O3 -fmax-errors=5 -I./nrlmsise-00

3rdparty = $(wildcard $(ORBDET_DIR)/3rdparty/*.cpp)
sofa = $(wildcard $(ORBDET_DIR)/sofa/*.cpp)


FILTER_SRC = $(wildcard $(FILTER_DIR)/*.cpp)
FILTER_OBJ = $(FILTER_SRC:$(FILTER_DIR)/%.cpp=$(OBJDIR)/%.o)

# only required files, so manually including them
ORBDET_SRC_FILES = constants.cpp coordTrans.cpp satRefSys.cpp jplEph.cpp gravity.cpp forceModels.cpp auxillaryData.cpp
ORBDET_SRC = $(addprefix $(ORBDET_DIR)/,$(ORBDET_SRC_FILES))
ORBDET_OBJ = $(ORBDET_SRC:$(ORBDET_DIR)/%.cpp=$(OBJDIR)/%.o)

OBJECTS = $(ORBDET_OBJ) $(FILTER_OBJ) $(OBJDIR)/house.a $(OBJDIR)/nrlmsise-00.o $(OBJDIR)/nrlmsise-00_data.o
# ---------------------------------------
# ----- COMPILE ORBIT DETERMINTAION -----
# ---------------------------------------
filter_testing: filter_testing.cpp filter_testing.hpp $(ORBDET_OBJ) $(OBJDIR)/house.a \
	$(OBJDIR)/nrlmsise-00.o $(OBJDIR)/nrlmsise-00_data.o 
	$(CC) $(CPPFLAGS2) $(INCLUDE2) $(OBJECTS) $(3rdparty) $(sofa) filter_testing.cpp -o filter_testing -lyaml-cpp 
	
# compile orbit determination object files
$(ORBDET_OBJ): $(OBJDIR)/%.o : $(ORBDET_DIR)/%.cpp
	$(CC) $(CPPFLAGS2) $(INCLUDE2) -c $< -o $@ 

# --------------------------------------------------	
# ----- COMPILE NRLMSISE-00 (atmosphere model) -----
# --------------------------------------------------

$(OBJDIR)/nrlmsise-00.o: nrlmsise-00/nrlmsise-00.c nrlmsise-00/nrlmsise-00.h
	gcc nrlmsise-00/nrlmsise-00.c -c -o $(OBJDIR)/nrlmsise-00.o

$(OBJDIR)/nrlmsise-00_data.o: nrlmsise-00/nrlmsise-00_data.c nrlmsise-00/nrlmsise-00.h
	gcc nrlmsise-00/nrlmsise-00_data.c -c -o $(OBJDIR)/nrlmsise-00_data.o


# ---------------------------	
# ----- COMPILE FILTERS -----
# ---------------------------	

# linking filters
$(OBJDIR)/house.a: $(FILTER_OBJ)
	ar rcs $(OBJDIR)/house.a $(FILTER_OBJ)

# compile filter object files
$(FILTER_OBJ): $(OBJDIR)/%.o : $(FILTER_DIR)/%.cpp
	$(CC) $(CPPFLAGS1) $(INCLUDE1) -c $< -o $@ 



# ----- CLEAN -----
.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.a *.exe