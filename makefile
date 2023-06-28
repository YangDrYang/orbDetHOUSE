CXX = g++
CC = gcc

# create directories
$(shell mkdir -p bin/orbmdl)
$(shell mkdir -p bin/orbmdl/3rdparty)
$(shell mkdir -p bin/orbmdl/sofa)
$(shell mkdir -p bin/orbmdl/nrlmsise-00)
$(shell mkdir -p bin/filter)
$(shell mkdir -p bin/scripts)

OBJDIR = bin
LocalDIR = usr/local
FUNDIR = scripts

# for compiling orbmdl and filter
CPPFLAGS = --std=c++11 -Wall -pedantic -g -O3 
CFLAGS = -Wall -pedantic -g -O3 
INCLUDE = -I/$(LocalDIR)/include/ -I/$(LocalDIR)/include/eigen3/ -I/$(LocalDIR)/boost_1_81_0 -I./orbmdl -I./orbmdl/3rdparty -I./orbmdl/sofa -I./orbmdl/nrlmsise-00 -I./filter

# orbmdl files
ORBMDL_SRC = $(wildcard orbmdl/*.cpp) $(wildcard orbmdl/3rdparty/*.cpp) $(wildcard orbmdl/sofa/*.cpp)
ORBMDL_OBJ = $(ORBMDL_SRC:orbmdl/%.cpp=$(OBJDIR)/orbmdl/%.o)
ORBMDL_AR = $(OBJDIR)/orbmdl/liborbmdl.a

# nrlmsise-00 (drag density model) files
NRLMSISE00_SRC = orbmdl/nrlmsise-00/nrlmsise-00.c orbmdl/nrlmsise-00/nrlmsise-00_data.c
NRLMSISE00_OBJ = $(NRLMSISE00_SRC:orbmdl/nrlmsise-00/%.c=$(OBJDIR)/orbmdl/nrlmsise-00/%.o)
NRLMSISE00_AR = $(OBJDIR)/orbmdl/libnrlmsise00.a

# filter files
FILTER_SRC = $(wildcard filter/*.cpp)
FILTER_OBJ = $(FILTER_SRC:filter/%.cpp=$(OBJDIR)/filter/%.o)
FILTER_AR = $(OBJDIR)/filter/libfilter.a

# ORBDET files
FILENAME ?= 
ORBDET_SRC = $(FUNDIR)/$(FILENAME).cpp
ORBDET_OBJ = $(ORBDET_SRC:$(FUNDIR)/%.cpp=$(OBJDIR)/$(FUNDIR)/%.o)

OBJECTS = $(ORBMDL_AR) $(FILTER_AR) $(NRLMSISE00_AR)

# ----- COMPILE ORBDET -----
$(OBJDIR)/$(FUNDIR)/$(FILENAME): $(ORBDET_OBJ) $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(OBJECTS) $(ORBDET_OBJ) -L$(OBJDIR)/filter -lfilter -lyaml-cpp -L$(OBJDIR)/orbmdl -lorbmdl -lnrlmsise00 -o $@

# compile ORBDET object file
$(ORBDET_OBJ): $(OBJDIR)/$(FUNDIR)/%.o : $(FUNDIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

# ----- COMPILE orbmdl -----
# compile orbmdl object files
$(ORBMDL_OBJ): $(OBJDIR)/orbmdl/%.o : orbmdl/%.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

# ----- COMPILE nrlmsise-00 -----
# compile nrlmsise-00 object files
$(NRLMSISE00_OBJ): $(OBJDIR)/orbmdl/nrlmsise-00/%.o : orbmdl/nrlmsise-00/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# create orbmdl library
$(ORBMDL_AR): $(ORBMDL_OBJ)
	ar rcs $(ORBMDL_AR) $(ORBMDL_OBJ)

# create nrlmsise-00 library
$(NRLMSISE00_AR): $(NRLMSISE00_OBJ)
	ar rcs $(NRLMSISE00_AR) $(NRLMSISE00_OBJ)

# ----- COMPILE filter -----
# create filter library
$(FILTER_AR): $(FILTER_OBJ)
	ar rcs $(FILTER_AR) $(FILTER_OBJ)

# compile filter object files
$(FILTER_OBJ): $(OBJDIR)/filter/%.o : filter/%.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

# ----- CLEAN -----
.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/orbmdl/*.o $(OBJDIR)/orbmdl/3rdparty/*.o $(OBJDIR)/orbmdl/sofa/*.o $(OBJDIR)/orbmdl/nrlmsise-00/*.o $(ORBMDL_AR) $(NRLMSISE00_AR) $(OBJDIR)/filter/*.o $(FILTER_AR) $(OBJDIR)/$(FUNDIR)/$(FILENAME)
