# Compiler settings
CXX = g++-11
CC = gcc-11

# create directories
$(shell mkdir -p bin/orbmdl)
$(shell mkdir -p bin/orbmdl/3rdparty)
$(shell mkdir -p bin/orbmdl/sofa)
$(shell mkdir -p bin/orbmdl/nrlmsise-00)
$(shell mkdir -p bin/filter)
$(shell mkdir -p bin/scripts)
$(shell mkdir -p lib)

# Directories
OBJDIR = bin
LIBDIR = lib
LocalDIR = /usr/local
FUNDIR = scripts

# Source files for building object files and libraries
ORBMDL_SRC = $(wildcard orbmdl/*.cpp) $(wildcard orbmdl/3rdparty/*.cpp) $(wildcard orbmdl/sofa/*.cpp)
NRLMSISE00_SRC = orbmdl/nrlmsise-00/nrlmsise-00.c orbmdl/nrlmsise-00/nrlmsise-00_data.c
FILTER_SRC = $(wildcard filter/*.cpp)
SCRIPT_SRC = $(wildcard scripts/*.cpp)

# Object files for libraries
ORBMDL_OBJ = $(ORBMDL_SRC:orbmdl/%.cpp=$(OBJDIR)/orbmdl/%.o)
NRLMSISE00_OBJ = $(NRLMSISE00_SRC:orbmdl/nrlmsise-00/%.c=$(OBJDIR)/orbmdl/nrlmsise-00/%.o)
FILTER_OBJ = $(FILTER_SRC:filter/%.cpp=$(OBJDIR)/filter/%.o)

# Libraries
ORBMDL_AR = $(LIBDIR)/liborbmdl.a
NRLMSISE00_AR = $(LIBDIR)/libnrlmsise00.a
FILTER_AR = $(LIBDIR)/libfilter.a

# Compilation flags
CPPFLAGS = --std=c++11 -Wall -pedantic -g -O0 -D_GLIBCXX_DEBUG 
CFLAGS = -Wall -pedantic -g -O3 
LDFLAGS = -L./yaml-cpp/build -L$(LocalDIR)/lib -L$(LIBDIR) -L$(OBJDIR)/orbmdl -L$(OBJDIR)/orbmdl/nrlmsise-00 -L$(OBJDIR)/filter -lorbmdl -lnrlmsise00 -lfilter -lyaml-cpp
INCLUDE = -I$(LocalDIR)/include/eigen3/ -I$(LocalDIR)/boost_1_81_0 -I./orbmdl -I./orbmdl/3rdparty -I./orbmdl/sofa -I./orbmdl/nrlmsise-00 -I./filter
YAML_INCLUDE = -I./yaml-cpp/include

# Default target
all: $(ORBMDL_AR) $(NRLMSISE00_AR) $(FILTER_AR) $(OBJDIR)/scripts/$(FILENAME)

# Build libraries
$(ORBMDL_AR): $(ORBMDL_OBJ) $(OBJDIR)/orbmdl/3rdparty/*.o $(OBJDIR)/orbmdl/sofa/*.o | $(LIBDIR)
	/usr/bin/ar rcs $@ $^

$(NRLMSISE00_AR): $(NRLMSISE00_OBJ) | $(LIBDIR)
	/usr/bin/ar rcs $@ $^

$(FILTER_AR): $(FILTER_OBJ) | $(LIBDIR)
	/usr/bin/ar rcs $@ $^

# Compile object files for libraries
# $(OBJDIR)/filter/%.o: filter/%.cpp | $(OBJDIR)/filter
#   $(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
$(OBJDIR)/filter/%.o: filter/%.cpp $(filter %.cpp,$(wildcard filter/*.cpp)) | $(OBJDIR)/filter
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@
    
$(OBJDIR)/orbmdl/%.o: orbmdl/%.cpp | $(OBJDIR)/orbmdl
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

$(OBJDIR)/orbmdl/3rdparty/%.o: orbmdl/3rdparty/%.cpp | $(OBJDIR)/orbmdl/3rdparty
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

$(OBJDIR)/orbmdl/sofa/%.o: orbmdl/sofa/%.cpp | $(OBJDIR)/orbmdl/sofa
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

$(OBJDIR)/orbmdl/nrlmsise-00/%.o: orbmdl/nrlmsise-00/%.c | $(OBJDIR)/orbmdl/nrlmsise-00
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Build executable
$(OBJDIR)/scripts/$(FILENAME): $(OBJDIR)/scripts/$(FILENAME).o $(ORBMDL_AR) $(NRLMSISE00_AR) $(FILTER_AR)
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(YAML_INCLUDE) $^ $(LDFLAGS) -o $@

$(OBJDIR)/scripts/%.o: scripts/%.cpp | $(OBJDIR)/scripts
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(YAML_INCLUDE) -c $< -o $@

# Clean up build files
clean:
    rm -rf $(OBJDIR)/orbmdl/*.o $(OBJDIR)/orbmdl/3rdparty/*.o $(OBJDIR)/orbmdl/sofa/*.o $(OBJDIR)/orbmdl/nrlmsise-00/*.o $(OBJDIR)/filter/*.o $(OBJDIR)/scripts/*.o $(ORBMDL_AR) $(NRLMSISE00_AR) $(FILTER_AR) $(OBJDIR)/scripts/*
    rm -rf $(OBJDIR)/orbmdl/*.a $(OBJDIR)/orbmdl/nrlmsise-00/*.a $(OBJDIR)/filter/*.a 
    rm -rf orbmdl/*.o orbmdl/3rdparty/*.o orbmdl/sofa/*.o orbmdl/nrlmsise-00/*.o  filter/*.o
