CXX = g++-11
CC = gcc-11

# create directories
$(shell mkdir -p bin/orbmdl)
$(shell mkdir -p bin/orbmdl/3rdparty)
$(shell mkdir -p bin/orbmdl/sofa)
$(shell mkdir -p bin/orbmdl/nrlmsise-00)
$(shell mkdir -p bin/filter)
$(shell mkdir -p bin/scripts)

OBJDIR = bin
LocalDIR = /usr/local
FUNDIR = scripts

# for compiling orbmdl and filter
CPPFLAGS = --std=c++11 -Wall -pedantic -g -O3 
CFLAGS = -Wall -pedantic -g -O3 
LDFLAGS = -L./yaml-cpp/build -L$(LocalDIR)/lib -L$(OBJDIR)/filter -L$(OBJDIR)/orbmdl -lfilter -lorbmdl -lnrlmsise00 -lyaml-cpp
INCLUDE = -I$(LocalDIR)/include/eigen3/ -I$(LocalDIR)/boost_1_81_0 -I./orbmdl -I./orbmdl/3rdparty -I./orbmdl/sofa -I./orbmdl/nrlmsise-00 -I./filter
YAML_INCLUDE = -I./yaml-cpp/include

# orbmdl files
ORBMDL_SRC = $(wildcard orbmdl/*.cpp) $(wildcard orbmdl/3rdparty/*.cpp) $(wildcard orbmdl/sofa/*.cpp)
ORBMDL_OBJ = $(ORBMDL_SRC:orbmdl/%.cpp=$(OBJDIR)/orbmdl/%.o)
ORBMDL_AR = $(OBJDIR)/orbmdl/liborbmdl.a

# nrlmsise-00 (drag density model) files
NRLMSISE00_SRC = orbmdl/nrlmsise-00/nrlmsise-00.c orbmdl/nrlmsise-00/nrlmsise-00_data.c
NRLMSISE00_OBJ = $(NRLMSISE00_SRC:orbmdl/nrlmsise-00/%.c=$(OBJDIR)/orbmdl/nrlmsise-00/%.o)
NRLMSISE00_AR = $(OBJDIR)/orbmdl/libnrlmsise00.a

# targets
all: $(ORBMDL_AR) $(NRLMSISE00_AR) $(OBJDIR)/scripts/$(FILENAME)

$(ORBMDL_AR): $(ORBMDL_OBJ)
	ar rcs $@ $^

$(NRLMSISE00_AR): $(NRLMSISE00_OBJ)
	ar rcs $@ $^

$(OBJDIR)/orbmdl/%.o: orbmdl/%.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

$(OBJDIR)/orbmdl/nrlmsise-00/%.o: orbmdl/nrlmsise-00/%.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# target for generating an executable
$(OBJDIR)/scripts/%: $(OBJDIR)/scripts/%.o $(ORBMDL_AR) $(NRLMSISE00_AR)
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(YAML_INCLUDE) $^ $(LDFLAGS) -o $@

$(OBJDIR)/scripts/%.o: scripts/%.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(YAML_INCLUDE) -c $< -o $@

# target for generating an executable based on FILENAME
$(OBJDIR)/scripts/$(FILENAME): $(OBJDIR)/scripts/$(FILENAME).o $(ORBMDL_AR) $(NRLMSISE00_AR)
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(YAML_INCLUDE) $^ $(LDFLAGS) -o $@

clean:
	rm -rf $(OBJDIR)/orbmdl/*.o $(OBJDIR)/orbmdl/nrlmsise-00/*.o $(OBJDIR)/scripts/*.o $(ORBMDL_AR) $(NRLMSISE00_AR) $(OBJDIR)/scripts/*