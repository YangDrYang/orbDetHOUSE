# Makefile for orbit_propagator project

# Compiler and flags
CXX = clang++

# Set SDKROOT
SDKROOT = $(shell xcrun --sdk macosx --show-sdk-path)

# Compiler and flags
CXXFLAGS = -std=c++11 -isystem $(SDKROOT)/usr/include/c++/v1 -isysroot $(SDKROOT) -Wall -pedantic -g -Wno-error -v
CFLAGS = -isysroot $(SDKROOT) -Wall -pedantic -g
LDFLAGS = -L$(SDKROOT)/usr/lib -lSystem

# Include directories
INCLUDES = -Iorbmdl -Iorbmdl/3rdparty -Iorbmdl/sofa -Iorbmdl/nrlmsise-00 -Ifilter \
           -I/usr/local/include -I/usr/local/include/eigen3 -I/usr/local/boost_1_81_0 \
           -I/Users/yangyang/Documents/GitHub/ordDetHOUSEPublished/orbDetHOUSE/odvenv/lib/python3.8/site-packages/pybind11/include \
           -I/Library/Frameworks/Python.framework/Versions/3.8/include/python3.8 \
           -I$(SDKROOT)/usr/include \
           -I$(SDKROOT)/usr/include/c++/v1

# Source files
ORBMDL_SRC = $(wildcard orbmdl/*.cpp orbmdl/3rdparty/*.cpp orbmdl/sofa/*.cpp)
NRLMSISE00_SRC = $(wildcard orbmdl/nrlmsise-00/*.c)
FILTER_SRC = $(wildcard filter/*.cpp)

# Object files
ORBMDL_OBJ = $(ORBMDL_SRC:.cpp=.o)
NRLMSISE00_OBJ = $(NRLMSISE00_SRC:.c=.o)
FILTER_OBJ = $(FILTER_SRC:.cpp=.o)

# Targets
all: orbit_propagator_wrapper.so

# Build libraries
liborbmdl.a: $(ORBMDL_OBJ)
	ar rcs $@ $^

libnrlmsise00.a: $(NRLMSISE00_OBJ)
	ar rcs $@ $^

libfilter.a: $(FILTER_OBJ)
	ar rcs $@ $^

# Build the main library
orbit_propagator_wrapper.so: orbmdl/orbit_propagator_wrapper_py.o orbmdl/orbit_propagator_wrapper.o liborbmdl.a libnrlmsise00.a libfilter.a
	$(CXX) -shared -o $@ $^ $(LDFLAGS) -lyaml-cpp

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean up
clean:
	rm -f $(ORBMDL_OBJ) $(NRLMSISE00_OBJ) $(FILTER_OBJ) orbmdl/orbit_propagator_wrapper_py.o orbmdl/orbit_propagator_wrapper.o liborbmdl.a libnrlmsise00.a libfilter.a orbit_propagator_wrapper.so