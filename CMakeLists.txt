cmake_minimum_required(VERSION 3.12)
project(orbit_propagator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set SDKROOT
execute_process(COMMAND xcrun --sdk macosx --show-sdk-path OUTPUT_VARIABLE SDKROOT OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get pybind11 include directories
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE PYBIND11_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)

# Split PYBIND11_INCLUDE_DIRS into a list
string(REPLACE " " ";" PYBIND11_INCLUDE_DIRS_LIST ${PYBIND11_INCLUDE_DIRS})

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -isystem ${SDKROOT}/usr/include/c++/v1 -isysroot ${SDKROOT} -nostdinc++ -Wall -pedantic -g -Wno-error -v")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -isysroot ${SDKROOT} -Wall -pedantic -g")

# Include directories
include_directories(
    orbmdl
    orbmdl/3rdparty
    orbmdl/sofa
    orbmdl/nrlmsise-00
    filter
    /usr/local/include
    /usr/local/include/eigen3
    /usr/local/boost_1_86_0
    ${PYBIND11_INCLUDE_DIRS_LIST}
    /Users/yangyang/Documents/GitHub/ordDetHOUSEPublished/orbDetHOUSE/odvenv/lib/python3.10/site-packages/pybind11/include
    /Library/Frameworks/Python.framework/Versions/3.10/include/python3.10
    ${SDKROOT}/usr/include
    ${SDKROOT}/usr/include/c++/v1
)

# Source files
file(GLOB ORBMDL_SRC "orbmdl/*.cpp" "orbmdl/3rdparty/*.cpp" "orbmdl/sofa/*.cpp")
file(GLOB NRLMSISE00_SRC "orbmdl/nrlmsise-00/*.c")
file(GLOB FILTER_SRC "filter/*.cpp")

# Create static libraries
add_library(orbmdl STATIC ${ORBMDL_SRC})
add_library(nrlmsise00 STATIC ${NRLMSISE00_SRC})
add_library(filter STATIC ${FILTER_SRC})

# Find yaml-cpp
find_package(yaml-cpp REQUIRED)

# Create shared library
add_library(orbit_propagator_wrapper SHARED orbmdl/orbit_propagator_wrapper_py.cpp orbmdl/orbit_propagator_wrapper.cpp)
target_link_libraries(orbit_propagator_wrapper orbmdl nrlmsise00 filter yaml-cpp::yaml-cpp python3.10)

# Set output directories
set_target_properties(orbmdl PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set_target_properties(nrlmsise00 PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set_target_properties(filter PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set_target_properties(orbit_propagator_wrapper PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)