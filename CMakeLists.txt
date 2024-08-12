# didn't work

cmake_minimum_required(VERSION 3.12)
project(orbit_propagator)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add pybind11 installation path to CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH "/Users/yangyang/Documents/GitHub/ordDetHOUSEPublished/orbDetHOUSE/odvenv/lib/python3.8/site-packages/pybind11/include")

# Set pybind11_DIR to the directory containing pybind11Config.cmake
set(pybind11_DIR "/Users/yangyang/Documents/GitHub/ordDetHOUSEPublished/orbDetHOUSE/odvenv/lib/python3.8/site-packages/pybind11/share/cmake/pybind11")

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

# Get SDK path using CMake's built-in command
execute_process(COMMAND xcrun --sdk macosx --show-sdk-path OUTPUT_VARIABLE SDK_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)

# Ensure correct paths for C++ standard library
include_directories(SYSTEM ${SDK_PATH}/usr/include/c++/v1)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/orbmdl
    ${CMAKE_SOURCE_DIR}/orbmdl/3rdparty
    ${CMAKE_SOURCE_DIR}/orbmdl/sofa
    ${CMAKE_SOURCE_DIR}/orbmdl/nrlmsise-00
    ${CMAKE_SOURCE_DIR}/filter
    /usr/local/include
    /usr/local/include/eigen3
    /usr/local/boost_1_81_0
    /Users/yangyang/Documents/GitHub/ordDetHOUSEPublished/orbDetHOUSE/odvenv/lib/python3.8/site-packages/pybind11/include
)

# Source files
file(GLOB ORBMDL_SRC "orbmdl/*.cpp" "orbmdl/3rdparty/*.cpp" "orbmdl/sofa/*.cpp")
file(GLOB NRLMSISE00_SRC "orbmdl/nrlmsise-00/*.c")
file(GLOB FILTER_SRC "filter/*.cpp")

# Add libraries
add_library(orbmdl STATIC ${ORBMDL_SRC})
add_library(nrlmsise00 STATIC ${NRLMSISE00_SRC})
add_library(filter STATIC ${FILTER_SRC})

# Add the main library
add_library(orbit_propagator_wrapper MODULE orbmdl/orbit_propagator_wrapper_py.cpp orbmdl/orbit_propagator_wrapper.cpp)
target_include_directories(orbit_propagator_wrapper PRIVATE ${PYBIND11_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(orbit_propagator_wrapper PRIVATE pybind11::module orbmdl nrlmsise00 filter yaml-cpp)

# Set compiler and linker flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${SDK_PATH}/usr/include/c++/v1 -isysroot ${SDK_PATH} -stdlib=libc++ -Wno-error")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -g -isysroot ${SDK_PATH}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${SDK_PATH}/usr/lib -lSystem")

# Add executable if needed
# add_executable(main_executable main.cpp)
# target_link_libraries(main_executable PRIVATE orbit_propagator_wrapper)