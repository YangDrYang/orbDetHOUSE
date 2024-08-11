cmake_minimum_required(VERSION 3.12)
project(orbit_propagator)

set(CMAKE_CXX_STANDARD 11)

# Add pybind11 installation path to CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH "/Users/yangyang/Documents/GitHub/ordDetHOUSEPublished/orbDetHOUSE/odvenv/lib/python3.8/site-packages/pybind11")

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(orbit_propagator MODULE orbit_propagator_wrapper_py.cpp orbit_propagator_wrapper.cpp)
target_include_directories(orbit_propagator PRIVATE ${PYBIND11_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(orbit_propagator PRIVATE pybind11::module)