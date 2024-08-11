#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "orbit_propagator_wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(orbit_propagator, m)
{
    py::class_<OrbitPropagatorWapper>(m, "OrbitPropagatorWapper")
        .def(py::init<const std::string &, double, int>())
        .def("propagate", &OrbitPropagatorWapper::propagate)
        .def("readConfigFile", &OrbitPropagatorWapper::readConfigFile)
        .def("initGlobalVariables", &OrbitPropagatorWapper::initGlobalVariables)
        .def("stdVec2EigenVec", &OrbitPropagatorWapper::stdVec2EigenVec)
        .def("initEGMCoef", &OrbitPropagatorWapper::initEGMCoef);
}
