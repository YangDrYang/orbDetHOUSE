#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "orbit_propagator_wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(orbit_propagator_wrapper, m)
{
    py::class_<OrbitPropagatorWrapper>(m, "OrbitPropagatorWrapper")
        .def(py::init<const std::string &>())
        .def("propagateOrbit", &OrbitPropagatorWrapper::propagateOrbit)
        .def("saveResults", &OrbitPropagatorWrapper::saveResults)
        .def("readConfigFile", &OrbitPropagatorWrapper::readConfigFile)
        .def("initGlobalVariables", &OrbitPropagatorWrapper::initGlobalVariables)
        .def("stdVec2EigenVec", &OrbitPropagatorWrapper::stdVec2EigenVec)
        .def("initEGMCoef", &OrbitPropagatorWrapper::initEGMCoef);
}
