// bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "orbit_propagator.h"

namespace py = pybind11;

PYBIND11_MODULE(orbit_propagator, m)
{
    py::class_<OrbitPropagator>(m, "OrbitPropagator")
        .def(py::init<const std::string &>())
        .def("propagate", &OrbitPropagator::propagate);
}