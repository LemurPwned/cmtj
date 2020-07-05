#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "cvector.hpp"

#include "junction.hpp"
#include <stdio.h>
#include <vector>

#define VERSION_INFO 1.0

using namespace pybind11::literals;

namespace py = pybind11;

PYBIND11_MODULE(spinpy, m)
{
    m.doc() = "Python binding of C++ PyMTJ Library";

    py::enum_<Axis>(m, "Axis")
        .value("xaxis", xaxis)
        .value("yaxis", yaxis)
        .value("zaxis", zaxis)
        .export_values();

    py::class_<CVector>(m, "CVector")
        .def(py::init<
             double, double, double>());

    py::implicitly_convertible<std::list<double>, CVector>();
    // py::implicitly_convertible<std::list<std::list<double>>, std::list<CVector>>();

    py::class_<Layer>(m, "Layer")
        .def(py::init<
                 std::string,
                 CVector,
                 CVector,
                 double,
                 double,
                 double,
                 double,
                 double,
                 std::vector<CVector>,
                 std::vector<CVector>>(), );

    py::class_<Junction>(m, "Junction")
        .def(py::init<
             std::vector<Layer>,
             std::string,
             double, double>())

        .def("runSimulation", &Junction::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "persist"_a = false,
             "log"_a = false)

        // set Layer parameters
        .def("setLayerAnisotropy", &Junction::setLayerAnisotropy)
        .def("setLayerCoupling", &Junction::setLayerCoupling)
        .def("setConstantExternalField", &Junction::setConstantExternalField)

        // set updates
        .def("setLayerStepUpdate", &Junction::setLayerStepUpdate)
        .def("setLayerIECUpdate", &Junction::setLayerIECUpdate)
        .def("setLayerStepUpdate", &Junction::setLayerStepUpdate)

        // junction calculations
        .def("calculateMagnetoresistance", &Junction::calculateMagnetoresistance)
        .def("calculateVoltageSpinDiode", &Junction::calculateVoltageSpinDiode,
             "frequency"_a,
             "power"_a = 10e-6,
             "minTime"_a)
        .def("calculateFFT", &Junction::calculateFFT,
             "minTime"_a,
             "timeStep"_a = 1e-11);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
