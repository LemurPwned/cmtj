#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/cvector.hpp"
#include "../core/drivers.hpp"
#include "../core/junction.hpp"
#include <stdio.h>
#include <vector>

using namespace pybind11::literals;

namespace py = pybind11;

#define USING_PY true
PYBIND11_MODULE(cmtj, m)
{
    m.doc() = "Python binding for C++ CMTJ Library";

    // helpers
    m.def("c_dot", &c_dot);
    // Driver Class
    py::class_<ScalarDriver>(m, "ScalarDriver")
        .def_static("getConstantDriver",
                    &ScalarDriver::getConstantDriver,
                    "constantValue"_a)
        .def_static("getPulseDriver",
                    &ScalarDriver::getPulseDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "period"_a,
                    "cycle"_a)
        .def_static("getSineDriver",
                    &ScalarDriver::getSineDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "frequency"_a,
                    "phase"_a)
        .def_static("getStepDriver",
                    &ScalarDriver::getStepDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "timeStart"_a,
                    "timeStop"_a);

    py::class_<NullDriver, ScalarDriver>(m, "NullDriver")
        .def(py::init<>());

    py::class_<AxialDriver>(m, "AxialDriver")
        .def(py::init<ScalarDriver, ScalarDriver, ScalarDriver>())
        .def(py::init<std::vector<ScalarDriver>>())
        .def("getVectorAxialDriver", AxialDriver::getVectorAxialDriver)
        .def("getCurrentAxialDrivers",
             &AxialDriver::getCurrentAxialDrivers)
        .def("applyMask", py::overload_cast<CVector>(&AxialDriver::applyMask))
        .def("applyMask", py::overload_cast<std::vector<unsigned int>>(&AxialDriver::applyMask));
    // CVector
    py::enum_<Axis>(m, "Axis")
        .value("xaxis", xaxis)
        .value("yaxis", yaxis)
        .value("zaxis", zaxis)
        .export_values();

    py::class_<CVector>(m, "CVector")
        .def(py::init<
             double, double, double>())
        .def_readwrite("x", &CVector::x)
        .def_readwrite("y", &CVector::y)
        .def_readwrite("z", &CVector::z)
        .def("length", &CVector::length);

    py::implicitly_convertible<std::list<double>, CVector>();
    py::implicitly_convertible<std::vector<double>, CVector>();

    py::class_<Layer>(m, "Layer")
        .def(py::init<
                 std::string,          // id
                 CVector,              // mag
                 CVector,              // anis
                 double,               // Ms
                 double,               // thickness
                 double,               // cellSurface
                 std::vector<CVector>, // demagTensor
                 std::vector<CVector>, // dipoleTensor
                 double,               // temperature
                 bool,                 // includeSTT
                 double,               // damping
                 double,               // SlonczewskiSpacerLayerParameter
                 double,               // beta
                 double,               // spinPolarisation
                 bool>(),
             "id"_a,
             "mag"_a,
             "anis"_a,
             "Ms"_a,
             "thickness"_a,
             "cellSurface"_a,
             "demagTensor"_a,
             "dipoleTensor"_a,
             "temperature"_a = 0.0,
             "includeSTT"_a = false,
             "damping"_a = 0.011,
             "SlonczewskiSpacerLayerParameter"_a = 1.0,
             "beta"_a = 0.0,
             "spinPolarisation"_a = 0.8,
             "silent"_a = true);

    py::class_<Junction>(m, "Junction")
        .def(py::init<std::vector<Layer>,
                      std::string>(),
             "layers"_a,
             "filename"_a = "")
        .def(py::init<
                 std::vector<Layer>,
                 std::string,
                 double, double>(),
             "layers"_a,
             "filename"_a,
             "Rp"_a = 100,
             "Rap"_a = 200)
        .def(py::init<
                 std::vector<Layer>,
                 std::string, 
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>>(),
             "layers"_a,
             "filename"_a,
             "Rx0"_a,
             "Ry0"_a,
             "AMR_X"_a,
             "AMR_Y"_a,
             "SMR_X"_a,
             "SMR_Y"_a,
             "AHE"_a)
        // log utils
        .def("getLog", &Junction::getLog)
        .def("clearLog", &Junction::clearLog)

        .def("runSimulation", &Junction::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "writeFrequency"_a = 1e-11,
             "persist"_a = false,
             "log"_a = false,
             "calculateEnergies"_a = false)

        // driver setters
        .def("setLayerOerstedFieldDriver", &Junction::setLayerOerstedFieldDriver)
        .def("setLayerExternalFieldDriver", &Junction::setLayerExternalFieldDriver)
        .def("setLayerCurrentDriver", &Junction::setLayerCurrentDriver)
        .def("setLayerAnisotropyDriver", &Junction::setLayerAnisotropyDriver)
        .def("setLayerIECDriver", &Junction::setLayerIECDriver)
        .def("setLayerOerstedFieldDriver", &Junction::setLayerOerstedFieldDriver)
        .def("setLayerMagnetisation", &Junction::setLayerMagnetisation) 
        // junction calculations
        .def("getMagnetoresistance", &Junction::getMagnetoresistance);
}
