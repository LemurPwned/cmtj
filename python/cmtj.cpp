#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/cvector.hpp"
#include "../core/drivers.hpp"
#include "../core/junction.hpp"
#include <stdio.h>
#include <vector>

using namespace pybind11::literals;
using DJunction = Junction<double>;
using DLayer = Layer<double>;
using DVector = CVector<double>;
using DScalarDriver = ScalarDriver<double>;
using DAxialDriver = AxialDriver<double>;
using DNullDriver = NullDriver<double>;

namespace py = pybind11;

#define USING_PY true
PYBIND11_MODULE(cmtj, m)
{
    // helpers
    m.def("c_dot", &c_dot<double>);
    m.doc() = "Python binding for C++ CMTJ Library";
    py::enum_<Axis>(m, "Axis")
        .value("xaxis", xaxis)
        .value("yaxis", yaxis)
        .value("zaxis", zaxis)
        .export_values();

    // CVector
    py::class_<DVector>(m, "CVector")
        .def(py::init<
             double, double, double>())
        .def_readwrite("x", &DVector::x)
        .def_readwrite("y", &DVector::y)
        .def_readwrite("z", &DVector::z)
        .def("length", &DVector::length);

    py::implicitly_convertible<std::list<double>, DVector>();
    py::implicitly_convertible<std::vector<double>, DVector>();

    // Driver Class
    py::class_<DScalarDriver>(m, "ScalarDriver")
        .def_static("getConstantDriver",
                    &DScalarDriver::getConstantDriver,
                    "constantValue"_a)
        .def_static("getPulseDriver",
                    &DScalarDriver::getPulseDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "period"_a,
                    "cycle"_a)
        .def_static("getSineDriver",
                    &DScalarDriver::getSineDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "frequency"_a,
                    "phase"_a)
        .def_static("getStepDriver",
                    &DScalarDriver::getStepDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "timeStart"_a,
                    "timeStop"_a);

    py::class_<DNullDriver, DScalarDriver>(m, "NullDriver")
        .def(py::init<>());

    py::class_<DAxialDriver>(m, "AxialDriver")
        .def(py::init<DScalarDriver, DScalarDriver, DScalarDriver>())
        .def(py::init<std::vector<ScalarDriver<double>>>())
        .def("getVectorAxialDriver", DAxialDriver::getVectorAxialDriver)
        .def("getCurrentAxialDrivers",
             &DAxialDriver::getCurrentAxialDrivers)
        .def("applyMask", py::overload_cast<DVector>(&DAxialDriver::applyMask))
        .def("applyMask", py::overload_cast<std::vector<unsigned int>>(&DAxialDriver::applyMask));

    py::class_<DLayer>(m, "Layer")
        .def(py::init<
                 std::string,          // id
                 DVector,              // mag
                 DVector,              // anis
                 double,               // Ms
                 double,               // thickness
                 double,               // cellSurface
                 std::vector<DVector>, // demagTensor
                 std::vector<DVector>, // dipoleTensor
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
             "silent"_a = true)
        .def("setMagnetisation", &DLayer::setMagnetisation)
        .def("setAnisotropyDriver", &DLayer::setAnisotropyDriver)
        .def("setExternalFieldDriver", &DLayer::setExternalFieldDriver)
        .def("setOerstedFieldDriver", &DLayer::setOerstedFieldDriver)
        .def("setReferenceLayer", &DLayer::setReferenceLayer);

    py::class_<DJunction>(m, "Junction")
        .def(py::init<std::vector<DLayer>,
                      std::string>(),
             "layers"_a,
             "filename"_a = "")
        .def(py::init<
                 std::vector<DLayer>,
                 std::string,
                 double, double>(),
             "layers"_a,
             "filename"_a,
             "Rp"_a = 100,
             "Rap"_a = 200)
        .def(py::init<
                 std::vector<DLayer>,
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
        .def("getLog", &DJunction::getLog)
        .def("clearLog", &DJunction::clearLog)

        .def("runSimulation", &DJunction::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "writeFrequency"_a = 1e-11,
             "persist"_a = false,
             "log"_a = false,
             "calculateEnergies"_a = false)

        // driver setters
        .def("setLayerOerstedFieldDriver", &DJunction::setLayerOerstedFieldDriver)
        .def("setLayerExternalFieldDriver", &DJunction::setLayerExternalFieldDriver)
        .def("setLayerCurrentDriver", &DJunction::setLayerCurrentDriver)
        .def("setLayerAnisotropyDriver", &DJunction::setLayerAnisotropyDriver)
        .def("setIECDriver", &DJunction::setIECDriver)
        .def("setLayerOerstedFieldDriver", &DJunction::setLayerOerstedFieldDriver)
        .def("setLayerMagnetisation", &DJunction::setLayerMagnetisation)
        // junction calculations
        .def("getLayerMagnetisation", &DJunction::getLayerMagnetisation)
        .def("getMagnetoresistance", &DJunction::getMagnetoresistance);
}
