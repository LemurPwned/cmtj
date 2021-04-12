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
    m.def("c_dot", &c_dot<double>);
    // Driver Class
    py::class_<ScalarDriver<double>>(m, "ScalarDriver")
        .def_static("getConstantDriver",
                    &ScalarDriver<double>::getConstantDriver,
                    "constantValue"_a)
        .def_static("getPulseDriver",
                    &ScalarDriver<double>::getPulseDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "period"_a,
                    "cycle"_a)
        .def_static("getSineDriver",
                    &ScalarDriver<double>::getSineDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "frequency"_a,
                    "phase"_a)
        .def_static("getStepDriver",
                    &ScalarDriver<double>::getStepDriver,
                    "constantValue"_a,
                    "amplitude"_a,
                    "timeStart"_a,
                    "timeStop"_a);

    py::class_<NullDriver<double>, ScalarDriver<double>>(m, "NullDriver")
        .def(py::init<>());

    py::class_<AxialDriver<double>>(m, "AxialDriver")
        .def(py::init<ScalarDriver<double>, ScalarDriver<double>, ScalarDriver<double>>())
        .def(py::init<std::vector<ScalarDriver<double>>>())
        .def("getVectorAxialDriver", AxialDriver<double>::getVectorAxialDriver)
        .def("getCurrentAxialDrivers",
             &AxialDriver<double>::getCurrentAxialDrivers)
        .def("applyMask", py::overload_cast<CVector<double>>(&AxialDriver<double>::applyMask))
        .def("applyMask", py::overload_cast<std::vector<unsigned int>>(&AxialDriver<double>::applyMask));
    // CVector
    py::enum_<Axis>(m, "Axis")
        .value("xaxis", xaxis)
        .value("yaxis", yaxis)
        .value("zaxis", zaxis)
        .export_values();

    py::class_<CVector<double>>(m, "CVector")
        .def(py::init<
             double, double, double>())
        .def_readwrite("x", &CVector<double>::x)
        .def_readwrite("y", &CVector<double>::y)
        .def_readwrite("z", &CVector<double>::z)
        .def("length", &CVector<double>::length);

    py::implicitly_convertible<std::list<double>, CVector<double>>();
    py::implicitly_convertible<std::vector<double>, CVector<double>>();

    py::class_<Layer<double>>(m, "Layer")
        .def(py::init<
                 std::string,                  // id
                 CVector<double>,              // mag
                 CVector<double>,              // anis
                 double,                       // Ms
                 double,                       // thickness
                 double,                       // cellSurface
                 std::vector<CVector<double>>, // demagTensor
                 std::vector<CVector<double>>, // dipoleTensor
                 double,                       // temperature
                 bool,                         // includeSTT
                 double,                       // damping
                 double,                       // SlonczewskiSpacerLayerParameter
                 double,                       // beta
                 double,                       // spinPolarisation
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
        .def("setMagnetisation", &Layer<double>::setMagnetisation)
        .def("setAnisotropyDriver", &Layer<double>::setAnisotropyDriver)
        .def("setExternalFieldDriver", &Layer<double>::setExternalFieldDriver)
        .def("setOerstedFieldDriver", &Layer<double>::setOerstedFieldDriver)
        .def("setReferenceLayer", &Layer<double>::setReferenceLayer);

    py::class_<Junction<double>>(m, "Junction")
        .def(py::init<std::vector<Layer<double>>,
                      std::string>(),
             "layers"_a,
             "filename"_a = "")
        .def(py::init<
                 std::vector<Layer<double>>,
                 std::string,
                 double, double>(),
             "layers"_a,
             "filename"_a,
             "Rp"_a = 100,
             "Rap"_a = 200)
        .def(py::init<
                 std::vector<Layer<double>>,
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
        .def("getLog", &Junction<double>::getLog)
        .def("clearLog", &Junction<double>::clearLog)

        .def("runSimulation", &Junction<double>::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "writeFrequency"_a = 1e-11,
             "persist"_a = false,
             "log"_a = false,
             "calculateEnergies"_a = false)

        // driver setters
        .def("setLayerOerstedFieldDriver", &Junction<double>::setLayerOerstedFieldDriver)
        .def("setLayerExternalFieldDriver", &Junction<double>::setLayerExternalFieldDriver)
        .def("setLayerCurrentDriver", &Junction<double>::setLayerCurrentDriver)
        .def("setLayerAnisotropyDriver", &Junction<double>::setLayerAnisotropyDriver)
        .def("setIECDriver", &Junction<double>::setIECDriver)
        .def("setLayerOerstedFieldDriver", &Junction<double>::setLayerOerstedFieldDriver)
        .def("setLayerMagnetisation", &Junction<double>::setLayerMagnetisation)
        // junction calculations
        .def("getLayerMagnetisation", &Junction<double>::getLayerMagnetisation)
        .def("getMagnetoresistance", &Junction<double>::getMagnetoresistance);
}
