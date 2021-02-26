#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/cvector.hpp"
#include "../core/drivers.hpp"
#include "../core/junction.hpp"
#include "../core/parallel.hpp"
#include <stdio.h>
#include <vector>

using namespace pybind11::literals;

namespace py = pybind11;

#define USING_PY true

static std::map<std::string, std::vector<double>> parallelGILWrapper(Junction &mtj,
                                                                     double minField,
                                                                     double maxField,
                                                                     int numberOfPoints,
                                                                     int threadNumber,
                                                                     const std::function<fnRes(Junction &mtj,
                                                                                               const double scanningParam)>
                                                                         runnableFunction)
{
    std::map<std::string, std::vector<double>> res;
    pybind11::gil_scoped_acquire acquire;
    {
        pybind11::gil_scoped_release release;
        res = ComputeUtil::parallelFieldScan(mtj, minField, maxField, numberOfPoints, threadNumber, runnableFunction);
    }

    return res;
}

PYBIND11_MODULE(cmtj, m)
{
    m.doc() = "Python binding for C++ CMTJ Library";

    // helpers
    m.def("c_dot", &c_dot);
    m.def("customResultMap", &ComputeUtil::customResultMap,
          "resultMap"_a,
          "filename"_a);

    m.def("parallelFieldScan", &parallelGILWrapper,
          pybind11::call_guard<pybind11::gil_scoped_release>(),
          "mtj"_a,
          "minField"_a,
          "maxField"_a,
          "numberOfPoints"_a,
          "numberOfThreads"_a,
          "runnableFunction"_a);

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
                 double,               // cekkSurface
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
        .def(py::init<
                 std::vector<Layer>,
                 std::string,
                 double, double>(),
             "layers"_a,
             "filename"_a,
             "Rp"_a,
             "Rap"_a)
        .def(py::init<
                 std::vector<Layer>,
                 std::string, std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>>(),
             "layers"_a,
             "filename"_a,
             "Rx0"_a,
             "Ry0"_a,
             "AMR"_a,
             "AHE"_a,
             "SMR"_a)
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

        // junction calculations
        .def("advancedMagnetoResistance", &Junction::advancedMagnetoResistance)
        .def("calculateMagnetoresistance", &Junction::calculateMagnetoresistance)
        .def("getMagnetoresistance", &Junction::getMagnetoresistance)
        .def("calculateVoltageSpinDiode", &Junction::calculateVoltageSpinDiode,
             "frequency"_a,
             "power"_a = 10e-6,
             "minTime"_a)
        .def("calculateFFT", &Junction::calculateFFT,
             "minTime"_a,
             "timeStep"_a = 1e-11);
}
