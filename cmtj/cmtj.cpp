#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "cvector.hpp"
#include "parallel.hpp"
#include "junction.hpp"
#include "drivers.hpp"
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

    m.def("RK45", &Junction::RK45);
    m.def("LLG", &Junction::LLG);

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

    py::class_<NullDriver>(m, "NullDriver")
        .def(py::init<>());

    py::class_<AxialDriver>(m, "AxialDriver")
        .def(py::init<std::vector<ScalarDriver>>())
        .def("getCurrentAxialDrivers",
             &AxialDriver::getCurrentAxialDrivers);

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

    py::implicitly_convertible<std::list<float>, CVector>();

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
                 double                // spinPolarisation
                 >(),
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
             "spinPolarisation"_a = 0.8)
        .def("calculateLayerCriticalSwitchingCurrent",
             &Layer::calculateLayerCriticalSwitchingCurrent);

    py::class_<Junction>(m, "Junction")
        .def(py::init<
                 std::vector<Layer>,
                 std::string,
                 double, double>(),
             "layers"_a,
             "filename"_a,
             "Rp"_a,
             "Rap"_a)

        // getters & setters

        // log utils
        .def("getLog", &Junction::getLog)
        .def("clearLog", &Junction::clearLog)

        .def("runSimulation", &Junction::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "persist"_a = false,
             "log"_a = false,
             "calculateEnergies"_a = false)

        // driver setters
        .def("setLayerExternalFieldDriver", &Junction::setLayerExternalFieldDriver)
        .def("setLayerCurrentDriver", &Junction::setLayerCurrentDriver)
        .def("setLayerAnisotropyDriver", &Junction::setLayerAnisotropyDriver)
        .def("setLayerIECDriver", &Junction::setLayerIECDriver)

        // junction calculations
        .def("calculateMagnetoresistance", &Junction::calculateMagnetoresistance)
        .def("calculateVoltageSpinDiode", &Junction::calculateVoltageSpinDiode,
             "frequency"_a,
             "power"_a = 10e-6,
             "minTime"_a)
        .def("calculateFFT", &Junction::calculateFFT,
             "minTime"_a,
             "timeStep"_a = 1e-11);
}
