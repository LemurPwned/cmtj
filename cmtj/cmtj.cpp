#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "cvector.hpp"
#include "parallel.hpp"
#include "junction.hpp"
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

    py::enum_<Axis>(m, "Axis")
        .value("xaxis", xaxis)
        .value("yaxis", yaxis)
        .value("zaxis", zaxis)
        .export_values();

    py::class_<CVector>(m, "CVector")
        .def(py::init<
             double, double, double>());

    py::implicitly_convertible<std::list<double>, CVector>();

    py::class_<Layer>(m, "Layer")
        .def(py::init<
                 std::string,          // id
                 CVector,              // mag
                 CVector,              // anis
                 double,               // K
                 double,               // Ms
                 double,               // J
                 double,               // thickness
                 double,               // cekkSurface
                 std::vector<CVector>, // demagTensor
                 std::vector<CVector>, // dipoleTensor
                 double,               // temperature
                 bool,                 // includeSTT
                 double,               // damping
                 double,               // currentDensity
                 double,               // SlonczewskiSpacerLayerParameter
                 double,                // beta
                 double                // spinPolarisation
                 >(),
             "id"_a,
             "mag"_a,
             "anis"_a,
             "K"_a,
             "Ms"_a,
             "J"_a,
             "thickness"_a,
             "cellSurface"_a,
             "demagTensor"_a,
             "dipoleTensor"_a,
             "temperature"_a = 0.0,
             "includeSTT"_a = false,
             "damping"_a = 0.011,
             "currentDensity"_a = 1,
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
             "log"_a = false)

        // set Layer parameters
        .def("setLayerAnisotropy", &Junction::setLayerAnisotropy)
        .def("setLayerCoupling", &Junction::setLayerCoupling)
        .def("setConstantExternalField", &Junction::setConstantExternalField)

        // set updates
        .def("setLayerStepUpdate", &Junction::setLayerStepUpdate)
        .def("setLayerIECUpdate", &Junction::setLayerIECUpdate)
        .def("setLayerStepUpdate", &Junction::setLayerStepUpdate)
        .def("setLayerCurrentDensity", &Junction::setLayerCurrentDensity)

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
