#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/stack.hpp"
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
using DStack = Stack<double>;
namespace py = pybind11;

#define USING_PY true
PYBIND11_MODULE(cmtj, m)
{
    // helpers
    m.def("c_dot", &c_dot<double>);
    m.doc() = "Python binding for C++ CMTJ Library.";

    // CVector
    py::class_<DVector>(m, "CVector")
        .def(py::init<
             double, double, double>())
        .def_readwrite("x", &DVector::x)
        .def_readwrite("y", &DVector::y)
        .def_readwrite("z", &DVector::z);
    // .def("length", &DVector::length);

    py::implicitly_convertible<std::list<double>, DVector>();
    py::implicitly_convertible<std::vector<double>, DVector>();

    py::enum_<Axis>(m, "Axis")
        .value("xaxis", xaxis)
        .value("yaxis", yaxis)
        .value("zaxis", zaxis)
        .export_values();

    py::enum_<Reference>(m, "Reference")
        .value("none", NONE)
        .value("fixed", FIXED)
        .value("top", TOP)
        .value("bottom", BOTTOM)
        .export_values();

    py::enum_<SolverMode>(m, "SolverMode")
        .value("RK4", RK4)
        .value("EulerHeun", EULER_HEUN)
        .value("DormandPrice", DORMAND_PRICE)
        .export_values();
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
        .def_static("getPosSineDriver",
                    &DScalarDriver::getPosSineDriver,
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
                 double                // damping
                 >(),
             "id"_a,
             "mag"_a,
             "anis"_a,
             "Ms"_a,
             "thickness"_a,
             "cellSurface"_a,
             "demagTensor"_a,
             "damping"_a = 0.011)
        .def_static("createSOTLayer", &DLayer::LayerSOT,
                    "id"_a,
                    "mag"_a,
                    "anis"_a,
                    "Ms"_a,
                    "thickness"_a,
                    "cellSurface"_a,
                    "demagTensor"_a,
                    "damping"_a = 0.011,
                    "fieldLikeTorque"_a = 0.0,
                    "dampingLikeTorque"_a = 0.0)
        .def_static("createSTTLayer", &DLayer::LayerSTT,
                    "id"_a,
                    "mag"_a,
                    "anis"_a,
                    "Ms"_a,
                    "thickness"_a,
                    "cellSurface"_a,
                    "demagTensor"_a,
                    "damping"_a = 0.011,
                    "SlonczewskiSpacerLayerParameter"_a = 1.0,
                    "beta"_a = 0.0,
                    "spinPolarisation"_a = 0.0)
        .def("setMagnetisation", &DLayer::setMagnetisation)
        .def("setAnisotropyDriver", &DLayer::setAnisotropyDriver)
        .def("setExternalFieldDriver", &DLayer::setExternalFieldDriver)
        .def("setOerstedFieldDriver", &DLayer::setOerstedFieldDriver)
        // reference layers
        .def("setReferenceLayer", py::overload_cast<DVector>(&DLayer::setReferenceLayer))
        .def("setReferenceLayer", py::overload_cast<Reference>(&DLayer::setReferenceLayer))

        .def("setFieldLikeTorqueDriver", &DLayer::setFieldLikeTorqueDriver)
        .def("setDampingLikeTorqueDriver", &DLayer::setDampingLikeTorqueDriver)
        .def("setTemperatureDriver", &DLayer::setTemperatureDriver)
        .def("setTopDipoleTensor", &DLayer::setTopDipoleTensor)
        .def("setBottomDipoleTensor", &DLayer::setBottomDipoleTensor);

    py::class_<DJunction>(m, "Junction")
        .def(py::init<std::vector<DLayer>>(),
             "layers"_a)
        .def(py::init<std::vector<DLayer>,
                      double, double>(),
             "layers"_a,
             "Rp"_a = 100,
             "Rap"_a = 200)
        .def(py::init<
                 std::vector<DLayer>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>,
                 std::vector<double>>(),
             "layers"_a,
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
        .def("saveLog", &DJunction::saveLogs, "filename"_a)

        .def("runSimulation", &DJunction::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "writeFrequency"_a = 1e-11,
             "log"_a = false,
             "calculateEnergies"_a = false,
             "solverMode"_a = RK4)

        // driver setters
        .def("setLayerOerstedFieldDriver", &DJunction::setLayerOerstedFieldDriver)
        .def("setLayerExternalFieldDriver", &DJunction::setLayerExternalFieldDriver)
        .def("setLayerCurrentDriver", &DJunction::setLayerCurrentDriver)
        .def("setLayerAnisotropyDriver", &DJunction::setLayerAnisotropyDriver)
        .def("setIECDriver", &DJunction::setIECDriver)
        .def("setQuadIECDriver", &DJunction::setQuadIECDriver)
        .def("setLayerOerstedFieldDriver", &DJunction::setLayerOerstedFieldDriver)
        .def("setLayerMagnetisation", &DJunction::setLayerMagnetisation)
        // temp
        .def("setLayerTemperatureDriver", &DJunction::setLayerTemperatureDriver)
        .def("setLayerNonStochasticLangevinDriver", &DJunction::setLayerNonStochasticLangevinDriver)
        // SOT setters
        .def("setLayerFieldLikeTorqueDriver", &DJunction::setLayerFieldLikeTorqueDriver)
        .def("setLayerDampingLikeTorqueDriver", &DJunction::setLayerDampingLikeTorqueDriver)
        // Reference setters
        .def("setLayerReferenceType", &DJunction::setLayerReferenceType)
        .def("setLayerReferenceLayer", &DJunction::setLayerReferenceLayer)
        // junction calculations
        .def("getLayerMagnetisation", &DJunction::getLayerMagnetisation)
        .def("getMagnetoresistance", &DJunction::getMagnetoresistance);

    // stack module
    py::module stack_module = m.def_submodule("stack", "A stack submodule for joining MTJ junctions");
    py::class_<SeriesStack<double>>(stack_module, "SeriesStack")
        .def(py::init<std::vector<DJunction>>(), "junctionList"_a)
        .def("runSimulation", &SeriesStack<double>::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "writeFrequency"_a = 1e-11)
        .def("setMagnetistation", &SeriesStack<double>::setMagnetisation, "juncionId"_a, "layerId"_a, "mag"_a)
        .def("setCoupledCurrentDriver", &SeriesStack<double>::setCoupledCurrentDriver, "driver"_a)
        .def("setExternalFieldDriver", &SeriesStack<double>::setExternalFieldDriver, "driver"_a)
        .def("setCouplingStrength", &SeriesStack<double>::setCouplingStrength, "coupling"_a)
        // logging
        .def("clearLogs", &ParallelStack<double>::clearLogs)
        .def("getLog", py::overload_cast<unsigned int>(&SeriesStack<double>::getLog))
        .def("getLog", py::overload_cast<>(&SeriesStack<double>::getLog));
    py::class_<ParallelStack<double>>(stack_module, "ParallelStack")
        .def(py::init<std::vector<DJunction>>(), "junctionList"_a)
        .def("runSimulation", &ParallelStack<double>::runSimulation,
             "totalTime"_a,
             "timeStep"_a = 1e-13,
             "writeFrequency"_a = 1e-11)
        .def("setMagnetistation", &ParallelStack<double>::setMagnetisation, "juncionId"_a, "layerId"_a, "mag"_a)
        .def("setCoupledCurrentDriver", &ParallelStack<double>::setCoupledCurrentDriver, "driver"_a)
        .def("setExternalFieldDriver", &ParallelStack<double>::setExternalFieldDriver, "driver"_a)
        .def("setCouplingStrength", &ParallelStack<double>::setCouplingStrength, "coupling"_a)
        // logging
        .def("clearLogs", &ParallelStack<double>::clearLogs)
        .def("getLog", py::overload_cast<unsigned int>(&ParallelStack<double>::getLog))
        .def("getLog", py::overload_cast<>(&ParallelStack<double>::getLog));
}
