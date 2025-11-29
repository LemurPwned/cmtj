#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/cvector.hpp"
#include "../core/drivers.hpp"
#include "../core/junction.hpp"
#include "../core/llgb.hpp"
#include "../core/noise.hpp"
#include "../core/reservoir.hpp"
#include "../core/stack.hpp"
#include "../core/constants.hpp"
#include <stdio.h>
#include <vector>

using namespace pybind11::literals;
namespace py = pybind11;

using DJunction = Junction<double>;
using DParallelStack = ParallelStack<double>;
using DSeriesStack = SeriesStack<double>;
using DLayer = Layer<double>;
using DVector = CVector<double>;
using DScalarDriver = ScalarDriver<double>;
using DAxialDriver = AxialDriver<double>;
using DNullDriver = NullDriver<double>;
using DStack = Stack<double>;
using DLayerMatrix = std::vector<std::vector<DLayer>>;
using DVectorMatrix = std::vector<std::vector<DVector>>;
using DLLGBLayer = LLGBLayer<double>;
using DLLGBJunction = LLGBJunction<double>;

#define USING_PY true
PYBIND11_MODULE(_cmtj, m) {
     // helpers
     m.def("c_dot", &c_dot<double>);
     m.doc() = "Python binding for C++ CMTJ Library.";

     // driver aliases
     m.def(
          "constantDriver",
          [](double value) { return DScalarDriver::getConstantDriver(value); },
          "value"_a);
     m.def(
          "pulseDriver",
          [](double constantValue, double amplitude, double period, double cycle) {
               return DScalarDriver::getPulseDriver(constantValue, amplitude, period,
                    cycle);
          },
          "constantValue"_a, "amplitude"_a, "period"_a, "cycle"_a);
     m.def(
          "sineDriver",
          [](double constantValue, double amplitude, double frequency,
               double phase) {
                    return DScalarDriver::getSineDriver(constantValue, amplitude, frequency,
                         phase);
          },
          "constantValue"_a, "amplitude"_a, "frequency"_a, "phase"_a);
     m.def(
          "posSineDriver",
          [](double constantValue, double amplitude, double frequency,
               double phase) {
                    return DScalarDriver::getPosSineDriver(constantValue, amplitude,
                         frequency, phase);
          },
          "constantValue"_a, "amplitude"_a, "frequency"_a, "phase"_a);
     m.def(
          "stepDriver",
          [](double constantValue, double amplitude, double timeStart,
               double timeStop) {
                    return DScalarDriver::getStepDriver(constantValue, amplitude, timeStart,
                         timeStop);
          },
          "constantValue"_a, "amplitude"_a, "timeStart"_a, "timeStop"_a);
     m.def(
          "trapezoidDriver",
          [](double constantValue, double amplitude, double timeStart,
               double edgeTime, double steadyTime) {
                    return DScalarDriver::getTrapezoidDriver(
                         constantValue, amplitude, timeStart, edgeTime, steadyTime);
          },
          "constantValue"_a, "amplitude"_a, "timeStart"_a, "edgeTime"_a,
          "steadyTime"_a);
     m.def(
          "gaussianImpulseDriver",
          [](double constantValue, double amplitude, double t0, double sigma) {
               return DScalarDriver::getGaussianImpulseDriver(constantValue, amplitude,
                    t0, sigma);
          },
          "constantValue"_a, "amplitude"_a, "t0"_a, "sigma"_a);
     m.def(
          "gaussianStepDriver",
          [](double constantValue, double amplitude, double t0, double sigma) {
               return DScalarDriver::getGaussianStepDriver(constantValue, amplitude,
                    t0, sigma);
          },
          "constantValue"_a, "amplitude"_a, "t0"_a, "sigma"_a);

     // CVector
     py::class_<DVector>(m, "CVector")
          .def(py::init<double, double, double>())
          .def_readwrite("x", &DVector::x)
          .def_readwrite("y", &DVector::y)
          .def_readwrite("z", &DVector::z)
          .def("length", [](const DVector& vec) { return vec.length(); })
          .def("normalize", &DVector::normalize)
          .def("tolist", &DVector::tolist)
          // operators
          .def(py::self + py::self)
          .def(py::self += py::self)
          .def(py::self - py::self)
          .def(py::self -= py::self)
          .def(py::self *= double())
          .def(py::self == py::self)
          .def(py::self != py::self)
          .def(double() * py::self)
          .def(py::self * double())
          .def("__getitem__",
               [](const DVector& v, const int key) { return v[key]; })
          .def("__len__", [](const DVector& v) { return 3; })
          .def("__str__", py::overload_cast<>(&DVector::toString))
          .def("__repr__", py::overload_cast<>(&DVector::toString))
          .def_static("fromSpherical", &DVector::fromSpherical, "theta"_a, "phi"_a, "r"_a = 1.0);

     py::implicitly_convertible<std::list<double>, DVector>();
     py::implicitly_convertible<std::vector<double>, DVector>();

     py::enum_<Axis>(m, "Axis")
          .value("xaxis", xaxis)
          .value("yaxis", yaxis)
          .value("zaxis", zaxis)
          .value("all", all)
          .value("none", none)
          .export_values();

     py::enum_<Reference>(m, "Reference")
          .value("none", NONE)
          .value("fixed", FIXED)
          .value("top", TOP)
          .value("bottom", BOTTOM)
          .export_values();

     py::enum_<SolverMode>(m, "SolverMode")
          .value("RK4", RK4)
          .value("Heun", HEUN)
          .value("EulerHeun", EULER_HEUN)
          .value("DormandPrince", DORMAND_PRINCE)
          .export_values();

     py::class_<AdaptiveIntegrationParams<double>>(m, "AdaptiveIntegrationParams")
          .def(py::init<>())
          .def_readwrite("abs_tol", &AdaptiveIntegrationParams<double>::abs_tol)
          .def_readwrite("rel_tol", &AdaptiveIntegrationParams<double>::rel_tol)
          .def_readwrite("max_factor", &AdaptiveIntegrationParams<double>::max_factor)
          .def_readwrite("min_factor", &AdaptiveIntegrationParams<double>::min_factor)
          .def_readwrite("safety_factor", &AdaptiveIntegrationParams<double>::safety_factor)
          .def_readwrite("use_pid_control", &AdaptiveIntegrationParams<double>::use_pid_control)
          .def_readwrite("ki", &AdaptiveIntegrationParams<double>::ki)
          .def_readwrite("kp", &AdaptiveIntegrationParams<double>::kp)
          .def_readwrite("kd", &AdaptiveIntegrationParams<double>::kd)
          .def_readwrite("prev_error_ratio", &AdaptiveIntegrationParams<double>::prev_error_ratio)
          .def_readwrite("integral_error", &AdaptiveIntegrationParams<double>::integral_error);

     py::enum_<UpdateType>(m, "UpdateType")
          .value("constant", constant)
          .value("pulse", pulse)
          .value("sine", sine)
          .value("step", step)
          .value("posine", posine)
          .value("halfsine", halfsine)
          .value("trapezoid", trapezoid)
          .value("gaussimpulse", gaussimpulse)
          .value("gaussstep", gaussstep)
          .value("custom", custom)
          .export_values();
     // Driver Class
     py::class_<DScalarDriver>(m, "ScalarDriver")
          .def(py::init<>())
          .def(py::self + double())
          .def(py::self += double())
          .def(py::self * double())
          .def(py::self *= double())
          .def("getCurrentScalarValue", &DScalarDriver::getCurrentScalarValue,
               "time"_a)
          .def_static("getConstantDriver", &DScalarDriver::getConstantDriver,
               "constantValue"_a)
          .def_static("getPulseDriver", &DScalarDriver::getPulseDriver,
               "constantValue"_a, "amplitude"_a, "period"_a, "cycle"_a)
          .def_static("getSineDriver", &DScalarDriver::getSineDriver,
               "constantValue"_a, "amplitude"_a, "frequency"_a, "phase"_a)
          .def_static("getPosSineDriver", &DScalarDriver::getPosSineDriver,
               "constantValue"_a, "amplitude"_a, "frequency"_a, "phase"_a)
          .def_static("getStepDriver", &DScalarDriver::getStepDriver,
               "constantValue"_a, "amplitude"_a, "timeStart"_a, "timeStop"_a)
          .def_static("getTrapezoidDriver", &DScalarDriver::getTrapezoidDriver,
               "constantValue"_a, "amplitude"_a, "timeStart"_a, "edgeTime"_a,
               "steadyTime"_a)
          .def_static("getGaussianImpulseDriver",
               &DScalarDriver::getGaussianImpulseDriver, "constantValue"_a,
               "amplitude"_a, "t0"_a, "sigma"_a)
          .def_static("getGaussianStepDriver",
               &DScalarDriver::getGaussianStepDriver, "constantValue"_a,
               "amplitude"_a, "t0"_a, "sigma"_a)
          .def_static("getCustomDriver", &DScalarDriver::getCustomDriver,
               "callback"_a);

     py::class_<DNullDriver, DScalarDriver>(m, "NullDriver")
          .def(py::init<>())
          .def("getCurrentScalarValue", &DScalarDriver::getCurrentScalarValue,
               "time"_a);

     py::class_<DAxialDriver>(m, "AxialDriver")
          .def(py::init<DScalarDriver, DScalarDriver, DScalarDriver>())
          .def(py::init<std::vector<ScalarDriver<double>>>())
          .def(py::init<double, double, double>())
          .def(py::init<DVector>())
          .def("getVectorAxialDriver", &DAxialDriver::getVectorAxialDriver)
          .def("getCurrentAxialDrivers", &DAxialDriver::getCurrentAxialDrivers,
               "time"_a)
          .def("applyMask",
               py::overload_cast<const DVector&>(&DAxialDriver::applyMask))
          .def("applyMask", py::overload_cast<const std::vector<unsigned int> &>(
               &DAxialDriver::applyMask));

     py::class_<DLayer>(m, "Layer")
          .def(py::init<std::string,          // id
               DVector,              // mag
               DVector,              // anis
               double,               // Ms
               double,               // thickness
               double,               // cellSurface
               std::vector<DVector>, // demagTensor
               double                // damping
          >(),
               "id"_a, "mag"_a, "anis"_a, "Ms"_a, "thickness"_a, "cellSurface"_a,
               "demagTensor"_a, "damping"_a = 0.011)
          .def_static("createSOTLayer", &DLayer::LayerSOT, "id"_a, "mag"_a,
               "anis"_a, "Ms"_a, "thickness"_a, "cellSurface"_a,
               "demagTensor"_a, "damping"_a = 0.011,
               "fieldLikeTorque"_a = 0.0, "dampingLikeTorque"_a = 0.0)
          .def_static("createSTTLayer", &DLayer::LayerSTT, "id"_a, "mag"_a,
               "anis"_a, "Ms"_a, "thickness"_a, "cellSurface"_a,
               "demagTensor"_a, "damping"_a = 0.011,
               "SlonczewskiSpacerLayerParameter"_a = 1.0, "beta"_a = 0.0,
               "spinPolarisation"_a = 0.0)
          .def("setMagnetisation", &DLayer::setMagnetisation)
          .def("setAnisotropyDriver", &DLayer::setAnisotropyDriver)
          .def("setSecondOrderAnisotropyDriver", &DLayer::setSecondOrderAnisotropyDriver)
          .def("setExternalFieldDriver", &DLayer::setExternalFieldDriver)
          .def("setOerstedFieldDriver", &DLayer::setOerstedFieldDriver)
          .def("setHdmiDriver", &DLayer::setHdmiDriver)
          // reference layers
          .def("setReferenceLayer",
               py::overload_cast<const DVector&>(&DLayer::setReferenceLayer))
          .def("setReferenceLayer",
               py::overload_cast<Reference>(&DLayer::setReferenceLayer))
          .def("setSecondaryReferenceLayer",
               &DLayer::setSecondaryReferenceLayer)
          // drivers
          .def("setFieldLikeTorqueDriver", &DLayer::setFieldLikeTorqueDriver)
          .def("setDampingLikeTorqueDriver", &DLayer::setDampingLikeTorqueDriver)
          .def("setSecondaryFieldLikeTorqueDriver", &DLayer::setSecondaryFieldLikeTorqueDriver)
          .def("setSecondaryDampingLikeTorqueDriver", &DLayer::setSecondaryDampingLikeTorqueDriver)
          .def("setPrimaryTorqueDrivers", &DLayer::setPrimaryTorqueDrivers, "fieldLikeTorque"_a, "dampingLikeTorque"_a)
          .def("setSecondaryTorqueDrivers", &DLayer::setSecondaryTorqueDrivers, "fieldLikeTorque"_a, "dampingLikeTorque"_a)

          .def("setTemperatureDriver", &DLayer::setTemperatureDriver)
          .def("setTopDipoleTensor", &DLayer::setTopDipoleTensor)
          .def("setBottomDipoleTensor", &DLayer::setBottomDipoleTensor)
          .def("setKappa", &DLayer::setKappa)
          .def("setAlternativeSTT", &DLayer::setAlternativeSTT)
          // readonly props
          .def_readonly("id", &DLayer::id)
          .def_readonly("Ms", &DLayer::Ms)
          .def_readonly("thickness", &DLayer::thickness)
          .def_readonly("damping", &DLayer::damping)
          .def_readonly("cellSurface", &DLayer::cellSurface)
          .def_readonly("demagTensor", &DLayer::demagTensor)
          // noise
          .def("setAlphaNoise", &DLayer::setAlphaNoise, "alpha"_a, "std"_a, "scale"_a, "axis"_a = Axis::all)
          .def("setOneFNoise", &DLayer::setOneFNoise)
          // getters
          .def("getReferenceLayer", &DLayer::getReferenceLayer, py::return_value_policy::reference)
          .def("getSecondaryReferenceLayer", &DLayer::getSecondaryReferenceLayer, py::return_value_policy::reference)
          .def("getId", &DLayer::getId)
          .def("getOneFVector", &DLayer::getOneFVector)
          .def("setAdaptiveParams", &DLayer::setAdaptiveParams, "params"_a)
          .def("createBufferedAlphaNoise", &DLayer::createBufferedAlphaNoise);

     py::class_<DJunction>(m, "Junction")
          .def(py::init<std::vector<DLayer>>(), "layers"_a)
          .def(py::init<std::vector<DLayer>, double, double>(), "layers"_a,
               "Rp"_a = 100, "Rap"_a = 200)
          .def(py::init<std::vector<DLayer>, std::vector<double>,
               std::vector<double>, std::vector<double>,
               std::vector<double>, std::vector<double>,
               std::vector<double>, std::vector<double>>(),
               "layers"_a, "Rx0"_a, "Ry0"_a, "AMR_X"_a, "AMR_Y"_a, "SMR_X"_a,
               "SMR_Y"_a, "AHE"_a)
          // log utils
          .def("getLog", &DJunction::getLog)
          .def("clearLog", &DJunction::clearLog)
          .def("saveLog", &DJunction::saveLogs, "filename"_a)
          // main run
          .def("runSimulation", &DJunction::runSimulation, "totalTime"_a,
               "timeStep"_a = 1e-13, "writeFrequency"_a = 1e-11, "verbose"_a = false,
               "calculateEnergies"_a = false, "solverMode"_a = RK4)

          // driver setters
          .def("setLayerOerstedFieldDriver", &DJunction::setLayerOerstedFieldDriver)
          .def("setLayerExternalFieldDriver",
               &DJunction::setLayerExternalFieldDriver)
          .def("setLayerCurrentDriver", &DJunction::setLayerCurrentDriver)
          .def("setLayerAnisotropyDriver", &DJunction::setLayerAnisotropyDriver)
          .def("setLayerSecondOrderAnisotropyDriver", &DJunction::setLayerSecondOrderAnisotropyDriver)
          .def("setLayerOerstedFieldDriver", &DJunction::setLayerOerstedFieldDriver)
          .def("setLayerMagnetisation", &DJunction::setLayerMagnetisation)
          .def("setLayerHdmiDriver", &DJunction::setLayerHdmiDriver)
          // interaction setters
          .def("setIECDriver", &DJunction::setIECDriver)
          .def("setQuadIECDriver", &DJunction::setQuadIECDriver)
          .def("setIDMIDriver", &DJunction::setIDMIDriver)
          // noise
          .def("setLayerTemperatureDriver", &DJunction::setLayerTemperatureDriver)
          .def("setLayerNonStochasticLangevinDriver",
               &DJunction::setLayerNonStochasticLangevinDriver)
          .def("setLayerOneFNoise", &DJunction::setLayerOneFNoise)
          // SOT setters
          .def("setLayerFieldLikeTorqueDriver",
               &DJunction::setLayerFieldLikeTorqueDriver)
          .def("setLayerDampingLikeTorqueDriver",
               &DJunction::setLayerDampingLikeTorqueDriver)
          .def("setLayerSecondaryFieldLikeTorqueDriver",
               &DJunction::setLayerSecondaryFieldLikeTorqueDriver)
          .def("setLayerSecondaryDampingLikeTorqueDriver",
               &DJunction::setLayerSecondaryDampingLikeTorqueDriver)
          .def("setLayerPrimaryTorqueDrivers",
               &DJunction::setLayerPrimaryTorqueDrivers, "layerId"_a, "fieldLikeTorque"_a, "dampingLikeTorque"_a)
          .def("setLayerSecondaryTorqueDrivers",
               &DJunction::setLayerSecondaryTorqueDrivers, "layerId"_a, "fieldLikeTorque"_a, "dampingLikeTorque"_a)
          // Reference setters
          .def("setLayerReferenceType", &DJunction::setLayerReferenceType)
          .def("setLayerReferenceLayer", &DJunction::setLayerReferenceLayer)
          // other setters
          .def("setLayerAlternativeSTT", &DJunction::setLayerAlternativeSTT)
          // junction calculations
          .def("getLayerMagnetisation", &DJunction::getLayerMagnetisation)
          .def("getMagnetoresistance", &DJunction::getMagnetoresistance)
          // getters
          .def("getLayerIds", &DJunction::getLayerIds)
          .def("getLayer", &DJunction::getLayer, "layerId"_a,
               py::return_value_policy::reference)
          // readonly props
          .def_readonly("layers", &DJunction::layers);

     // stack module
     py::module stack_module =
          m.def_submodule("stack", "A stack submodule for joining MTJ junctions");

     py::class_<DSeriesStack>(stack_module, "SeriesStack")
          .def(py::init<std::vector<DJunction>, std::string, std::string, double, bool>(),
               "junctionList"_a, "topId_a"_a = "free", "bottomId"_a = "bottom",
               "phaseOffset"_a = 0.0, "useKCL"_a = true)
          .def("runSimulation", &DSeriesStack::runSimulation, "totalTime"_a,
               "timeStep"_a = 1e-13, "writeFrequency"_a = 1e-11)
          .def("setMagnetisation", &DSeriesStack::setMagnetisation, "junction"_a,
               "layerId"_a, "mag"_a)
          .def("getMagnetisation", &DSeriesStack::getMagnetisation, "junction"_a,
               "layerId"_a)
          .def("setCoupledCurrentDriver", &DSeriesStack::setCoupledCurrentDriver,
               "driver"_a)
          .def("setExternalFieldDriver", &DSeriesStack::setExternalFieldDriver,
               "driver"_a)
          .def(
               "setCouplingStrength",
               py::overload_cast<const double&>(&DSeriesStack::setCouplingStrength),
               "coupling"_a)
          .def("setCouplingStrength",
               py::overload_cast<const std::vector<double> &>(
                    &DSeriesStack::setCouplingStrength),
               "coupling"_a)
          .def("setDelayed", &DSeriesStack::setDelayed, "delayed"_a)
          .def("getJunction", &DParallelStack::getJunction, "junctionId"_a,
               py::return_value_policy::reference)
          .def("setJunctionAnisotropyDriver",
               &DSeriesStack::setJunctionAnisotropyDriver, "junctionId"_a,
               "layerId"_a, "k"_a)
          // logging
          .def("clearLogs", &DSeriesStack::clearLogs)
          .def("getLog", py::overload_cast<unsigned int>(&DSeriesStack::getLog))
          .def("getLog", py::overload_cast<>(&DSeriesStack::getLog));

     py::class_<DParallelStack>(stack_module, "ParallelStack")
          .def(py::init<std::vector<DJunction>, std::string, std::string, double, bool>(),
               "junctionList"_a, "topId_a"_a = "free", "bottomId"_a = "bottom",
               "phaseOffset"_a = 0.0, "useKCL"_a = true)
          .def("runSimulation", &DParallelStack::runSimulation, "totalTime"_a,
               "timeStep"_a = 1e-13, "writeFrequency"_a = 1e-11)
          .def("setMagnetisation", &DParallelStack::setMagnetisation, "junction"_a,
               "layerId"_a, "mag"_a)
          .def("getMagnetisation", &DParallelStack::getMagnetisation, "junction"_a,
               "layerId"_a)
          .def("setCoupledCurrentDriver", &DParallelStack::setCoupledCurrentDriver,
               "driver"_a)
          .def("setExternalFieldDriver", &DParallelStack::setExternalFieldDriver,
               "driver"_a)
          .def("setCouplingStrength",
               py::overload_cast<const double&>(
                    &DParallelStack::setCouplingStrength),
               "coupling"_a)
          .def("setCouplingStrength",
               py::overload_cast<const std::vector<double> &>(
                    &DParallelStack::setCouplingStrength),
               "coupling"_a)
          .def("setDelayed", &DParallelStack::setDelayed, "delayed"_a)
          .def("getJunction", &DParallelStack::getJunction, "junctionId"_a,
               py::return_value_policy::reference)
          .def("setJunctionAnisotropyDriver",
               &DSeriesStack::setJunctionAnisotropyDriver, "junctionId"_a,
               "layerId"_a, "k"_a)
          // logging
          .def("clearLogs", &ParallelStack<double>::clearLogs)
          .def("getLog",
               py::overload_cast<unsigned int>(&ParallelStack<double>::getLog))
          .def("getLog", py::overload_cast<>(&ParallelStack<double>::getLog));

     // reservoir module
     py::module reservoir_module = m.def_submodule(
          "reservoir", "A reservoir submodule for joining MTJ junctions");
     reservoir_module.def("nullDipoleInteraction", &nullDipoleInteraction, "r1"_a, "r2"_a, "layer1"_a, "layer2"_a);
     reservoir_module.def("computeDipoleInteraction", &computeDipoleInteraction, "r1"_a, "r2"_a, "layer1"_a, "layer2"_a);
     reservoir_module.def("computeDipoleInteractionNoumra", &computeDipoleInteractionNoumra, "r1"_a, "r2"_a, "layer1"_a, "layer2"_a);

     py::class_<GroupInteraction>(reservoir_module, "GroupInteraction")
          .def(py::init<const std::vector<DVector>&, const std::vector<DJunction>&, const std::string&>(),
               "coordinateMatrix"_a, "junctionList"_a, "topId"_a = "free")
          .def("setInteractionFunction", &GroupInteraction::setInteractionFunction)
          .def("runSimulation", &GroupInteraction::runSimulation, "totalTime"_a,
               "timeStep"_a = 1e-13, "writeFrequency"_a = 1e-11)
          .def("clearLogs", &GroupInteraction::clearLogs)
          .def("getLog",
               py::overload_cast<unsigned int>(&GroupInteraction::getLog))
          .def("getLog", py::overload_cast<unsigned int>(&GroupInteraction::getLog),
               py::return_value_policy::reference);


     // generator module
     py::module generator_module =
          m.def_submodule("noise", "Submodule with noise generation functions");
     py::class_<BufferedAlphaNoise<double>>(generator_module, "BufferedAlphaNoise")
          .def(py::init<unsigned int, double, double, double>(), "bufferSize"_a,
               "alpha"_a, "std"_a, "scale"_a)
          .def("fillBuffer", &BufferedAlphaNoise<double>::fillBuffer)
          .def("tick", &BufferedAlphaNoise<double>::tick);
     py::class_<VectorAlphaNoise<double>>(generator_module, "VectorAlphaNoise")
          .def(py::init<unsigned int, double, double, double, Axis>(),
               "bufferSize"_a, "alpha"_a, "std"_a, "scale"_a, "axis"_a = Axis::all)
          .def("tickVector", &VectorAlphaNoise<double>::tickVector)
          .def("tick", &VectorAlphaNoise<double>::tick)
          .def("getPrevSample", &VectorAlphaNoise<double>::getPrevSample)
          .def("getScale", &VectorAlphaNoise<double>::getScale);

     // LLGB module
     auto llgb_module = m.def_submodule("llgb", "A submodule for LLGB junctions");
     llgb_module.def("MFAWeissCurie", &LLGB::MFAWeissCurie<double>, "me"_a, "T"_a,
          "J0"_a, "relax"_a = 0.2, "tolerance"_a = 1e-6,
          "maxIter"_a = 1000);
     llgb_module.def("langevin", &LLGB::langevin<double>);
     llgb_module.def("langevinDerivative", &LLGB::langevinDerivative<double>);

     py::class_<DLLGBLayer>(llgb_module, "LLGBLayer")
          .def(py::init<const std::string&, DVector, DVector, double, double,
               double, const std::vector<DVector> &, double, double,
               double, double>(),
               "id"_a, "mag"_a, "anis"_a, "Ms"_a, "thickness"_a, "cellSurface"_a,
               "demagTensor"_a, "damping"_a, "Tc"_a, "susceptibility"_a, "me"_a)
          // setters
          .def("setTemperatureDriver", &DLLGBLayer::setTemperatureDriver)
          .def("setExternalFieldDriver", &DLLGBLayer::setExternalFieldDriver)
          .def("setAnisotropyDriver", &DLLGBLayer::setAnisotropyDriver);

     py::class_<DLLGBJunction>(llgb_module, "LLGBJunction")
          .def(py::init<std::vector<DLLGBLayer>>(), "layers"_a)
          .def("runSimulation", &DLLGBJunction::runSimulation, "totalTime"_a,
               "timeStep"_a = 1e-13, "writeFrequency"_a = 1e-11, "verbose"_a = false,
               "solverMode"_a = HEUN)
          .def("setLayerTemperatureDriver",
               &DLLGBJunction::setLayerTemperatureDriver)
          .def("setLayerExternalFieldDriver",
               &DLLGBJunction::setLayerExternalFieldDriver)
          .def("saveLogs", &DLLGBJunction::saveLogs)
          .def("getLog", &DLLGBJunction::getLog)
          .def("clearLog", &DLLGBJunction::clearLog);

     // constants module
     py::module constants_module = m.def_submodule("constants", "A submodule for physical constants");
     py::class_<PhysicalConstants>(constants_module, "PhysicalConstants")
          .def_static("set_magnetic_permeability", &PhysicalConstants::setMagneticPermeability, "value"_a)
          .def_static("set_gyromagnetic_ratio", &PhysicalConstants::setGyro, "value"_a)
          .def_static("set_TtoAm", &PhysicalConstants::setTtoAm, "value"_a)
          .def_static("set_hbar", &PhysicalConstants::setHbar, "value"_a)
          .def_static("set_elementary_charge", &PhysicalConstants::setElectronCharge, "value"_a)
          .def_static("set_boltzmann_constant", &PhysicalConstants::setBoltzmannConst, "value"_a)
          .def_static("magnetic_permeability", &PhysicalConstants::getMagneticPermeability)
          .def_static("gyromagnetic_ratio", &PhysicalConstants::getGyro)
          .def_static("TtoAm", &PhysicalConstants::getTtoAm)
          .def_static("hbar", &PhysicalConstants::getHbar)
          .def_static("elementary_charge", &PhysicalConstants::getElectronCharge)
          .def_static("boltzmann_constant", &PhysicalConstants::getBoltzmannConst)
          .def_static("resetToDefaults", &PhysicalConstants::resetToDefaults);

}
