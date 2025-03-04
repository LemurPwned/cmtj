#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/2.0/driver.hpp"
#include "../core/2.0/fm.hpp"
#include "../core/2.0/junction.hpp"
#include "../core/cvector.hpp"
#include <stdio.h>
#include <vector>

using namespace pybind11::literals;
namespace py = pybind11;

using DJunction = FMJunction<double>;
using DLayer = Layer<double>;
using DVector = CVector<double>;
using DScalarDriver = ScalarDriver<double>;
using DAxialDriver = AxialDriver<double>;
using DNullDriver = NullDriver<double>;

#define USING_PY true
PYBIND11_MODULE(pycmtj, m) {
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
      "stepDriver",
      [](double constantValue, double amplitude, double timeStart,
         double timeStop) {
        return DScalarDriver::getStepDriver(constantValue, amplitude, timeStart,
                                            timeStop);
      },
      "constantValue"_a, "amplitude"_a, "timeStart"_a, "timeStop"_a);
  // CVector
  py::class_<DVector>(m, "CVector")
      .def(py::init<double, double, double>())
      .def("__getitem__",
           [](const DVector &v, size_t i) -> double {
             if (i >= 3)
               throw py::index_error();
             return v[i];
           })
      .def("__setitem__",
           [](DVector &v, size_t i, double val) {
             if (i >= 3)
               throw py::index_error();
             if (i == 0)
               v.x = val;
             else if (i == 1)
               v.y = val;
             else
               v.z = val;
           })
      .def("normalize", &DVector::normalize)
      .def("length", py::overload_cast<>(&DVector::length))
      .def("tolist", &DVector::tolist)
      // operators
      .def(py::self + py::self)
      .def(py::self += py::self)
      .def(py::self - py::self)
      .def(py::self -= py::self)
      .def(py::self * double())
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(double() * py::self)
      .def(py::self * double())
      .def("__len__", [](const DVector &) { return 3; })
      .def("__str__", py::overload_cast<>(&DVector::toString))
      .def("__repr__", py::overload_cast<>(&DVector::toString))
      .def_static("fromSpherical", &DVector::fromSpherical, "theta"_a, "phi"_a,
                  "r"_a = 1.0);

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
      .value("DormandPrice", DORMAND_PRICE)
      .export_values();

  // Driver Class
  py::class_<DScalarDriver>(m, "ScalarDriver")
      .def(py::init<>())
      .def("getCurrentScalarValue", &DScalarDriver::getCurrentScalarValue,
           "time"_a)
      .def_static("getConstantDriver", &DScalarDriver::getConstantDriver,
                  "constantValue"_a)
      .def_static("getPulseDriver", &DScalarDriver::getPulseDriver,
                  "constantValue"_a, "amplitude"_a, "period"_a, "cycle"_a)
      .def_static("getSineDriver", &DScalarDriver::getSineDriver,
                  "constantValue"_a, "amplitude"_a, "frequency"_a, "phase"_a)
      .def_static("getStepDriver", &DScalarDriver::getStepDriver,
                  "constantValue"_a, "amplitude"_a, "timeStart"_a,
                  "timeStop"_a);

  py::class_<DNullDriver, std::shared_ptr<DNullDriver>>(m, "NullDriver")
      .def(py::init<double>())
      .def("getCurrentScalarValue", &DNullDriver::getCurrentScalarValue,
           "time"_a);

  py::class_<DAxialDriver>(m, "AxialDriver")
      .def(py::init<DVector>())
      .def(py::init<std::shared_ptr<Driver<double>>,
                    std::shared_ptr<Driver<double>>,
                    std::shared_ptr<Driver<double>>>())
      .def("getCurrentAxialDrivers", &DAxialDriver::getCurrentAxialDrivers,
           "time"_a);

  // AbstractLayer class (base class)
  py::class_<AbstractLayer<double>, std::shared_ptr<AbstractLayer<double>>>(
      m, "AbstractLayer")
      .def("getMagnetisation", &AbstractLayer<double>::getMagnetisation)
      .def("setMagnetisation", &AbstractLayer<double>::setMagnetisation)
      .def("setReferenceLayer", &AbstractLayer<double>::setReferenceLayer)
      .def("setReferenceType", &AbstractLayer<double>::setReferenceType)
      .def("setAnisotropyDriver", &AbstractLayer<double>::setAnisotropyDriver)
      .def("setTemperatureDriver", &AbstractLayer<double>::setTemperatureDriver)
      .def("setCurrentDriver", &AbstractLayer<double>::setCurrentDriver)
      .def("setExternalFieldDriver",
           &AbstractLayer<double>::setExternalFieldDriver)
      .def("setHdmiDriver", &AbstractLayer<double>::setHdmiDriver)
      .def("setReservedInteractionFieldDriver",
           &AbstractLayer<double>::setReservedInteractionFieldDriver)
      .def("getBufferedNoiseParameters",
           &AbstractLayer<double>::getBufferedNoiseParameters)
      .def("getId", &AbstractLayer<double>::getId);

  // Layer class (standard FM layer)
  py::class_<Layer<double>, AbstractLayer<double>,
             std::shared_ptr<Layer<double>>>(m, "Layer")
      .def(py::init<const std::string &, const DVector &, const DVector &,
                    double, double, double, const std::vector<DVector> &,
                    double>())
      .def("setMagnetisation", &Layer<double>::setMagnetisation)
      .def("getMagnetisation", &Layer<double>::getMagnetisation)
      .def("setReferenceLayer", &Layer<double>::setReferenceLayer)
      .def("setReferenceType", &Layer<double>::setReferenceType)
      .def("setTemperatureDriver", &Layer<double>::setTemperatureDriver)
      .def("setCurrentDriver", &Layer<double>::setCurrentDriver)
      .def("setAnisotropyDriver", &Layer<double>::setAnisotropyDriver)
      .def("setExternalFieldDriver", &Layer<double>::setExternalFieldDriver)
      .def("setIECDriverBottom", &Layer<double>::setIECDriverBottom)
      .def("setIECDriverTop", &Layer<double>::setIECDriverTop)
      .def("setQuadIECDriverTop", &Layer<double>::setQuadIECDriverTop)
      .def("setQuadIECDriverBottom", &Layer<double>::setQuadIECDriverBottom)
      .def("setIDMIDriver", &Layer<double>::setIDMIDriver);

  // LayerSTT class (for Spin Transfer Torque)
  py::class_<LayerSTT<double>, Layer<double>,
             std::shared_ptr<LayerSTT<double>>>(m, "LayerSTT")
      .def(py::init<const std::string &, const DVector &, const DVector &,
                    double, double, double, const std::vector<DVector> &,
                    double, double, double, double>())
      .def_readwrite("SlonczewskiParameter",
                     &LayerSTT<double>::SlonczewskiParameter)
      .def_readwrite("beta", &LayerSTT<double>::beta)
      .def_readwrite("spinPolarisation", &LayerSTT<double>::spinPolarisation)
      .def_readwrite("kappa", &LayerSTT<double>::kappa)
      .def_readwrite("alternativeSTTSet", &LayerSTT<double>::alternativeSTTSet);

  // LayerSOT class (for Spin Orbit Torque)
  py::class_<LayerSOT<double>, Layer<double>,
             std::shared_ptr<LayerSOT<double>>>(m, "LayerSOT")
      .def(py::init<const std::string &, const DVector &, const DVector &,
                    double, double, double, const std::vector<DVector> &,
                    double, double, double, bool>())
      .def("setFieldLikeTorqueDriver",
           &LayerSOT<double>::setFieldLikeTorqueDriver)
      .def("setDampingLikeTorqueDriver",
           &LayerSOT<double>::setDampingLikeTorqueDriver)
      .def_readwrite("directTorques", &LayerSOT<double>::directTorques);

  // AbstractJunction (base class)
  py::class_<AbstractJunction<double>,
             std::shared_ptr<AbstractJunction<double>>>(m, "AbstractJunction")
      .def("getLayerIds", &AbstractJunction<double>::getLayerIds)
      .def("getLayer", &AbstractJunction<double>::getLayer,
           py::return_value_policy::reference)
      .def("runSimulation", &AbstractJunction<double>::runSimulation,
           "totalTime"_a, "timeStep"_a = 1e-13, "writeFrequency"_a = 1e-13,
           "mode"_a = HEUN)
      .def("saveLogs", &AbstractJunction<double>::saveLogs)
      .def("getLog", &AbstractJunction<double>::getLog)
      .def("clearLog", &AbstractJunction<double>::clearLog);

  // FMJunction class
  py::class_<FMJunction<double>, AbstractJunction<double>,
             std::shared_ptr<FMJunction<double>>>(m, "FMJunction")
      .def(py::init<
           const std::vector<std::shared_ptr<AbstractLayer<double>>> &>())
      .def("setLayerAnisotropyDriver",
           &FMJunction<double>::setLayerAnisotropyDriver)
      .def("setLayerTemperatureDriver",
           &FMJunction<double>::setLayerTemperatureDriver)
      .def("setLayerCurrentDriver", &FMJunction<double>::setLayerCurrentDriver)
      .def("setLayerExternalFieldDriver",
           &FMJunction<double>::setLayerExternalFieldDriver)
      .def("setMRMode", &FMJunction<double>::setMRMode)
      .def("setRp", &FMJunction<double>::setRp)
      .def("setRap", &FMJunction<double>::setRap)
      .def("getMagnetoresistance", &FMJunction<double>::getMagnetoresistance)
      .def("getLayer", &FMJunction<double>::getLayer, "layerId"_a,
           py::return_value_policy::reference);
}
