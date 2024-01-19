#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/reservoir.hpp"
#include "../core/stack.hpp"
#include "../core/cvector.hpp"
#include "../core/drivers.hpp"
#include "../core/junction.hpp"
#include "../core/noise.hpp"
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
        .def("__getitem__", [](const DVector& v, const int key) { return v[key]; })
        .def("__len__", [](const DVector& v) { return 3; })
        .def("__str__", py::overload_cast<>(&DVector::toString))
        .def("__repr__", py::overload_cast<>(&DVector::toString));

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
        .value("Heun", HEUN)
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
            "timeStop"_a)
        .def_static("getTrapezoidDriver",
            &DScalarDriver::getTrapezoidDriver,
            "constantValue"_a,
            "amplitude"_a,
            "timeStart"_a,
            "edgeTime"_a,
            "steadyTime"_a)
        .def_static("getGaussianImpulseDriver",
            &DScalarDriver::getGaussianImpulseDriver,
            "constantValue"_a,
            "amplitude"_a,
            "t0"_a,
            "sigma"_a)
        .def_static("getGaussianStepDriver",
            &DScalarDriver::getGaussianStepDriver,
            "constantValue"_a,
            "amplitude"_a,
            "t0"_a,
            "sigma"_a);

    py::class_<DNullDriver, DScalarDriver>(m, "NullDriver")
        .def(py::init<>());

    py::class_<DAxialDriver>(m, "AxialDriver")
        .def(py::init<DScalarDriver, DScalarDriver, DScalarDriver>())
        .def(py::init<std::vector<ScalarDriver<double>>>())
        .def(py::init<double, double, double>())
        .def(py::init<DVector>())
        .def("getVectorAxialDriver", &DAxialDriver::getVectorAxialDriver)
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
        .def("setReferenceLayer", py::overload_cast<const DVector&>(&DLayer::setReferenceLayer))
        .def("setReferenceLayer", py::overload_cast<Reference>(&DLayer::setReferenceLayer))

        .def("setFieldLikeTorqueDriver", &DLayer::setFieldLikeTorqueDriver)
        .def("setDampingLikeTorqueDriver", &DLayer::setDampingLikeTorqueDriver)
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
        .def("setAlphaNoise", &DLayer::setAlphaNoise)
        .def("setOneFNoise", &DLayer::setOneFNoise)
        // getters
        .def("getId", &DLayer::getId)
        .def("getOneFVector", &DLayer::getOneFVector)
        .def("createBufferedAlphaNoise", &DLayer::createBufferedAlphaNoise);

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
        // main run
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
        // noise
        .def("setLayerTemperatureDriver", &DJunction::setLayerTemperatureDriver)
        .def("setLayerNonStochasticLangevinDriver", &DJunction::setLayerNonStochasticLangevinDriver)
        .def("setLayerOneFNoise", &DJunction::setLayerOneFNoise)
        // SOT setters
        .def("setLayerFieldLikeTorqueDriver", &DJunction::setLayerFieldLikeTorqueDriver)
        .def("setLayerDampingLikeTorqueDriver", &DJunction::setLayerDampingLikeTorqueDriver)
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
        // readonly props
        .def_readonly("layers", &DJunction::layers);

    // stack module
    py::module stack_module = m.def_submodule("stack", "A stack submodule for joining MTJ junctions");
    py::class_<DSeriesStack>(stack_module, "SeriesStack")
        .def(py::init<std::vector<DJunction>,
            std::string,
            std::string>(),
            "junctionList"_a,
            "topId_a"_a = "free",
            "bottomId"_a = "bottom")
        .def("runSimulation", &DSeriesStack::runSimulation,
            "totalTime"_a,
            "timeStep"_a = 1e-13,
            "writeFrequency"_a = 1e-11)
        .def("setMagnetisation", &DSeriesStack::setMagnetisation, "junction"_a, "layerId"_a, "mag"_a)
        .def("getMagnetisation", &DSeriesStack::getMagnetisation, "junction"_a, "layerId"_a)
        .def("setCoupledCurrentDriver", &DSeriesStack::setCoupledCurrentDriver, "driver"_a)
        .def("setExternalFieldDriver", &DSeriesStack::setExternalFieldDriver, "driver"_a)
        .def("setCouplingStrength", &DSeriesStack::setCouplingStrength, "coupling"_a)
        .def("setDelayed", &DParallelStack::setDelayed, "delayed"_a)
        // logging
        .def("clearLogs", &DSeriesStack::clearLogs)
        .def("getLog", py::overload_cast<unsigned int>(&DSeriesStack::getLog))
        .def("getLog", py::overload_cast<>(&DSeriesStack::getLog));

    py::class_<DParallelStack>(stack_module, "ParallelStack")
        .def(py::init<std::vector<DJunction>,
            std::string,
            std::string>(),
            "junctionList"_a,
            "topId_a"_a = "free",
            "bottomId"_a = "bottom")
        .def("runSimulation", &DParallelStack::runSimulation,
            "totalTime"_a,
            "timeStep"_a = 1e-13,
            "writeFrequency"_a = 1e-11)
        .def("setMagnetisation", &DParallelStack::setMagnetisation, "junction"_a, "layerId"_a, "mag"_a)
        .def("getMagnetisation", &DParallelStack::getMagnetisation, "junction"_a, "layerId"_a)
        .def("setCoupledCurrentDriver", &DParallelStack::setCoupledCurrentDriver, "driver"_a)
        .def("setExternalFieldDriver", &DParallelStack::setExternalFieldDriver, "driver"_a)
        .def("setCouplingStrength", &DParallelStack::setCouplingStrength, "coupling"_a)
        .def("setDelayed", &DParallelStack::setDelayed, "delayed"_a)
        // logging
        .def("clearLogs", &ParallelStack<double>::clearLogs)
        .def("getLog", py::overload_cast<unsigned int>(&ParallelStack<double>::getLog))
        .def("getLog", py::overload_cast<>(&ParallelStack<double>::getLog));


    // reservoir module
    py::module reservoir_module = m.def_submodule("reservoir", "A reservoir submodule for joining MTJ junctions");
    py::class_<Reservoir>(reservoir_module, "Reservoir")
        .def(py::init<DVectorMatrix, DLayerMatrix >(),
            "coordinateMatrix"_a,
            "layerMatrix"_a)
        .def("runSimulation", &Reservoir::runSimulation)
        .def("clearLogs", &Reservoir::clearLogs)
        .def("saveLogs", &Reservoir::saveLogs)
        .def("getLayer", &Reservoir::getLayer)
        .def("setAllExternalField", &Reservoir::setAllExternalField)
        .def("setLayerAnisotropy", &Reservoir::setLayerAnisotropy)
        .def("setLayerExternalField", &Reservoir::setLayerExternalField)
        .def("getMagnetisation", &Reservoir::getMagnetisation);

    // generator module
    py::module generator_module = m.def_submodule("noise", "Submodule with noise generation functions");
    py::class_<BufferedAlphaNoise<double>>(generator_module, "BufferedAlphaNoise")
        .def(py::init<unsigned int, double, double, double>(),
            "bufferSize"_a,
            "alpha"_a,
            "std"_a,
            "scale"_a)
        .def("fillBuffer", &BufferedAlphaNoise<double>::fillBuffer)
        .def("tick", &BufferedAlphaNoise<double>::tick);
    py::class_<VectorAlphaNoise<double>>(generator_module, "VectorAlphaNoise")
        .def(py::init<unsigned int, double, double, double>(),
            "bufferSize"_a,
            "alpha"_a,
            "std"_a,
            "scale"_a)
        .def("tickVector", &VectorAlphaNoise<double>::tickVector)
        .def("tick", &VectorAlphaNoise<double>::tick)
        .def("getPrevSample", &VectorAlphaNoise<double>::getPrevSample)
        .def("getScale", &VectorAlphaNoise<double>::getScale);
}
