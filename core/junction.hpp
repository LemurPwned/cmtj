/**
 * @file junction.hpp
 * @author Jakub
 * @brief CMTJ Junction model
 * @version 1.0
 * @date 2022-03-22
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef CORE_JUNCTION_HPP_
#define CORE_JUNCTION_HPP_

#define _USE_MATH_DEFINES
#include "constants.hpp"
#include "cvector.hpp"   // for CVector
#include "drivers.hpp"   // for ScalarDriver, AxialDriver
#include "noise.hpp"     // for OneFNoise
#include <algorithm>     // for find_if
#include <array>         // for array, array<>::value_type
#include <chrono>        // for seconds, steady_clock, duration
#include <cmath>         // for isnan, M_PI
#include <fstream>       // for file save
#include <functional>    // for bind, function
#include <iostream>      // for string, operator<<, basic_ostream
#include <random>        // for mt19937, normal_distribution
#include <stdexcept>     // for runtime_error, invalid_argument
#include <string>        // for operator+, operator==, basic_string
#include <type_traits>   // for enable_if<>::type
#include <unordered_map> // for unordered_map
#include <vector>        // for vector, __vector_base<>::value_type

#define MAGNETIC_PERMEABILITY PhysicalConstants::MAGNETIC_PERMEABILITY
#define GYRO PhysicalConstants::GYRO
#define TtoAm PhysicalConstants::TtoAm
#define HBAR PhysicalConstants::HBAR
#define ELECTRON_CHARGE PhysicalConstants::ELECTRON_CHARGE
#define BOLTZMANN_CONST PhysicalConstants::BOLTZMANN_CONST

typedef CVector<double> DVector;
typedef CVector<float> FVector;

double operator"" _ns(unsigned long long timeUnit) {
  return ((double)timeUnit) / 1e9;
}
double operator"" _ns(long double timeUnit) { return ((double)timeUnit) / 1e9; }

double operator"" _mT(unsigned long long tesla) {
  return ((double)tesla) / 1000.0;
}

double operator"" _mT(long double tesla) { return ((double)tesla) / 1000.0; }

template <typename T>
inline CVector<T> calculate_tensor_interaction(
    const CVector<T> &m, const std::vector<CVector<T>> &tensor, const T &Ms) {
  // Cache magnetization components for better register usage
  const T m0 = m[0], m1 = m[1], m2 = m[2];
  const T scale = Ms / MAGNETIC_PERMEABILITY;
  
  // Compute matrix-vector product with fewer temporary objects
  return CVector<T>(
      (tensor[0][0] * m0 + tensor[0][1] * m1 + tensor[0][2] * m2) * scale,
      (tensor[1][0] * m0 + tensor[1][1] * m1 + tensor[1][2] * m2) * scale,
      (tensor[2][0] * m0 + tensor[2][1] * m1 + tensor[2][2] * m2) * scale);
}

template <typename T>
inline CVector<T> calculate_tensor_interaction(
    const CVector<T> &m, const std::array<CVector<T>, 3> &tensor, const T &Ms) {
  // Cache magnetization components for better register usage
  const T m0 = m[0], m1 = m[1], m2 = m[2];
  const T scale = Ms / MAGNETIC_PERMEABILITY;
  
  // Compute matrix-vector product with fewer temporary objects
  return CVector<T>(
      (tensor[0][0] * m0 + tensor[0][1] * m1 + tensor[0][2] * m2) * scale,
      (tensor[1][0] * m0 + tensor[1][1] * m1 + tensor[1][2] * m2) * scale,
      (tensor[2][0] * m0 + tensor[2][1] * m1 + tensor[2][2] * m2) * scale);
}

template <typename T>
inline CVector<T> c_cross(const CVector<T> &a, const CVector<T> &b) {
  // Cache array accesses for better performance
  const T a0 = a[0], a1 = a[1], a2 = a[2];
  const T b0 = b[0], b1 = b[1], b2 = b[2];
  
  return CVector<T>(a1 * b2 - a2 * b1, a2 * b0 - a0 * b2, a0 * b1 - a1 * b0);
}

template <typename T> 
constexpr inline T c_dot(const CVector<T> &a, const CVector<T> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T> class EnergyDriver {
public:
  static inline T calculateZeemanEnergy(const CVector<T> &mag, const CVector<T> &Hext, T cellVolume,
                                 T Ms) {
    return -MAGNETIC_PERMEABILITY * Ms * c_dot<T>(mag, Hext) * cellVolume;
  }

  static inline T calculateAnisotropyEnergy(const CVector<T> &mag, const CVector<T> &anis, T K,
                                     T cellVolume) {
    const T dot = c_dot<T>(mag, anis);
    const T normProd = anis.length() * mag.length();
    const T sinSq = 1.0 - (dot * dot) / (normProd * normProd);
    return K * sinSq * cellVolume;
  }

  static inline T calculateIECEnergy(const CVector<T> &mag, const CVector<T> &other, T J,
                              T cellSurface) {
    return -c_dot<T>(mag, other) * J * cellSurface;
  }

  static inline T calculateDemagEnergy(const CVector<T> &mag, const CVector<T> &Hdemag, T Ms,
                                T cellVolume) {
    return -0.5 * MAGNETIC_PERMEABILITY * Ms * c_dot<T>(mag, Hdemag) *
           cellVolume;
  }
};

// Forward declarations if needed
template <typename T> class Layer;
template <typename T> class Junction;

// Typedef declarations outside the class
template <typename T>
using SolverFn = bool (Layer<T>::*)(T t, T &timeStep, const CVector<T> &bottom,
                                    const CVector<T> &top);

template <typename T>
using RunnerFn = void (Junction<T>::*)(SolverFn<T> &functor, T &t, T &timeStep,
                                       bool &step_accepted);

template <typename T> struct AdaptiveIntegrationParams {
  T abs_tol = 1e-6;             // Absolute error tolerance
  T rel_tol = 1e-3;             // Relative error tolerance
  T max_factor = 5.0;           // Maximum allowed increase in step size
  T min_factor = 0.1;           // Minimum allowed decrease in step size
  T safety_factor = 0.9;        // Safety factor for step size adjustment
  bool use_pid_control = false; // Whether to use PID control for step size
  T ki = 0.0;                   // Integral gain for PID controller
  T kp = 0.2;                   // Proportional gain for PID controller
  T kd = 0.0;                   // Derivative gain for PID controller

  // Previous error for PID controller
  T prev_error_ratio = 1.0;
  T integral_error = 0.0;
};

enum Reference { NONE = 0, FIXED, TOP, BOTTOM };

enum SolverMode { EULER_HEUN = 0, RK4 = 1, DORMAND_PRINCE = 2, HEUN = 3 };

template <typename T = double> class Layer {
private:
  ScalarDriver<T> temperatureDriver;

  // CMTJ interaction drivers
  ScalarDriver<T> IECDriverTop;
  ScalarDriver<T> IECDriverBottom;
  ScalarDriver<T> IECQuadDriverTop;
  ScalarDriver<T> IECQuadDriverBottom;
  AxialDriver<T> IDMIDriverTop;
  AxialDriver<T> IDMIDriverBottom;
  AxialDriver<T> HreservedInteractionFieldDriver;
  // CMTJ Torque & Field drivers
  ScalarDriver<T> currentDriver;
  ScalarDriver<T> anisotropyDriver, secondOrderAnisotropyDriver;
  ScalarDriver<T> fieldLikeTorqueDriver, fieldLikeTorqueDriver2;
  ScalarDriver<T> dampingLikeTorqueDriver, dampingLikeTorqueDriver2;
  AxialDriver<T> externalFieldDriver;
  AxialDriver<T> HoeDriver, HdmiDriver;

  bool nonStochasticTempSet = false;
  bool nonStochasticOneFSet = true;
  bool temperatureSet = false;
  bool pinkNoiseSet = false;
  bool alternativeSTTSet = false;
  Reference referenceType = NONE;

  // the distribution is binded for faster generation
  // is also shared between 1/f and Gaussian noise.
  std::function<T()> distribution = std::bind(
      std::normal_distribution<T>(0, 1), std::mt19937(std::random_device{}()));

  CVector<T> dWn, dWn2; // one for thermal, one for OneF
  Layer(const std::string &id, CVector<T> mag, CVector<T> anis, T Ms,
        T thickness, T cellSurface, const std::vector<CVector<T>> &demagTensor,
        T damping, T fieldLikeTorque, T dampingLikeTorque,
        T SlonczewskiSpacerLayerParameter, T beta, T spinPolarisation)
      : id(id), mag(mag), anis(anis), Ms(Ms), thickness(thickness),
        cellSurface(cellSurface), demagTensor(demagTensor), damping(damping),
        dampingSq(damping * damping),
        fieldLikeTorque(fieldLikeTorque), dampingLikeTorque(dampingLikeTorque),
        SlonczewskiSpacerLayerParameter(SlonczewskiSpacerLayerParameter),
        SlonczewskiSpacerLayerParameterSq(SlonczewskiSpacerLayerParameter * SlonczewskiSpacerLayerParameter),
        beta(beta), spinPolarisation(spinPolarisation) {
    if (mag.length() == 0) {
      throw std::runtime_error(
          "Initial magnetisation was set to a zero vector!");
    }
    if (anis.length() == 0) {
      throw std::runtime_error("Anisotropy was set to a zero vector!");
    }
    // normalise magnetisation
    mag.normalize();
    dWn = CVector<T>(this->distribution);
    dWn.normalize();
    this->cellVolume = this->cellSurface * this->thickness;
    this->ofn = std::shared_ptr<OneFNoise<T>>(new OneFNoise<T>(0, 0., 0.));
  }

public:
  struct BufferedNoiseParameters {
    /* data */
    T alphaNoise = 1.0;
    T scaleNoise = 0.0;
    T stdNoise = 0.0;
    Axis axis = Axis::all;
  };
  BufferedNoiseParameters noiseParams;
  AdaptiveIntegrationParams<T> adaptiveParams;
  std::shared_ptr<OneFNoise<T>> ofn;
  std::shared_ptr<VectorAlphaNoise<T>> bfn;
  bool includeSTT = false;
  bool includeSOT = false;

  std::string id;
  T Ms = 0.0;

  // geometric parameters
  T thickness = 0.0;
  T cellVolume = 0.0, cellSurface = 0.0;

  CVector<T> H_log, Hoe_log, Hconst, mag, anis, referenceLayer,
      secondaryReferenceLayer;
  CVector<T> Hext, Hdipole, Hdemag, Hoe, HAnis, Hthermal, Hfluctuation, Hdmi,
      Hidmi;

  CVector<T> Hfl_v, Hdl_v;

  CVector<T> HIEC, HIECtop, HIECbottom;
  T Jbottom_log = 0.0, Jtop_log = 0.0;
  T J2bottom_log = 0.0, J2top_log = 0.0;
  T K_log = 0.0;
  T K2_log = 0.0;
  T I_log = 0.0;

  // dipole and demag tensors
  std::vector<CVector<T>> demagTensor;
  std::vector<CVector<T>> dipoleBottom =
      std::vector<CVector<T>>{CVector<T>(), CVector<T>(), CVector<T>()};
  std::vector<CVector<T>> dipoleTop =
      std::vector<CVector<T>>{CVector<T>(), CVector<T>(), CVector<T>()};

  // LLG params
  T damping;
  T dampingSq; // cached damping^2 for performance

  // SOT params
  bool dynamicSOT = true;
  T fieldLikeTorque;
  T dampingLikeTorque;

  // STT params
  T SlonczewskiSpacerLayerParameter;
  T SlonczewskiSpacerLayerParameterSq; // cached for performance
  T beta;      // usually either set to 0 or to damping
  T kappa = 1; // for damping-like off -turning torque
  T spinPolarisation;

  T hopt = -1.0;

  Layer() {}
  explicit Layer(const std::string &id, const CVector<T> &mag,
                 const CVector<T> &anis, T Ms, T thickness, T cellSurface,
                 const std::vector<CVector<T>> &demagTensor, T damping)
      : Layer(id, mag, anis, Ms, thickness, cellSurface, demagTensor, damping,
              0, 0, 0, 0, 0) {}

  /**
   * The basic structure is a magnetic layer.
   * Its parameters are defined by the constructor and may be altered
   * by the drivers during the simulation time.
   * If you want STT, remember to set the reference vector for the polarisation
   * of the layer. Use `setReferenceLayer` function to do that.
   * @param id: identifiable name for a layer -- e.g. "bottom" or "free".
   * @param mag: initial magnetisation. Must be normalised (norm of 1). Used for
   * quicker convergence.
   * @param anis: anisotropy of the layer. A normalised vector
   * @param Ms: magnetisation saturation. Unit: Tesla [T].
   * @param thickness: thickness of the layer. Unit: meter [m].
   * @param cellSurface: surface of the layer, for volume calculation. Unit:
   * meter^2 [m^2].
   * @param demagTensor: demagnetisation tensor of the layer.
   * @param damping: often marked as alpha in the LLG equation. Damping of the
   * layer. Default 0.011. Dimensionless.
   * @param fieldLikeTorque: [SOT] effective spin Hall angle (spin
   * effectiveness) for Hfl.
   * @param dampingLikeTorque: [SOT] effective spin Hall angle (spin
   * effectiveness) for Hdl.
   */
  explicit Layer(const std::string &id, const CVector<T> &mag,
                 const CVector<T> &anis, T Ms, T thickness, T cellSurface,
                 const std::vector<CVector<T>> &demagTensor, T damping,
                 T fieldLikeTorque, T dampingLikeTorque)
      : Layer(id, mag, anis, Ms, thickness, cellSurface, demagTensor, damping,
              fieldLikeTorque, dampingLikeTorque, 0, 0, 0) {
    this->includeSTT = false;
    this->includeSOT = true;
    this->dynamicSOT = false;
  }

  /**
   * The basic structure is a magnetic layer.
   * Its parameters are defined by the constructor and may be altered
   * by the drivers during the simulation time.
   * If you want STT, remember to set the reference vector for the polarisation
   * of the layer. Use `setReferenceLayer` function to do that.
   * @param id: identifiable name for a layer -- e.g. "bottom" or "free".
   * @param mag: initial magnetisation. Must be normalised (norm of 1). Used for
   * quicker convergence.
   * @param anis: anisotropy of the layer. A normalised vector
   * @param Ms: magnetisation saturation. Unit: Tesla [T].
   * @param thickness: thickness of the layer. Unit: meter [m].
   * @param cellSurface: surface of the layer, for volume calculation. Unit:
   * meter^2 [m^2].
   * @param demagTensor: demagnetisation tensor of the layer.
   * @param damping: often marked as alpha in the LLG equation. Damping of the
   * layer. Default 0.011. Dimensionless.
   * @param SlomczewskiSpacerLayerParameter: [STT] Slomczewski parameter.
   * Default 1.0. Dimensionless.
   * @param beta: [STT] beta parameter for the STT. Default 0.0. Dimensionless.
   * @param spinPolarisation: [STT] polarisation ratio while passing through
   * reference layer.
   */
  explicit Layer(const std::string &id, const CVector<T> &mag,
                 const CVector<T> &anis, T Ms, T thickness, T cellSurface,
                 const std::vector<CVector<T>> &demagTensor, T damping,
                 T SlonczewskiSpacerLayerParameter, T beta, T spinPolarisation)
      : Layer(id, mag, anis, Ms, thickness, cellSurface, demagTensor, damping,
              0, 0, SlonczewskiSpacerLayerParameter, beta, spinPolarisation) {
    this->includeSTT = true;
    this->includeSOT = false;
  }

  inline static Layer<T> LayerSTT(const std::string &id, const CVector<T> &mag,
                                  const CVector<T> &anis, T Ms, T thickness,
                                  T cellSurface,
                                  const std::vector<CVector<T>> &demagTensor,
                                  T damping, T SlonczewskiSpacerLayerParameter,
                                  T beta, T spinPolarisation) {
    return Layer<T>(id, mag, anis, Ms, thickness, cellSurface, demagTensor,
                    damping, SlonczewskiSpacerLayerParameter, beta,
                    spinPolarisation);
  }

  inline static Layer<T> LayerSOT(const std::string &id, const CVector<T> &mag,
                                  const CVector<T> &anis, T Ms, T thickness,
                                  T cellSurface,
                                  const std::vector<CVector<T>> &demagTensor,
                                  T damping, T fieldLikeTorque,
                                  T dampingLikeTorque) {
    return Layer<T>(id, mag, anis, Ms, thickness, cellSurface, demagTensor,
                    damping, fieldLikeTorque, dampingLikeTorque);
  }

  void setAdaptiveParams(const AdaptiveIntegrationParams<T> &params) {
    this->adaptiveParams = params;
  }

  /**
   * @brief Get the Id object
   *
   * @return const std::string
   */
  const std::string &getId() const { return id; }
  /**
   * @brief Set the Alternative STT formulation
   *
   * @param alternativeSTT: True if you want to use the alternative STT
   * formulation.
   */
  void setAlternativeSTT(bool alternativeSTT) {
    this->alternativeSTTSet = alternativeSTT;
  }
  void setKappa(T kappa) { this->kappa = kappa; }
  void setTopDipoleTensor(const std::vector<CVector<T>> &dipoleTensor) {
    this->dipoleTop = dipoleTensor;
  }

  void setBottomDipoleTensor(const std::vector<CVector<T>> &dipoleTensor) {
    this->dipoleBottom = dipoleTensor;
  }

  const bool hasTemperature() { return this->temperatureSet; }

  void setTemperatureDriver(const ScalarDriver<T> &driver) {
    this->temperatureDriver = driver;
    this->temperatureSet = true;
  }

  void setNonStochasticLangevinDriver(const ScalarDriver<T> &driver) {
    this->temperatureDriver = driver;
    // do not set the SDE flag here
    this->temperatureSet = false;
    this->nonStochasticTempSet = true;
  }

  void setOneFNoise(unsigned int sources, T bias, T scale) {
    this->ofn =
        std::shared_ptr<OneFNoise<T>>(new OneFNoise<T>(sources, bias, scale));
    this->pinkNoiseSet = true;
    // by default turn it on, but in the stochastic sims, we will have to turn
    // it off
    this->nonStochasticOneFSet = true;
  }

  void setAlphaNoise(T alpha, T std, T scale, Axis axis = Axis::all) {
    if ((alpha < 0) || (alpha > 2))
      throw std::runtime_error("alpha must be between 0 and 2");
    this->noiseParams.alphaNoise = alpha;
    this->noiseParams.stdNoise = std;
    this->noiseParams.scaleNoise = scale;
    this->noiseParams.axis = axis;
    this->pinkNoiseSet = true;
  }

  void createBufferedAlphaNoise(unsigned int bufferSize) {
    if (this->noiseParams.alphaNoise < 0)
      throw std::runtime_error(
          "alpha must be set before creating the noise!"
          " Use setAlphaNoise function to set the alpha parameter.");

    this->bfn = std::shared_ptr<VectorAlphaNoise<T>>(new VectorAlphaNoise<T>(
        bufferSize, this->noiseParams.alphaNoise, this->noiseParams.stdNoise,
        this->noiseParams.scaleNoise, this->noiseParams.axis));
  }

  void setCurrentDriver(const ScalarDriver<T> &driver) {
    this->currentDriver = driver;
  }

  void setTorqueParameters() {
    this->includeSOT = true;
    if (this->includeSTT)
      throw std::runtime_error(
          "includeSTT was on and now setting SOT interaction!");
    if (!this->dynamicSOT)
      throw std::runtime_error(
          "used a static SOT definition, now trying to set it dynamically!");
  }

  void setFieldLikeTorqueDriver(const ScalarDriver<T> &driver) {
    setTorqueParameters();
    this->fieldLikeTorqueDriver = driver;
  }

  void setDampingLikeTorqueDriver(const ScalarDriver<T> &driver) {
    setTorqueParameters();
    this->dampingLikeTorqueDriver = driver;
  }

  void setSecondaryFieldLikeTorqueDriver(const ScalarDriver<T> &driver) {
    setTorqueParameters();
    this->fieldLikeTorqueDriver2 = driver;
  }

  void setSecondaryDampingLikeTorqueDriver(const ScalarDriver<T> &driver) {
    setTorqueParameters();
    this->dampingLikeTorqueDriver2 = driver;
  }

  void setPrimaryTorqueDrivers(const ScalarDriver<T> &fieldLikeTorqueDriver,
                               const ScalarDriver<T> &dampingLikeTorqueDriver) {
    setTorqueParameters();
    this->fieldLikeTorqueDriver = fieldLikeTorqueDriver;
    this->dampingLikeTorqueDriver = dampingLikeTorqueDriver;
  }

  void
  setSecondaryTorqueDrivers(const ScalarDriver<T> &fieldLikeTorqueDriver,
                            const ScalarDriver<T> &dampingLikeTorqueDriver) {
    setTorqueParameters();
    this->fieldLikeTorqueDriver2 = fieldLikeTorqueDriver;
    this->dampingLikeTorqueDriver2 = dampingLikeTorqueDriver;
  }

  void setAnisotropyDriver(const ScalarDriver<T> &driver) {
    this->anisotropyDriver = driver;
  }

  void setSecondOrderAnisotropyDriver(const ScalarDriver<T> &driver) {
    this->secondOrderAnisotropyDriver = driver;
  }

  void setExternalFieldDriver(const AxialDriver<T> &driver) {
    this->externalFieldDriver = driver;
  }
  void setOerstedFieldDriver(const AxialDriver<T> &driver) {
    this->HoeDriver = driver;
  }

  void setMagnetisation(const CVector<T> &newMag) {
    if (newMag.length() == 0) {
      throw std::runtime_error(
          "Initial magnetisation was set to a zero vector!");
    }
    this->mag = newMag;
    this->mag.normalize();
  }

  void setIECDriverBottom(const ScalarDriver<T> &driver) {
    this->IECDriverBottom = driver;
  }

  void setIECDriverTop(const ScalarDriver<T> &driver) {
    this->IECDriverTop = driver;
  }

  void setQuadIECDriverTop(const ScalarDriver<T> &driver) {
    this->IECQuadDriverTop = driver;
  }

  void setQuadIECDriverBottom(const ScalarDriver<T> &driver) {
    this->IECQuadDriverBottom = driver;
  }

  void setIDMIDriverTop(const AxialDriver<T> &driver) {
    this->IDMIDriverTop = driver;
  }

  void setIDMIDriverBottom(const AxialDriver<T> &driver) {
    this->IDMIDriverBottom = driver;
  }

  void setHdmiDriver(const AxialDriver<T> &driver) {
    this->HdmiDriver = driver;
  }

  void setReservedInteractionField(const AxialDriver<T> &driver) {
    this->HreservedInteractionFieldDriver = driver;
  }

  /**
   * @brief Sets reference layer with a custom vector
   * Set reference layer parameter. This is for calculating the spin current
   * polarisation if `includeSTT` is true.
   * @param reference: CVector describing the reference layer.
   */
  void setReferenceLayer(const CVector<T> &reference) {
    this->referenceLayer = reference;
    this->referenceType = FIXED;
  }

  void setSecondaryReferenceLayer(const CVector<T> &reference) {
    this->secondaryReferenceLayer = reference;
    this->referenceType = FIXED;
  }

  /**
   * @brief Set reference layer with enum
   * Can be used to refer to other layers in stack as reference
   * for this layer.
   * @param reference: an enum: FIXED, TOP, BOTTOM, or CUSTOM
   */
  void setReferenceLayer(Reference reference) {
    if ((reference == FIXED) && (!this->referenceLayer.length())) {
      throw std::runtime_error("Cannot set fixed polarisation layer to 0!"
                               " Set reference to NONE to disable reference.");
    }
    this->referenceType = reference;
  }

  /**
   * @brief Get the Reference Layer object
   */
  CVector<T> getReferenceLayer() {
    // TODO: return other mags when the reference layer is not fixed.
    return this->referenceLayer;
  }

  CVector<T> getSecondaryReferenceLayer() {
    return this->secondaryReferenceLayer;
  }

  /**
   * @brief Get the Reference Layer Type object (enum type is returned)
   */
  Reference getReferenceType() { return this->referenceType; }

  inline const CVector<T>
  calculateHeff(T time, T timeStep, const CVector<T> &stepMag,
                const CVector<T> &bottom, const CVector<T> &top,
                const CVector<T> &Hfluctuation = CVector<T>()) {
    this->Hdipole =
        calculate_tensor_interaction(bottom, this->dipoleBottom, this->Ms) +
        calculate_tensor_interaction(top, this->dipoleTop, this->Ms);
    return calculateHeffDipoleInjection(time, timeStep, stepMag, bottom, top,
                                        this->Hdipole, Hfluctuation);
  }

  const CVector<T>
  calculateHeffDipoleInjection(T time, T timeStep, const CVector<T> &stepMag,
                               const CVector<T> &bottom, const CVector<T> &top,
                               const CVector<T> &dipole,
                               const CVector<T> &Hfluctuation) {
    this->Hext = calculateExternalField(time);
    this->Hoe = calculateHOeField(time);

    this->Hdemag =
        calculate_tensor_interaction(stepMag, this->demagTensor, this->Ms);
    this->HIEC = calculateIEC(time, stepMag, bottom, top);
    this->Hidmi = calculateIDMI(time, stepMag, bottom, top);
    this->HAnis = calculateAnisotropy(stepMag, time);
    this->Hdmi = calculateHdmiField(time);
    CVector<T> HreservedInteractionField =
        this->HreservedInteractionFieldDriver.getCurrentAxialDrivers(time);
    const CVector<T> HAnis2 = calculateSecondOrderAnisotropy(stepMag, time);
    const CVector<T> Heff = this->Hext     // external
                            + this->HAnis  // anistotropy
                            + HAnis2       // second order anisotropy
                            + this->HIEC   // IEC
                            + this->Hidmi  // IDMI
                            + this->Hoe    // Oersted field
                            + this->Hdmi   // regular DMI
                            + Hfluctuation // fluctuations
                            // demag -- negative contribution
                            - this->Hdemag
                            // dipole -- negative contribution
                            - dipole
                            // reserved interaction field
                            + HreservedInteractionField;
    return Heff;
  }

  inline CVector<T> calculateHOeField(const T &time) {
    this->Hoe_log = this->HoeDriver.getCurrentAxialDrivers(time);
    return this->Hoe_log;
  }

  inline CVector<T> calculateHdmiField(const T &time) {
    return this->HdmiDriver.getCurrentAxialDrivers(time);
  }

  inline CVector<T> calculateExternalField(const T &time) {
    this->H_log = this->externalFieldDriver.getCurrentAxialDrivers(time);
    return this->H_log;
  }

  inline CVector<T> calculateAnisotropy(const CVector<T> &stepMag, T &time) {
    this->K_log = this->anisotropyDriver.getCurrentScalarValue(time);
    const T nom = (2 * this->K_log) * c_dot<T>(this->anis, stepMag) / this->Ms;
    return this->anis * nom;
  }

  inline CVector<T> calculateSecondOrderAnisotropy(const CVector<T> &stepMag,
                                            T &time) {
    this->K2_log =
        this->secondOrderAnisotropyDriver.getCurrentScalarValue(time);
    const T dot = c_dot<T>(this->anis, stepMag);
    const T nom = (4 * this->K2_log) * dot * dot * dot / this->Ms;
    return this->anis * nom;
  }

  inline CVector<T> calculateIEC_(const T J, const T J2, const CVector<T> &stepMag,
                           const CVector<T> &coupledMag) {
    // below an alternative method for computing J -- it's here for reference
    // only. const T nom = J / (this->Ms * this->thickness); return (coupledMag
    // - stepMag) * nom; // alternative form return (coupledMag + coupledMag * 2
    // * J2 * c_dot(coupledMag, stepMag)) * nom;
    const T scale = (J + 2 * J2 * c_dot(coupledMag, stepMag)) / (this->Ms * this->thickness);
    return coupledMag * scale;
  }

  inline CVector<T> calculateIEC(T time, const CVector<T> &stepMag,
                          const CVector<T> &bottom, const CVector<T> &top) {
    this->Jbottom_log = this->IECDriverBottom.getCurrentScalarValue(time);
    this->Jtop_log = this->IECDriverTop.getCurrentScalarValue(time);

    this->J2bottom_log = this->IECQuadDriverBottom.getCurrentScalarValue(time);
    this->J2top_log = this->IECQuadDriverTop.getCurrentScalarValue(time);

    return calculateIEC_(this->Jbottom_log, this->J2bottom_log, stepMag,
                         bottom) +
           calculateIEC_(this->Jtop_log, this->J2top_log, stepMag, top);
  }

  inline CVector<T> calculateIDMI_(const CVector<T> &Dvector,
                            const CVector<T> &stepMag,
                            const CVector<T> &coupledMag) {
    // D * [(dm1/dm1x x m2) + (m1 x dm2/dm2x)]
    // dm1/dm1x x m2 = (0, -mz, my)
    // dm1/dm1y x m2 = (mz, 0, -mx)
    // dm1/dm1z x m2 = (-my, mx, 0)
    // E = D z * (m1 x m2) == D m1 (m2 x z)
    // dE/dm1 = D m2 x z
    const T scale = -1.0 / (this->Ms * this->thickness);
    return c_cross<T>(Dvector, coupledMag) * scale;
    // const CVector<T> dm1crossm2(
    //     c_dot(Dvector, CVector<T>(0, -coupledMag.z, coupledMag.y)),
    //     c_dot(Dvector, CVector<T>(coupledMag.z, 0, -coupledMag.x)),
    //     c_dot(Dvector, CVector<T>(-coupledMag.y, coupledMag.x, 0)));
    // return dm1crossm2 / (this->Ms * this->thickness);
  }

  inline CVector<T> calculateIDMI(T time, const CVector<T> &stepMag,
                           const CVector<T> &bottom, const CVector<T> &top) {
    return calculateIDMI_(this->IDMIDriverBottom.getCurrentAxialDrivers(time),
                          stepMag, bottom) +
           calculateIDMI_(this->IDMIDriverTop.getCurrentAxialDrivers(time),
                          stepMag, top);
  }

  /**
   * @brief Main solver function. It is solver-independent (all solvers use this
   * function). This function is called by the solver to calculate the next step
   * of the magnetisation. It computes implicitly, all torques, given the
   * current magnetisation and effective field.
   * @param time the time at which the solver is currently at.
   * @param m the current magnetisation (from the solver, may be a semi-step)
   * @param timeStep integration time
   * @param bottom magnetisation of the layer below
   * @param top magnetisation of the layer above
   * @param heff the effective field
   * @return const CVector<T> magnetisation after the step
   */
  const CVector<T> solveLLG(T time, const CVector<T> &m, T timeStep,
                            const CVector<T> &bottom, const CVector<T> &top,
                            const CVector<T> &heff) {
    const CVector<T> prod = c_cross<T>(m, heff);
    const CVector<T> prod2 = c_cross<T>(m, prod);
    const T convTerm = 1 / (1 + this->dampingSq); // LLGS -> LL form
    const CVector<T> dmdt = prod + prod2 * this->damping;
    CVector<T> reference;

    // decide what is to be the reference for (s)LLG-STT
    // dynamically substitute other active layers
    switch (this->referenceType) {
      // TODO: add the warning if reference layer is top/bottom and empty
    case FIXED:
      reference = this->referenceLayer;
      break;
    case TOP:
      reference = top;
      break;
    case BOTTOM:
      reference = bottom;
      break;
    default:
      break;
    }

    // extra terms
    if (this->includeSTT) {
      this->I_log = this->currentDriver.getCurrentScalarValue(time);
      // use standard STT formulation
      // see that literature reports Ms/MAGNETIC_PERMEABILITY
      // but then the units don't match, we use Ms [T] which works
      const T aJ =
          HBAR * this->I_log / (ELECTRON_CHARGE * this->Ms * this->thickness);
      // field like
      T eta = 0;
      if (this->alternativeSTTSet) {
        // this is simplified
        eta = (this->spinPolarisation) /
              (1 +
               this->SlonczewskiSpacerLayerParameter * c_dot<T>(m, reference));
      } else {
        // this is more complex model (classical STT)
        const T slonSq = this->SlonczewskiSpacerLayerParameterSq;
        eta = (this->spinPolarisation * slonSq) /
              (slonSq + 1 + (slonSq - 1) * c_dot<T>(m, reference));
      }
      const T sttTerm = GYRO * aJ * eta;
      const CVector<T> fieldLike = c_cross<T>(m, reference);
      // damping like
      const CVector<T> dampingLike = c_cross<T>(m, fieldLike);
      return (dmdt * -GYRO + dampingLike * -sttTerm * this->kappa +
              fieldLike * sttTerm * this->beta) *
             convTerm;
    } else if (this->includeSOT) {
      T Hdl = 0, Hfl = 0, Hdl2 = 0, Hfl2 = 0;

      // Get SOT field values - either from drivers or using current
      if (this->dynamicSOT) {
        Hdl = this->dampingLikeTorqueDriver.getCurrentScalarValue(time);
        Hfl = this->fieldLikeTorqueDriver.getCurrentScalarValue(time);
        Hdl2 = this->dampingLikeTorqueDriver2.getCurrentScalarValue(time);
        Hfl2 = this->fieldLikeTorqueDriver2.getCurrentScalarValue(time);
      } else {
        this->I_log = this->currentDriver.getCurrentScalarValue(time);
        Hdl = this->dampingLikeTorque * this->I_log;
        Hfl = this->fieldLikeTorque * this->I_log;
      }

      // Calculate field vectors
      this->Hfl_v = reference * (Hfl - this->damping * Hdl);
      this->Hdl_v = reference * (Hdl + this->damping * Hfl);
      const CVector<T> Hfl2_v =
          secondaryReferenceLayer * (Hfl2 - this->damping * Hdl2);
      const CVector<T> Hdl2_v =
          secondaryReferenceLayer * (Hdl2 + this->damping * Hfl2);

      // Calculate torques
      const CVector<T> cm_primary = c_cross<T>(m, reference);
      const CVector<T> ccm_primary = c_cross<T>(m, cm_primary);
      const CVector<T> cm_secondary = c_cross<T>(m, secondaryReferenceLayer);
      const CVector<T> ccm_secondary = c_cross<T>(m, cm_secondary);

      // Primary and secondary torque components
      const CVector<T> flTorque_primary =
          cm_primary * (Hfl - this->damping * Hdl);
      const CVector<T> dlTorque_primary =
          ccm_primary * (Hdl + this->damping * Hfl);
      const CVector<T> flTorque_secondary =
          cm_secondary * (Hfl2 - this->damping * Hdl2);
      const CVector<T> dlTorque_secondary =
          ccm_secondary * (Hdl2 + this->damping * Hfl2);

      return (dmdt + flTorque_primary + dlTorque_primary + flTorque_secondary +
              dlTorque_secondary) *
             -GYRO * convTerm;
    }

    return dmdt * -GYRO * convTerm;
  }

  /**
   * @brief Assumes the dW has the scale of sqrt(timeStep).
   *
   * @param currentMag
   * @param dW - stochastic vector already scaled properly
   * @return CVector<T>
   */
  inline CVector<T> stochasticTorque(const CVector<T> &currentMag,
                              const CVector<T> &dW) {

    const T convTerm = -GYRO / (1. + this->dampingSq);
    const CVector<T> thcross = c_cross(currentMag, dW);
    const CVector<T> thcross2 = c_cross(currentMag, thcross);
    return (thcross + thcross2 * this->damping) * convTerm;
  }

  const CVector<T> calculateLLGWithFieldTorqueDipoleInjection(
      T time, const CVector<T> &m, const CVector<T> &bottom,
      const CVector<T> &top, const CVector<T> &dipole, T timeStep,
      const CVector<T> &Hfluctuation = CVector<T>()) {
    // classic LLG first
    const CVector<T> heff = calculateHeffDipoleInjection(
        time, timeStep, m, bottom, top, dipole, Hfluctuation);
    return solveLLG(time, m, timeStep, bottom, top, heff);
  }

  /**
   * Compute the LLG time step. The efficient field vectors is calculated
   * implicitly here. Use the effective spin hall angles formulation for SOT
   * interaction.
   * @param time: current simulation time.
   * @param m: current RK45 magnetisation.
   * @param bottom: layer below the current layer (current layer's magnetisation
   * is m). For IEC interaction.
   * @param top: layer above the current layer (current layer's magnetisation is
   * m). For IEC interaction.
   * @param timeStep: RK45 integration step.
   */
  inline const CVector<T>
  calculateLLGWithFieldTorque(T time, const CVector<T> &m,
                              const CVector<T> &bottom, const CVector<T> &top,
                              T timeStep,
                              const CVector<T> &Hfluctuation = CVector<T>()) {
    // classic LLG first
    const CVector<T> heff =
        calculateHeff(time, timeStep, m, bottom, top, Hfluctuation);
    return solveLLG(time, m, timeStep, bottom, top, heff);
  }

  bool dormand_prince_step(T time, T &timeStep, const CVector<T> &bottom,
                           const CVector<T> &top) {
    CVector<T> m_t = this->mag;

    // Constants for Dormand-Prince method
    constexpr T c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0,
                c6 = 1.0, c7 = 1.0;

    // First stage - same as RK4
    const CVector<T> k1 =
        calculateLLGWithFieldTorque(time, m_t, bottom, top, timeStep) *
        timeStep;

    // Second stage
    const CVector<T> k2 =
        calculateLLGWithFieldTorque(time + c2 * timeStep, m_t + k1 * (c2),
                                    bottom, top, timeStep) *
        timeStep;

    // Third stage
    const CVector<T> k3 =
        calculateLLGWithFieldTorque(time + c3 * timeStep,
                                    m_t + k1 * (3.0 / 40.0) + k2 * (9.0 / 40.0),
                                    bottom, top, timeStep) *
        timeStep;

    // Fourth stage
    const CVector<T> k4 =
        calculateLLGWithFieldTorque(time + c4 * timeStep,
                                    m_t + k1 * (44.0 / 45.0) +
                                        k2 * (-56.0 / 15.0) + k3 * (32.0 / 9.0),
                                    bottom, top, timeStep) *
        timeStep;

    // Fifth stage
    const CVector<T> k5 =
        calculateLLGWithFieldTorque(
            time + c5 * timeStep,
            m_t + k1 * (19372.0 / 6561.0) + k2 * (-25360.0 / 2187.0) +
                k3 * (64448.0 / 6561.0) + k4 * (-212.0 / 729.0),
            bottom, top, timeStep) *
        timeStep;

    // Sixth stage
    const CVector<T> k6 =
        calculateLLGWithFieldTorque(
            time + c6 * timeStep,
            m_t + k1 * (9017.0 / 3168.0) + k2 * (-355.0 / 33.0) +
                k3 * (46732.0 / 5247.0) + k4 * (49.0 / 176.0) +
                k5 * (-5103.0 / 18656.0),
            bottom, top, timeStep) *
        timeStep;

    // Seventh stage (k7) - will be used as k1 for next step (FSAL property)
    const CVector<T> k7 = calculateLLGWithFieldTorque(
                              time + c7 * timeStep,
                              m_t + k1 * (35.0 / 384.0) +
                                  k3 * (500.0 / 1113.0) + k4 * (125.0 / 192.0) +
                                  k5 * (-2187.0 / 6784.0) + k6 * (11.0 / 84.0),
                              bottom, top, timeStep) *
                          timeStep;

    // 5th order solution
    CVector<T> m_t5 = m_t + k1 * (35.0 / 384.0) + k3 * (500.0 / 1113.0) +
                      k4 * (125.0 / 192.0) + k5 * (-2187.0 / 6784.0) +
                      k6 * (11.0 / 84.0);

    // 4th order solution
    CVector<T> m_t4 = m_t + k1 * (5179.0 / 57600.0) + k3 * (7571.0 / 16695.0) +
                      k4 * (393.0 / 640.0) + k5 * (-92097.0 / 339200.0) +
                      k6 * (187.0 / 2100.0) + k7 * (1.0 / 40.0);

    // Error estimation
    CVector<T> error = m_t5 - m_t4;
    T error_norm = error.length();

    // Calculate tolerance based on absolute and relative tolerances
    T tol = adaptiveParams.abs_tol + adaptiveParams.rel_tol * m_t.length();
    T error_ratio = error_norm / tol;

    // Prevent division by zero or other numerical issues
    if (tol < 1e-15 || isnan(error_ratio) || isinf(error_ratio)) {
      error_ratio = 10.0; // Force a smaller timestep but avoid NaN
    }

    // Step size control
    T factor;
    if (adaptiveParams.use_pid_control) {
      // PID controller for step size
      // TODO: Test it a bit more
      T e_n = log(error_ratio);
      adaptiveParams.integral_error += e_n;
      factor = exp(adaptiveParams.kp * e_n +
                   adaptiveParams.ki * adaptiveParams.integral_error +
                   adaptiveParams.kd * (e_n - adaptiveParams.prev_error_ratio));
      adaptiveParams.prev_error_ratio = e_n;
    } else {
      // Standard controller (error-based)
      factor = adaptiveParams.safety_factor * pow(1.0 / error_ratio, 0.2);
    }

    // Check for numerical errors in factor calculation
    if (isnan(factor) || isinf(factor)) {
      factor = adaptiveParams.min_factor; // Default to reducing timestep
    }

    // Limit the factor to avoid too large/small steps
    factor = std::min(factor, adaptiveParams.max_factor);
    factor = std::max(factor, adaptiveParams.min_factor);

    // Calculate new step size
    T new_timestep = timeStep * factor;

    // Ensure new timestep is valid
    if (isnan(new_timestep) || isinf(new_timestep) || new_timestep <= 0) {
      new_timestep = timeStep * adaptiveParams.min_factor;
    }

    // Accept or reject the step
    if (error_ratio <= 1.0) {
      // Accept step - update magnetization (using 5th order solution)
      m_t5.normalize();
      this->mag = m_t5;

      if (isnan(this->mag.x)) {
        throw std::runtime_error("NAN magnetisation");
      }

      // Update the timeStep for next iteration
      timeStep = new_timestep;
      return true; // Step was accepted
    } else {
      // Reject step - don't update magnetization, just reduce the step size
      timeStep = new_timestep;
      return false; // Step was rejected
    }
  }

  /**
   * @brief RK4 step of the LLG equation.
   * Compute the LLG time step. The efficient field vectors is calculated
   * implicitly here. Use the effective spin hall angles formulation for SOT
   * interaction.
   * @param time: current simulation time.
   * @param m: current RK45 magnetisation.
   * @param bottom: layer below the current layer (current layer's magnetisation
   * is m). For IEC interaction.
   * @param top: layer above the current layer (current layer's magnetisation is
   * m). For IEC interaction.
   * @param timeStep: RK45 integration step.
   */
  bool rk4_step(T time, T &timeStep, const CVector<T> &bottom,
                const CVector<T> &top) {
    CVector<T> m_t = this->mag;
    const CVector<T> k1 =
        calculateLLGWithFieldTorque(time, m_t, bottom, top, timeStep) *
        timeStep;
    const CVector<T> k2 =
        calculateLLGWithFieldTorque(time + 0.5 * timeStep, m_t + k1 * 0.5,
                                    bottom, top, timeStep) *
        timeStep;
    const CVector<T> k3 =
        calculateLLGWithFieldTorque(time + 0.5 * timeStep, m_t + k2 * 0.5,
                                    bottom, top, timeStep) *
        timeStep;
    const CVector<T> k4 = calculateLLGWithFieldTorque(time + timeStep, m_t + k3,
                                                      bottom, top, timeStep) *
                          timeStep;
    m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
    m_t.normalize();
    this->mag = m_t;
    if (isnan(this->mag.x)) {
      throw std::runtime_error("NAN magnetisation");
    }
    return true;
  }

  /**
   * @brief RK4 step of the LLG equation if dipole injection is present.
   * Compute the LLG time step. The efficient field vectors is calculated
   * implicitly here. Use the effective spin hall angles formulation for SOT
   * interaction.
   * @param time: current simulation time.
   * @param m: current RK45 magnetisation.
   * @param bottom: layer below the current layer (current layer's magnetisation
   * is m). For IEC interaction.
   * @param top: layer above the current layer (current layer's magnetisation is
   * m). For IEC interaction.
   * @param timeStep: RK45 integration step.
   */
  void rk4_stepDipoleInjection(T time, T &timeStep, const CVector<T> &bottom,
                               const CVector<T> &top,
                               const CVector<T> &dipole) {
    CVector<T> m_t = this->mag;
    const CVector<T> k1 = calculateLLGWithFieldTorqueDipoleInjection(
                              time, m_t, bottom, top, dipole, timeStep) *
                          timeStep;
    const CVector<T> k2 = calculateLLGWithFieldTorqueDipoleInjection(
                              time + 0.5 * timeStep, m_t + k1 * 0.5, bottom,
                              top, dipole, timeStep) *
                          timeStep;
    const CVector<T> k3 = calculateLLGWithFieldTorqueDipoleInjection(
                              time + 0.5 * timeStep, m_t + k2 * 0.5, bottom,
                              top, dipole, timeStep) *
                          timeStep;
    const CVector<T> k4 =
        calculateLLGWithFieldTorqueDipoleInjection(
            time + timeStep, m_t + k3, bottom, top, dipole, timeStep) *
        timeStep;
    m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
    m_t.normalize();
    this->mag = m_t;
  }

  inline CVector<T> stochastic_llg(const CVector<T> &cm, T time, T timeStep,
                            const CVector<T> &bottom, const CVector<T> &top,
                            const CVector<T> &dW, const CVector<T> &dW2,
                            const T &HoneF) {
    // compute the Langevin fluctuations -- this is the sigma
    const T convTerm = -GYRO / (1 + this->dampingSq);
    const T Hthermal_temp =
        this->getLangevinStochasticStandardDeviation(time, timeStep);
    const CVector<T> thcross = c_cross(cm, dW);
    const CVector<T> thcross2 = c_cross(cm, thcross);
    const T scalingTh = Hthermal_temp * convTerm;

    // compute 1/f noise term
    const CVector<T> onefcross = c_cross(cm, dW2);
    const CVector<T> onefcross2 = c_cross(cm, onefcross);
    const T scalingOneF = HoneF * convTerm;

    return (thcross + thcross2 * this->damping) * scalingTh +
           (onefcross + onefcross2 * this->damping) * scalingOneF;
  }

  inline const T getStochasticOneFNoise(T time) {
    if (!this->pinkNoiseSet)
      return 0;
    else if (this->noiseParams.scaleNoise != 0) {
      // use buffered noise if available
      return this->bfn->tick();
    }
    return this->ofn->tick();
  }

  inline T getLangevinStochasticStandardDeviation(T time, T timeStep) {
    if (this->cellVolume == 0.0)
      throw std::runtime_error(
          "Cell surface cannot be 0 during temp. calculations!");
    const T currentTemp = this->temperatureDriver.getCurrentScalarValue(time);
    const T mainFactor = (2 * this->damping * BOLTZMANN_CONST * currentTemp) /
                         (this->Ms * this->cellVolume * GYRO);
    return sqrt(mainFactor);
  }

  inline CVector<T> getStochasticLangevinVector(const T &time, const T &timeStep) {
    if (!this->temperatureSet)
      return CVector<T>();
    const T Hthermal_temp =
        this->getLangevinStochasticStandardDeviation(time, timeStep);
    return CVector<T>(this->distribution) * Hthermal_temp;
  }

  inline CVector<T> getOneFVector() {
    if (this->noiseParams.scaleNoise != 0) {
      // use buffered noise if available
      return this->bfn->tickVector();
    }
    return CVector<T>();
  }
};

template <typename T> class Junction {
  friend class Layer<T>;
  const std::vector<std::string> vectorNames = {"x", "y", "z"};

public:
  enum MRmode { NONE = 0, CLASSIC = 1, STRIP = 2 };

  MRmode MR_mode;
  std::vector<Layer<T>> layers;
  T Rp = 0.0, Rap = 0.0;

  std::vector<T> Rx0, Ry0, AMR_X, AMR_Y, SMR_X, SMR_Y, AHE;
  std::unordered_map<std::string, std::vector<T>> log;

  unsigned int logLength = 0;
  unsigned int layerNo;
  std::string Rtag = "R";

  Junction() {}

  /**
   * @brief Create a plain junction.
   * No magnetoresistance is calculated.
   * @param layersToSet: layers that compose the junction
   */
  explicit Junction(std::vector<Layer<T>> layersToSet) 
      : layers(std::move(layersToSet)), MR_mode(NONE), layerNo(layers.size()) {
    if (this->layerNo == 0) {
      throw std::invalid_argument("Passed a zero length Layer vector!");
    }
    // verify that all layers have unique ids
    std::unordered_set<std::string> _ids;
    for (const auto &layer : this->layers) {
      if (_ids.find(layer.id) != _ids.end()) {
        throw std::invalid_argument("Layers must have unique ids!");
      }
      _ids.insert(layer.id);
    }
  }
  explicit Junction(std::vector<Layer<T>> layersToSet, T Rp, T Rap)
      : Junction(std::move(layersToSet)) {
    if (this->layerNo == 1) {
      // we need to check if this layer has a reference layer.
      if (!this->layers[0].referenceLayer.length()) {
        throw std::invalid_argument("MTJ with a single layer must have"
                                    " a pinning (referenceLayer) set!");
      }
    }
    if (this->layerNo > 2) {
      throw std::invalid_argument(
          "This constructor supports only bilayers!"
          " Choose the other one with the strip resistance!");
    }
    this->Rp = Rp;
    this->Rap = Rap;
    this->MR_mode = CLASSIC;
    // A string representing the tag for the junction's resistance value.
    if (this->layerNo == 2)
      this->Rtag = "R_" + this->layers[0].id + "_" + this->layers[1].id;
  }

  /**
   * Creates a junction with a STRIP magnetoresistance.
   * Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
   * length of the layers passed (they directly correspond to each layer).
   * Calculates the magnetoresistance as per: __see reference__:
   * Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
   * @param Rx0
   * @param Ry0
   * @param AMR_X
   * @param AMR_Y
   * @param SMR_X
   * @param SMR_Y
   * @param AHE
   */
  explicit Junction(std::vector<Layer<T>> layersToSet,
                    std::vector<T> Rx0, std::vector<T> Ry0,
                    std::vector<T> AMR_X, std::vector<T> AMR_Y,
                    std::vector<T> SMR_X, std::vector<T> SMR_Y,
                    std::vector<T> AHE)
      : layers(std::move(layersToSet)), 
        Rx0(std::move(Rx0)), Ry0(std::move(Ry0)), AMR_X(std::move(AMR_X)),
        AMR_Y(std::move(AMR_Y)), SMR_X(std::move(SMR_X)),
        SMR_Y(std::move(SMR_Y)), AHE(std::move(AHE))

  {
    this->layerNo = this->layers.size();
    if (this->layerNo == 0) {
      throw std::invalid_argument("Passed a zero length Layer vector!");
    }
    if ((this->layerNo != (unsigned int)this->Rx0.size()) ||
        (this->layerNo != (unsigned int)this->Ry0.size()) ||
        (this->layerNo != (unsigned int)this->AMR_X.size()) ||
        (this->layerNo != (unsigned int)this->AMR_Y.size()) ||
        (this->layerNo != (unsigned int)this->AHE.size()) ||
        (this->layerNo != (unsigned int)this->SMR_X.size()) ||
        (this->layerNo != (unsigned int)this->SMR_Y.size())) {
      throw std::invalid_argument(
          "Layers and Rx0, Ry, AMR, AMR and SMR must be of the same size!");
    }
    // this->fileSave = std::move(filename);
    this->MR_mode = STRIP;
  }

  /**
   * @brief Get Ids of the layers in the junction.
   * @return vector of layer ids.
   */
  std::vector<std::string> getLayerIds() const {
    std::vector<std::string> ids;
    ids.reserve(this->layers.size());
    std::transform(this->layers.begin(), this->layers.end(),
                   std::back_inserter(ids),
                   [](const Layer<T> &layer) { return layer.id; });
    return ids;
  }

  /**
   * Clears the simulation log.
   **/
  void clearLog() {
    this->log.clear();
    this->logLength = 0;
  }

  std::unordered_map<std::string, std::vector<T>> &getLog() {
    return this->log;
  }

  typedef void (Layer<T>::*scalarDriverSetter)(const ScalarDriver<T> &driver);
  typedef void (Layer<T>::*axialDriverSetter)(const AxialDriver<T> &driver);
  void scalarlayerSetter(const std::string &layerID, scalarDriverSetter functor,
                         ScalarDriver<T> driver) {
    bool found = false;
    for (auto &l : this->layers) {
      if (l.id == layerID || layerID == "all") {
        (l.*functor)(driver);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  void axiallayerSetter(const std::string &layerID, axialDriverSetter functor,
                        AxialDriver<T> driver) {
    bool found = false;
    for (auto &l : this->layers) {
      if (l.id == layerID || layerID == "all") {
        (l.*functor)(driver);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  /**
   * Set coupling between two layers.
   * The names of the params are only for convention. The coupling will be set
   * between bottomLayer or topLayer, order is irrelevant.
   * @param bottomLayer: the first layer id
   * @param topLayer: the second layer id
   */
  void setCouplingDriver(
      const std::string &bottomLayer, const std::string &topLayer,
      const ScalarDriver<T> &driver,
      void (Layer<T>::*setDriverFuncTop)(const ScalarDriver<T> &),
      void (Layer<T>::*setDriverFuncBottom)(const ScalarDriver<T> &)) {
    bool found = false;
    for (unsigned int i = 0; i < this->layerNo - 1; i++) {
      // check if the layer above is actually top layer the user specified
      if ((this->layers[i].id == bottomLayer) &&
          (this->layers[i + 1].id == topLayer)) {
        (this->layers[i].*setDriverFuncTop)(driver);
        (this->layers[i + 1].*setDriverFuncBottom)(driver);
        found = true;
        break;
      } else if ((this->layers[i].id == topLayer) &&
                 (this->layers[i + 1].id == bottomLayer)) {
        (this->layers[i].*setDriverFuncTop)(driver);
        (this->layers[i + 1].*setDriverFuncBottom)(driver);
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to match the layer order or find layer ids: " + bottomLayer +
          " and " + topLayer + "!");
    }
  }

  /**
   * Set coupling between two layers with an AxialDriver
   * The names of the params are only for convention. The coupling will be set
   * between bottomLayer or topLayer, order is irrelevant.
   * @param bottomLayer: the first layer id
   * @param topLayer: the second layer id
   */
  void setCouplingDriverAxial(
      const std::string &bottomLayer, const std::string &topLayer,
      const AxialDriver<T> &driver,
      void (Layer<T>::*setDriverFuncTop)(const AxialDriver<T> &),
      void (Layer<T>::*setDriverFuncBottom)(const AxialDriver<T> &)) {
    bool found = false;
    for (unsigned int i = 0; i < this->layerNo - 1; i++) {
      // check if the layer above is actually top layer the user specified
      if ((this->layers[i].id == bottomLayer) &&
          (this->layers[i + 1].id == topLayer)) {
        (this->layers[i].*setDriverFuncTop)(driver);
        (this->layers[i + 1].*setDriverFuncBottom)(driver);
        found = true;
        break;
      } else if ((this->layers[i].id == topLayer) &&
                 (this->layers[i + 1].id == bottomLayer)) {
        (this->layers[i].*setDriverFuncTop)(driver);
        (this->layers[i + 1].*setDriverFuncBottom)(driver);
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to match the layer order or find layer ids: " + bottomLayer +
          " and " + topLayer + "!");
    }
  }

  void setLayerTemperatureDriver(const std::string &layerID,
                                 const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setTemperatureDriver, driver);
  }
  void setLayerNonStochasticLangevinDriver(const std::string &layerID,
                                           const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setNonStochasticLangevinDriver,
                      driver);
  }
  void setLayerAnisotropyDriver(const std::string &layerID,
                                const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setAnisotropyDriver, driver);
  }
  void setLayerSecondOrderAnisotropyDriver(const std::string &layerID,
                                           const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setSecondOrderAnisotropyDriver,
                      driver);
  }
  void setLayerExternalFieldDriver(const std::string &layerID,
                                   const AxialDriver<T> &driver) {
    axiallayerSetter(layerID, &Layer<T>::setExternalFieldDriver, driver);
  }
  void setLayerOerstedFieldDriver(const std::string &layerID,
                                  const AxialDriver<T> &driver) {
    axiallayerSetter(layerID, &Layer<T>::setOerstedFieldDriver, driver);
  }
  void setLayerCurrentDriver(const std::string &layerID,
                             const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setCurrentDriver, driver);
  }
  void setLayerDampingLikeTorqueDriver(const std::string &layerID,
                                       const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setDampingLikeTorqueDriver, driver);
  }
  void setLayerFieldLikeTorqueDriver(const std::string &layerID,
                                     const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setFieldLikeTorqueDriver, driver);
  }

  void setLayerSecondaryFieldLikeTorqueDriver(const std::string &layerID,
                                              const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setSecondaryFieldLikeTorqueDriver,
                      driver);
  }

  void setLayerSecondaryDampingLikeTorqueDriver(const std::string &layerID,
                                                const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &Layer<T>::setSecondaryDampingLikeTorqueDriver,
                      driver);
  }

  void
  setLayerPrimaryTorqueDrivers(const std::string &layerID,
                               const ScalarDriver<T> &fieldLikeTorqueDriver,
                               const ScalarDriver<T> &dampingLikeTorqueDriver) {
    scalarlayerSetter(layerID, &Layer<T>::setFieldLikeTorqueDriver,
                      fieldLikeTorqueDriver);
    scalarlayerSetter(layerID, &Layer<T>::setDampingLikeTorqueDriver,
                      dampingLikeTorqueDriver);
  }

  void setLayerSecondaryTorqueDrivers(
      const std::string &layerID, const ScalarDriver<T> &fieldLikeTorqueDriver,
      const ScalarDriver<T> &dampingLikeTorqueDriver) {
    scalarlayerSetter(layerID, &Layer<T>::setSecondaryFieldLikeTorqueDriver,
                      fieldLikeTorqueDriver);
    scalarlayerSetter(layerID, &Layer<T>::setSecondaryDampingLikeTorqueDriver,
                      dampingLikeTorqueDriver);
  }

  void setLayerReservedInteractionField(const std::string &layerID,
                                        const AxialDriver<T> &driver) {
    axiallayerSetter(layerID, &Layer<T>::setReservedInteractionField, driver);
  }

  void setLayerHdmiDriver(const std::string &layerID,
                          const AxialDriver<T> &driver) {
    axiallayerSetter(layerID, &Layer<T>::setHdmiDriver, driver);
  }

  void setLayerAlternativeSTT(const std::string &layerID,
                              const bool alternative) {
    if (layerID == "all") {
      for (auto &l : this->layers) {
        l.setAlternativeSTT(alternative);
      }
    } else
      getLayer(layerID).setAlternativeSTT(alternative);
  }

  void setLayerOneFNoise(const std::string &layerID, unsigned int sources,
                         T bias, T scale) {

    if (layerID == "all") {
      for (auto &l : this->layers) {
        l.setOneFNoise(sources, bias, scale);
      }
    } else
      getLayer(layerID).setOneFNoise(sources, bias, scale);
  }

  /**
   * Set IDMI interaction between two layers.
   * The names of the params are only for convention. The IDMI will be set
   * between bottomLayer or topLayer, order is irrelevant.
   * See Arregi et al, Nat. Comm. 2022: Large interlayer Dzyaloshinskii-Moriya
   * interactions across Ag-layers
   * @param bottomLayer: the first layer id
   * @param topLayer: the second layer id
   */
  void setIDMIDriver(const std::string &bottomLayer,
                     const std::string &topLayer,
                     const AxialDriver<T> &driver) {
    setCouplingDriverAxial(bottomLayer, topLayer, driver,
                           &Layer<T>::setIDMIDriverTop,
                           &Layer<T>::setIDMIDriverBottom);
  }

  /**
   * Set biquadratic IEC interaction between two layers.
   * The names of the params are only for convention. The IEC will be set
   * between bottomLayer or topLayer, order is irrelevant.
   * @param bottomLayer: the first layer id
   * @param topLayer: the second layer id
   */
  void setQuadIECDriver(const std::string &bottomLayer,
                        const std::string &topLayer,
                        const ScalarDriver<T> &driver) {
    setCouplingDriver(bottomLayer, topLayer, driver,
                      &Layer<T>::setQuadIECDriverTop,
                      &Layer<T>::setQuadIECDriverBottom);
  }

  /**
   * Set blilinear IEC interaction between two layers.
   * The names of the params are only for convention. The IEC will be set
   * between bottomLayer or topLayer, order is irrelevant.
   * @param bottomLayer: the first layer id
   * @param topLayer: the second layer id
   */
  void setIECDriver(const std::string &bottomLayer, const std::string &topLayer,
                    const ScalarDriver<T> &driver) {
    setCouplingDriver(bottomLayer, topLayer, driver, &Layer<T>::setIECDriverTop,
                      &Layer<T>::setIECDriverBottom);
  }

  void setLayerMagnetisation(const std::string &layerID, CVector<T> &mag) {
    bool found = false;
    for (auto &l : this->layers) {
      if (l.id == layerID || layerID == "all") {
        l.setMagnetisation(mag);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  CVector<T> getLayerMagnetisation(const std::string &layerID) {
    return getLayer(layerID).mag;
  }

  Reference getLayerReferenceType(const std::string &layerID) {
    return getLayer(layerID).referenceType;
  }

  void setLayerReferenceLayer(const std::string &layerID,
                              const CVector<T> &referenceLayer) {
    if (layerID == "all") {
      for (auto &l : this->layers) {
        l.setReferenceLayer(referenceLayer);
      }
    } else
      getLayer(layerID).setReferenceLayer(referenceLayer);
  }

  void setLayerReferenceType(const std::string &layerID,
                             Reference referenceType) {
    if (layerID == "all") {
      for (auto &l : this->layers) {
        l.setReferenceLayer(referenceType);
      }
    } else
      getLayer(layerID).setReferenceLayer(referenceType);
  }

  Layer<T> &getLayer(const std::string &layerID) {
    const auto res = std::find_if(
        this->layers.begin(), this->layers.end(),
        [&layerID](const auto &l) -> bool { return (l.id == layerID); });
    if (res != this->layers.end()) {
      return *res;
    }
    throw std::runtime_error("Failed to find a layer with a given id " +
                             layerID + "!");
  }

  /**
   * @brief Log computed layer parameters.
   * This function logs all the necessayr parameters of the layers.
   * @param t: current time
   * @param timeStep: timeStep of the simulation (unsued for now)
   * @param calculateEnergies: if true, also include fields for energy
   * computation.
   */
  void logLayerParams(T &t, T timeStep, bool calculateEnergies = false) {
    for (const auto &layer : this->layers) {
      const std::string lId = layer.id;

      if (calculateEnergies) {
        // TODO: avoid recomputation at a cost of a slight error
        // recompute the current Heff to avoid shadow persistence of the layer
        // parameters const CVector<T> heff = calculateHeff(t, timeStep,
        // layer.m, layer.bottom, layer.top);
        this->log[lId + "_K"].emplace_back(layer.K_log);
        this->log[lId + "_K2"].emplace_back(layer.K2_log);
        this->log[lId + "_Jbottom"].emplace_back(layer.Jbottom_log);
        this->log[lId + "_Jtop"].emplace_back(layer.Jtop_log);
        this->log[lId + "_I"].emplace_back(layer.I_log);
        for (int i = 0; i < 3; i++) {
          this->log[lId + "_Hext" + vectorNames[i]].emplace_back(layer.Hext[i]);
          this->log[lId + "_Hiec" + vectorNames[i]].emplace_back(layer.HIEC[i]);
          this->log[lId + "_Hanis" + vectorNames[i]].emplace_back(
              layer.HAnis[i]);

          this->log[lId + "_Hdemag" + vectorNames[i]].emplace_back(
              layer.Hdemag[i]);
          this->log[lId + "_Hth" + vectorNames[i]].emplace_back(
              layer.Hfluctuation[i]);
          if (layer.includeSOT) {
            this->log[lId + "_Hfl" + vectorNames[i]].emplace_back(
                layer.Hfl_v[i]);
            this->log[lId + "_Hdl" + vectorNames[i]].emplace_back(
                layer.Hdl_v[i]);
          }
        }
        if (layer.includeSTT | layer.includeSOT)
          this->log[lId + "_I"].emplace_back(layer.I_log);
      }
      // always save magnetisation
      for (int i = 0; i < 3; i++) {
        this->log[lId + "_m" + vectorNames[i]].emplace_back(layer.mag[i]);
      }
    }
    if (this->MR_mode == CLASSIC && this->layerNo == 1) {
      this->log["R"].emplace_back(calculateMagnetoresistance(
          c_dot<T>(layers[0].mag, layers[0].referenceLayer)));
    } else if (MR_mode == CLASSIC && this->layerNo > 1) {
      const auto magnetoresistance = calculateMagnetoresistance(
          c_dot<T>(this->layers[0].mag, this->layers[1].mag));
      this->log[this->Rtag].emplace_back(magnetoresistance);
    } else if (MR_mode == STRIP) {
      const auto magnetoresistance =
          stripMagnetoResistance(this->Rx0, this->Ry0, this->AMR_X, this->SMR_X,
                                 this->AMR_Y, this->SMR_Y, this->AHE);
      this->log["Rx"].emplace_back(magnetoresistance[0]);
      this->log["Ry"].emplace_back(magnetoresistance[1]);
      this->log["Rz"].emplace_back(magnetoresistance[2]);
    }
    this->log["time"].emplace_back(t);
    this->logLength++;
  }

  void saveLogs(const std::string &filename) {
    if (filename == "") {
      // if there's an empty fn, don't save
      throw std::runtime_error("The filename may not be empty!");
    }
    std::ofstream logFile;
    logFile.open(filename);
    for (const auto &keyPair : this->log) {
      logFile << keyPair.first << ";";
    }
    logFile << "\n";
    for (unsigned int i = 0; i < logLength; i++) {
      for (const auto &keyPair : this->log) {
        logFile << keyPair.second[i] << ";";
      }
      logFile << "\n";
    }
    logFile.close();
  }

  /**
   * @brief Run Euler-Heun or RK4 method for a single layer.
   *
   * The Euler-Heun method should only be used
   * for stochastic simulations where the temperature
   * driver is set.
   * @param functor: solver function.
   * @param t: current time
   * @param timeStep: integration step
   */
  void runSingleLayerSolver(SolverFn<T> &functor, T &t, T &timeStep,
                            bool &step_accepted) {
    const CVector<T> dummy(0, 0, 0);
    step_accepted = (layers[0].*functor)(t, timeStep, dummy, dummy);
  }

  /**
   * @brief Select a solver based on the setup.
   *
   * Multilayer layer solver iteration.
   * @param functor: solver function.
   * @param t: current time
   * @param timeStep: integration step
   * */
  void runMultiLayerSolver(SolverFn<T> &functor, T &t, T &timeStep,
                           bool &step_accepted) {
    // Run solver for each layer and check if all steps were accepted
    step_accepted = true;
    std::vector<CVector<T>> magCopies;
    magCopies.reserve(this->layerNo + 2);
    magCopies.emplace_back(); // First layer gets 0 vector coupled
    // the first and the last layer get 0 vector coupled
    for (unsigned int i = 0; i < this->layerNo; i++)
      magCopies.emplace_back(this->layers[i].mag);
    magCopies.emplace_back(); // Last layer gets 0 vector coupled

    for (unsigned int i = 0; i < layerNo; i++) {
      // If any layer rejects the step, the whole step is rejected
      if (!(layers[i].*functor)(t, timeStep, magCopies[i], magCopies[i + 2])) {
        step_accepted = false;
      }
    }

    // If any layer rejected, we need the smallest timestep among all layers
    if (!step_accepted) {
      T min_timestep = timeStep;
      for (unsigned int i = 0; i < layerNo; i++) {
        // Assuming each layer has its own timeStep that was updated
        min_timestep = std::min(min_timestep, layers[i].hopt);
      }
      timeStep = min_timestep;
    }
  }

  void eulerHeunSolverStep(SolverFn<T> &functor, T &t, T &timeStep,
                           bool &step_accepted) {
    /*
        Euler Heun method (stochastic heun)

        y_np = y + g(y,t,dW)*dt
        g_sp = g(y_np,t+1,dW)
        y(t+1) = y + dt*f(y,t) + .5*(g(y,t,dW)+g_sp)*sqrt(dt)

        with f being the non-stochastic part and g the stochastic part
    */
    // draw the noise for each layer, dW
    std::vector<CVector<T>> mPrime;
    mPrime.reserve(this->layerNo);
    for (unsigned int i = 0; i < this->layerNo; i++) {
      // todo: after you're done, double check the thermal magnitude and dt
      // scaling there
      const CVector<T> dW =
          this->layers[i].getStochasticLangevinVector(t, timeStep) +
          this->layers[i].getOneFVector();
      const CVector<T> bottom =
          (i == 0) ? CVector<T>() : this->layers[i - 1].mag;
      const CVector<T> top =
          (i == this->layerNo - 1) ? CVector<T>() : this->layers[i + 1].mag;

      const CVector<T> fnApprox = this->layers[i].calculateLLGWithFieldTorque(
          t, this->layers[i].mag, bottom, top, timeStep);
      const CVector<T> gnApprox =
          this->layers[i].stochasticTorque(this->layers[i].mag, dW);

      // theoretically we have 2 options
      // 1. calculate only the stochastic part with the second approximation
      // 2. calculate the second approximation of m with the stochastic and
      // non-stochastic
      //    part and then use if for torque est.
      const CVector<T> mNext = this->layers[i].mag + gnApprox * sqrt(timeStep);
      const CVector<T> gnPrimeApprox =
          this->layers[i].stochasticTorque(mNext, dW);
      mPrime.emplace_back(this->layers[i].mag + fnApprox * timeStep +
                  0.5 * (gnApprox + gnPrimeApprox) * sqrt(timeStep));
    }

    for (unsigned int i = 0; i < this->layerNo; i++) {
      this->layers[i].mag = mPrime[i];
      this->layers[i].mag.normalize();
    }
  }

  void heunSolverStep(SolverFn<T> &functor, T &t, T &timeStep,
                      bool &step_accepted) {
    /*
        Heun method
        y'(t+1) = y(t) + dy(y, t)
        y(t+1) = y(t) + 0.5 * (dy(y, t) + dy(y'(t+1), t+1))
    */
    /*
        Stochastic Heun method
        y_np = y + g(y,t,dW)*dt
        g_sp = g(y_np,t+1,dW)
        y' = y_n + f_n * dt + g_n * dt
        f' = f(y, )
        y(t+1) = y + dt*f(y,t) + .5*(g(y,t,dW)+g_sp)*sqrt(dt)
    */
    std::vector<CVector<T>> fn, gn, dW, mNext;
    fn.reserve(this->layerNo);
    gn.reserve(this->layerNo);
    dW.reserve(this->layerNo);
    mNext.reserve(this->layerNo);
    // first approximation

    // make sure that
    // 1. Thermal field is added if needed
    // 2. One/f noise is added if needed
    // 3. The timestep is correctly multiplied

    for (unsigned int i = 0; i < this->layerNo; i++) {
      const CVector<T> bottom =
          (i == 0) ? CVector<T>() : this->layers[i - 1].mag;
      const CVector<T> top =
          (i == this->layerNo - 1) ? CVector<T>() : this->layers[i + 1].mag;

      fn.emplace_back(this->layers[i].calculateLLGWithFieldTorque(
          t, this->layers[i].mag, bottom, top, timeStep));

      // draw the noise for each layer, dW
      dW.emplace_back(this->layers[i].getStochasticLangevinVector(t, timeStep) +
              this->layers[i].getOneFVector());
      gn.emplace_back(this->layers[i].stochasticTorque(this->layers[i].mag, dW[i]));

      mNext.emplace_back(
          this->layers[i].mag + fn[i] * timeStep + gn[i] * sqrt(timeStep));
    }
    // second approximation
    for (unsigned int i = 0; i < this->layerNo; i++) {
      const CVector<T> bottom = (i == 0) ? CVector<T>() : mNext[i - 1];
      const CVector<T> top =
          (i == this->layerNo - 1) ? CVector<T>() : mNext[i + 1];
      // first approximation is already multiplied by timeStep
      this->layers[i].mag =
          this->layers[i].mag +
          0.5 * timeStep *
              (fn[i] + this->layers[i].calculateLLGWithFieldTorque(
                           t + timeStep, mNext[i], bottom, top, timeStep)) +
          0.5 * (gn[i] + this->layers[i].stochasticTorque(mNext[i], dW[i])) *
              sqrt(timeStep);
      // normalise
      this->layers[i].mag.normalize();
    }
  }

  /**
   * @brief Calculate strip magnetoresistance for multilayer.
   *
   * Used when MR_MODE == STRIP
   * Magnetoresistance as per:
   * Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
   * Each of the Rx0, Ry, AMR, AMR and SMR is list matching the
   * length of the layers passed (they directly correspond to each layer).
   * Calculates the magnetoresistance as per: __see reference__:
   * Spin Hall magnetoresistance in metallic bilayers by Kim, J. et al.
   * @param Rx0
   * @param Ry0
   * @param AMR_X
   * @param AMR_Y
   * @param SMR_X
   * @param SMR_Y
   * @param AHE
   */
  inline std::vector<T> stripMagnetoResistance(const std::vector<T> &Rx0,
                                        const std::vector<T> &Ry0,
                                        const std::vector<T> &AMR_X,
                                        const std::vector<T> &SMR_X,
                                        const std::vector<T> &AMR_Y,
                                        const std::vector<T> &SMR_Y,
                                        const std::vector<T> &AHE) const {
    T Rx_acc = 0.0;
    T Ry_acc = 0.0;

    for (unsigned int i = 0; i < this->layers.size(); i++) {
      const T mx = this->layers[i].mag.x;
      const T my = this->layers[i].mag.y;
      const T mz = this->layers[i].mag.z;
      const T Rx = Rx0[i] + AMR_X[i] * mx * mx + SMR_X[i] * my * my;
      const T Ry = Ry0[i] + 0.5 * AHE[i] * mz + (AMR_Y[i] + SMR_Y[i]) * mx * my;
      Rx_acc += 1. / Rx;
      Ry_acc += 1. / Ry;
    }

    return {1. / Rx_acc, 1. / Ry_acc, 0.};
  }

  /**
   * Calculate classic magnetoresistance.
   * Only for bilayer structures.
   * used when MR_MODE == CLASSIC
   * @param cosTheta: cosine between two layers.
   */
  inline T calculateMagnetoresistance(T cosTheta) const {
    return this->Rp + (((this->Rap - this->Rp) * 0.5) * (1.0 - cosTheta));
  }

  std::vector<T> getMagnetoresistance() {
    // this is classical bilayer case
    if (this->MR_mode == CLASSIC && this->layerNo == 2) {
      return {
          calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[1].mag))};
    }
    // this is the case when we use the pinning layer
    else if (this->MR_mode == CLASSIC && this->layerNo == 1) {
      return {calculateMagnetoresistance(
          c_dot<T>(layers[0].mag, layers[0].referenceLayer))};
    }
    // this is strip magnetoresistance
    else if (this->MR_mode == STRIP) {
      return stripMagnetoResistance(this->Rx0, this->Ry0, this->AMR_X,
                                    this->SMR_X, this->AMR_Y, this->SMR_Y,
                                    this->AHE);
    } else {
      throw std::runtime_error(
          "Magnetisation calculation is not supported for this structure!");
    }
  }

  std::tuple<RunnerFn<T>, SolverFn<T>, SolverMode>
  getSolver(SolverMode mode, unsigned int totalIterations) {
    SolverMode localMode = mode;
    for (auto &l : this->layers) {
      if (l.hasTemperature()) {
        // if at least one temp. driver is set
        // then use heun for consistency
        if (localMode != HEUN && localMode != EULER_HEUN) {
          std::cout << "[WARNING] Solver automatically changed to Euler Heun "
                       "for stochastic calculation."
                    << std::endl;
          localMode = EULER_HEUN;
        }
      }
      if (l.noiseParams.scaleNoise != 0) {
        // if at least one temp. driver is set
        // then use heun for consistency
        if (localMode != HEUN && localMode != EULER_HEUN) {
          std::cout << "[WARNING] Solver automatically changed to Euler Heun "
                       "for stochastic calculation."
                    << std::endl;
          localMode = EULER_HEUN;
        }
        // create a buffer
        l.createBufferedAlphaNoise(totalIterations);
      }
    }
    auto solver = &Layer<T>::rk4_step;
    if (localMode == DORMAND_PRINCE)
      solver = &Layer<T>::dormand_prince_step;
    // assign a runner function pointer from junction
    auto runner = &Junction<T>::runMultiLayerSolver;
    if (this->layerNo == 1)
      runner = &Junction<T>::runSingleLayerSolver;
    if (localMode == HEUN)
      runner = &Junction<T>::heunSolverStep;
    else if (localMode == EULER_HEUN)
      runner = &Junction<T>::eulerHeunSolverStep;

    return std::make_tuple(runner, solver, localMode);
  }

  /**
   * Main run simulation function. Use it to run the simulation.
   * @param totalTime: total time of a simulation, give it in seconds. Typical
   * length is in ~couple ns.
   * @param timeStep: the integration step of the RK45 method. Default is 1e-13
   * @param writeFrequency: how often is the log saved to? Must be no smaller
   * than `timeStep`. Default is 1e-11.
   * @param persist: whether to save to the filename specified in the Junction
   * constructor. Default is true
   * @param log: if you want some verbosity like timing the simulation. Default
   * is false
   * @param calculateEnergies: [WORK IN PROGRESS] log energy values to the log.
   * Default is false.
   * @param mode: Solver mode EULER_HEUN, RK4 or DORMAND_PRINCE
   */
  void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11,
                     bool verbose = false, bool calculateEnergies = false,
                     SolverMode mode = RK4) {

    if (timeStep > writeFrequency) {
      throw std::runtime_error(
          "The time step cannot be larger than write frequency!");
    }
    const unsigned int totalIterations =
        static_cast<unsigned int>(totalTime / timeStep);
    const unsigned int writeEvery =
        static_cast<unsigned int>(writeFrequency / timeStep);
    std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    // pick a solver based on drivers
    auto [runner, solver, solver_mode] = getSolver(mode, totalIterations);
    T t = 0.0;
    T next_write_time = 0.0;
    if (solver_mode != DORMAND_PRINCE) {
      // assign parameters
      const unsigned int totalIterations =
          static_cast<unsigned int>(totalTime / timeStep);
      const unsigned int writeEvery =
          static_cast<unsigned int>(writeFrequency / timeStep);

      for (unsigned int i = 0; i < totalIterations; i++) {
        t = i * timeStep;
        bool step_accepted = true;
        (*this.*runner)(solver, t, timeStep, step_accepted);

        if (!(i % writeEvery)) {
          logLayerParams(t, timeStep, calculateEnergies);
        }
      }
    } else {
      T t = 0.0;
      T next_write_time = 0.0;
      T current_timestep = timeStep;
      int consecutive_rejections = 0;
      const unsigned int max_rejections = 10;
      const T min_timestep_threshold = 1e-15;
      while (t < totalTime) {
        // Ensure we don't exceed totalTime
        if (t + current_timestep > totalTime) {
          current_timestep = totalTime - t;
        }
        consecutive_rejections = 0;
        // Check for NaN in timestep and fix it
        if (isnan(current_timestep) || isinf(current_timestep) ||
            current_timestep <= 0) {
          std::cout << "Warning: Invalid timestep detected (NaN/Inf). "
                       "Resetting to initial timestep at t="
                    << t << std::endl;
          current_timestep = timeStep; // Reset to initial timestep
          consecutive_rejections++;    // Count this as a rejection
        }

        // Run one step with adaptive timestep and check if accepted
        bool step_accepted = false;

        // Minimum timestep threshold (can be adjusted based on needs)

        if (layerNo > 1) {
          // For multi-layer systems
          runMultiLayerSolver(solver, t, current_timestep, step_accepted);
        } else {
          // For single layer systems
          runSingleLayerSolver(solver, t, current_timestep, step_accepted);
        }

        // Force acceptance if timestep is too small or too many consecutive
        // rejections
        if (!step_accepted && (current_timestep < min_timestep_threshold ||
                               consecutive_rejections >= max_rejections)) {
          step_accepted = true;
          consecutive_rejections = 0;
          std::cout << "Warning: Forcing step acceptance after multiple "
                       "rejections or tiny timestep at t="
                    << t << ", dt=" << current_timestep << std::endl;
        }

        // If step was accepted, update time and possibly log
        if (step_accepted) {
          t += current_timestep;
          consecutive_rejections = 0;

          // Log data at regular intervals based on writeFrequency
          if (t >= next_write_time) {
            logLayerParams(t, current_timestep, calculateEnergies);
            next_write_time =
                ((static_cast<unsigned int>(t / writeFrequency) + 1) *
                 writeFrequency);
          }
        } else {
          consecutive_rejections++;
        }
      }
    }
    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    if (verbose) {
      std::cout << "Steps in simulation: " << totalIterations << std::endl;
      std::cout << "Write every: " << writeEvery << std::endl;
      std::cout << "Simulation time = "
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin)
                       .count()
                << "[s]" << std::endl;
    }
  }
};

#endif // CORE_JUNCTION_HPP_
