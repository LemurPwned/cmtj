#ifndef CORE_LAYER_HPP_
#define CORE_LAYER_HPP_

#include "abstract.hpp"
#include "cvector.hpp"
#include "drivers.hpp"
#include "noise.hpp"
#include <functional>
#include <memory>
#include <random>
#include <string>

template <typename T> class Layer : public AbstractLayer<T> {
private:
  // CMTJ interaction drivers
  std::shared_ptr<Driver<T>> IECDriverTop =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> IECDriverBottom =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> IECQuadDriverTop =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> IECQuadDriverBottom =
      ScalarDriver<T>::getConstantDriver(0.0);
  AxialDriver<T> IDMIDriverTop = AxialDriver<T>::getVectorAxialDriver(0, 0, 0);
  AxialDriver<T> IDMIDriverBottom =
      AxialDriver<T>::getVectorAxialDriver(0, 0, 0);

  // CMTJ Torque & Field drivers
  std::shared_ptr<Driver<T>> fieldLikeTorqueDriver =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> dampingLikeTorqueDriver =
      ScalarDriver<T>::getConstantDriver(0.0);

  bool nonStochasticTempSet = false;
  bool nonStochasticOneFSet = true;
  bool temperatureSet = false;
  bool pinkNoiseSet = false;
  bool alternativeSTTSet = false;

  // Distribution for noise generation
  std::function<T()> distribution = std::bind(
      std::normal_distribution<T>(0, 1), std::mt19937(std::random_device{}()));

  CVector<T> dWn, dWn2; // one for thermal, one for OneF
  Reference referenceType = Reference::NONE;

public:
  // Added to match LLGB implementation
  T K_log = 0.0;

  struct BufferedNoiseParameters {
    T alphaNoise = 1.0;
    T scaleNoise = 0.0;
    T stdNoise = 0.0;
    Axis axis = Axis::all;
  };
  BufferedNoiseParameters noiseParams;
  std::shared_ptr<OneFNoise<T>> ofn;
  std::shared_ptr<VectorAlphaNoise<T>> bfn;
  bool includeSTT = false;
  bool includeSOT = false;
  bool dynamicSOT = false;

  T Ms = 0.0;
  T thickness = 0.0;
  T cellVolume = 0.0, cellSurface = 0.0;
  T damping = 0.1;

  CVector<T> H_log, Hoe_log, Hconst, mag, anis, referenceLayer;
  CVector<T> Hext, Hdipole, Hdemag, Hoe, HAnis, Hthermal, Hfluctuation, Hdmi,
      Hidmi;
  CVector<T> Hfl_v, Hdl_v;
  CVector<T> HIEC, HIECtop, HIECbottom;

  // Fixed constructor chain
  Layer() = default;

  explicit Layer(const std::string &id, const CVector<T> &mag,
                 const CVector<T> &anis, T Ms, T thickness, T cellSurface,
                 const std::vector<CVector<T>> &demagTensor, T damping)
      : AbstractLayer<T>(id), mag(mag), anis(anis), Ms(Ms),
        thickness(thickness), cellSurface(cellSurface), damping(damping) {
    this->demagTensor = demagTensor;
    this->cellVolume = cellSurface * thickness;
  }

  explicit Layer(const std::string &id, const CVector<T> &mag,
                 const CVector<T> &anis, T Ms, T thickness, T cellSurface,
                 const std::vector<CVector<T>> &demagTensor, T damping,
                 T fieldLikeTorque, T dampingLikeTorque) {
    Layer(id, mag, anis, Ms, thickness, cellSurface, demagTensor, damping);
    this->includeSTT = false;
    this->includeSOT = true;
    this->dynamicSOT = false;
  }

  explicit Layer(const std::string &id, const CVector<T> &mag,
                 const CVector<T> &anis, T Ms, T thickness, T cellSurface,
                 const std::vector<CVector<T>> &demagTensor, T damping,
                 T SlonczewskiSpacerLayerParameter, T beta,
                 T spinPolarisation) {
    Layer(id, mag, anis, Ms, thickness, cellSurface, demagTensor, damping);
    this->includeSTT = true;
    this->includeSOT = false;
    // Add STT parameters as needed
  }

  // Accessor methods
  CVector<T> getMagnetisation() const override { return this->mag; }

  void setMagnetisation(const CVector<T> &newMag) override {
    if (newMag.length() < 1e-10) {
      throw std::runtime_error("Magnetization vector cannot be zero!");
    }
    this->mag = newMag;
    this->mag.normalize();
  }

  // Fixed field calculation methods
  CVector<T> calculateExternalField(const T &time) {
    return this->externalFieldDriver.getCurrentAxialDrivers(time);
  }

  CVector<T> calculateHOeField(const T &time) {
    return this->HoeDriver.getCurrentAxialDrivers(time);
  }

  CVector<T> calculateHdmiField(const T &time) {
    return this->HdmiDriver.getCurrentAxialDrivers(time);
  }

  CVector<T> calculateAnisotropy(const CVector<T> &stepMag, const T &time) {
    const T K = this->anisotropyDriver->getCurrentScalarValue(time);
    const T nom = (2 * K) * c_dot<T>(this->anis, stepMag) / (this->Ms);
    return this->anis * nom;
  }
  CVector<T> calculateIEC_(const T J, const T J2, const CVector<T> &stepMag,
                           const CVector<T> &coupledMag) {
    // below an alternative method for computing J -- it's here for reference
    // only. const T nom = J / (this->Ms * this->thickness); return (coupledMag
    // - stepMag) * nom; // alternative form return (coupledMag + coupledMag * 2
    // * J2 * c_dot(coupledMag, stepMag)) * nom;
    return coupledMag * (J + 2 * J2 * c_dot(coupledMag, stepMag)) /
           (this->Ms * this->thickness);
  }

  CVector<T> calculateIEC(T time, const CVector<T> &stepMag,
                          const CVector<T> &bottom, const CVector<T> &top) {
    const T Jbottom_log = this->IECDriverBottom->getCurrentScalarValue(time);
    const T Jtop_log = this->IECDriverTop->getCurrentScalarValue(time);

    const T J2bottom_log =
        this->IECQuadDriverBottom->getCurrentScalarValue(time);
    const T J2top_log = this->IECQuadDriverTop->getCurrentScalarValue(time);

    return calculateIEC_(Jbottom_log, J2bottom_log, stepMag, bottom) +
           calculateIEC_(Jtop_log, J2top_log, stepMag, top);
  }

  CVector<T> calculateIDMI_(const CVector<T> &Dvector,
                            const CVector<T> &stepMag,
                            const CVector<T> &coupledMag) {
    // D * [(dm1/dm1x x m2) + (m1 x dm2/dm2x)]
    // dm1/dm1x x m2 = (0, -mz, my)
    // dm1/dm1y x m2 = (mz, 0, -mx)
    // dm1/dm1z x m2 = (-my, mx, 0)
    // E = D z * (m1 x m2) == D m1 (m2 x z)
    // dE/dm1 = D m2 x z
    const CVector<T> dm1crossm2 = -1.0 * c_cross<T>(Dvector, coupledMag);
    return dm1crossm2 / (this->Ms * this->thickness);
    // const CVector<T> dm1crossm2(
    //     c_dot(Dvector, CVector<T>(0, -coupledMag.z, coupledMag.y)),
    //     c_dot(Dvector, CVector<T>(coupledMag.z, 0, -coupledMag.x)),
    //     c_dot(Dvector, CVector<T>(-coupledMag.y, coupledMag.x, 0)));
    // return dm1crossm2 / (this->Ms * this->thickness);
  }

  CVector<T> calculateIDMI(T time, const CVector<T> &stepMag,
                           const CVector<T> &bottom, const CVector<T> &top) {
    return calculateIDMI_(this->IDMIDriverBottom.getCurrentAxialDrivers(time),
                          stepMag, bottom) +
           calculateIDMI_(this->IDMIDriverTop.getCurrentAxialDrivers(time),
                          stepMag, top);
  }

  // Fixed LLG calculation
  const CVector<T> calculateLLG(const T &time, const T &timeStep,
                                const CVector<T> &m,
                                const CVector<T> &bottom = CVector<T>(),
                                const CVector<T> &top = CVector<T>()) override {
    // Get effective field
    CVector<T> Heff = calculateHeff(time, timeStep, m, bottom, top);

    const CVector<T> mxh = c_cross<T>(m, Heff);
    const CVector<T> mxmxh = c_cross<T>(m, mxh);
    const T convTerm = 1 / (1 + pow(this->damping, 2)); // LLGS -> LL form
    const CVector<T> dmdt = mxh + mxmxh * this->damping;
    // Calculate LLG
    return (dmdt * -GYRO * convTerm);
  }

  // Fixed effective field calculation
  const CVector<T>
  calculateHeff(const T &time, const T &timeStep, const CVector<T> &stepMag,
                const CVector<T> &bottom, const CVector<T> &top,
                const CVector<T> &Hfluctuation = CVector<T>(),
                const CVector<T> &Hdipole = CVector<T>()) override {
    this->Hext = calculateExternalField(time);
    this->Hoe = calculateHOeField(time);

    // Calculate demagnetization field
    this->Hdemag =
        calculate_tensor_interaction(stepMag, this->demagTensor, this->Ms);

    // Other field calculations
    this->HIEC = calculateIEC(time, stepMag, bottom, top);
    this->Hidmi = calculateIDMI(time, stepMag, bottom, top);
    this->HAnis = calculateAnisotropy(stepMag, time);
    this->Hdmi = calculateHdmiField(time);

    // Get reserved interaction field
    CVector<T> HreservedInteractionField =
        this->HreservedInteractionFieldDriver.getCurrentAxialDrivers(time);

    // Sum all field contributions
    const CVector<T> Heff =
        this->Hext                   // external
        + this->HAnis                // anisotropy
        + this->HIEC                 // IEC
        + this->Hidmi                // IDMI
        + this->Hoe                  // Oersted field
        + this->Hdmi                 // regular DMI
        + Hfluctuation               // fluctuations
        - this->Hdemag               // demag -- negative contribution
        - Hdipole                    // dipole -- negative contribution
        + HreservedInteractionField; // reserved interaction field
    return Heff;
  }

  void setIECDriverBottom(const std::shared_ptr<Driver<T>> &driver) {
    this->IECDriverBottom = driver;
  }

  void setIECDriverTop(const std::shared_ptr<Driver<T>> &driver) {
    this->IECDriverTop = driver;
  }

  void setQuadIECDriverTop(const std::shared_ptr<Driver<T>> &driver) {
    this->IECQuadDriverTop = driver;
  }

  void setQuadIECDriverBottom(const std::shared_ptr<Driver<T>> &driver) {
    this->IECQuadDriverBottom = driver;
  }

  void setHoeDriver(const AxialDriver<T> &driver) { this->HoeDriver = driver; }

  void setIDMIDriver(const AxialDriver<T> &driver, bool top) {
    if (top) {
      this->IDMIDriverTop = driver;
    } else {
      this->IDMIDriverBottom = driver;
    }
  }

  void setFieldLikeTorqueDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->fieldLikeTorqueDriver = driver;
  }

  void setDampingLikeTorqueDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->dampingLikeTorqueDriver = driver;
  }

  // RK4 step implementation (simplified for now)
  void rk4_step(const T &time, const T &timeStep, const CVector<T> &bottom,
                const CVector<T> &top) {
    // Basic RK4 implementation for testing
    CVector<T> m_t = this->mag;
    CVector<T> k1 = calculateLLG(time, timeStep, m_t, bottom, top) * timeStep;
    CVector<T> k2 =
        calculateLLG(time + timeStep / 2, timeStep, m_t + k1 / 2, bottom, top) *
        timeStep;
    CVector<T> k3 =
        calculateLLG(time + timeStep / 2, timeStep, m_t + k2 / 2, bottom, top) *
        timeStep;
    CVector<T> k4 =
        calculateLLG(time + timeStep, timeStep, m_t + k3, bottom, top) *
        timeStep;

    m_t += (k1 + k2 * 2 + k3 * 2 + k4) / 6;
    m_t.normalize();
    this->mag = m_t;
    if (isnan(this->mag.x)) {
      throw std::runtime_error("NAN magnetisation");
    }
  }
};

#endif // CORE_LAYER_HPP_
