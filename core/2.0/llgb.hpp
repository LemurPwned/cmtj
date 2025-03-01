#ifndef CORE_JUNCTION_IMPL_HPP
#define CORE_JUNCTION_IMPL_HPP

#include "abstract.hpp"

template <typename T = double> class LLGBLayer : public AbstractLayer<T> {
protected:
  // the distribution is binded (bound?) for faster generation
  // we need two distributions for the two types of noise in the LLB
  std::function<T()> distributionA = std::bind(
      std::normal_distribution<T>(0, 1), std::mt19937(std::random_device{}()));
  std::function<T()> distributionB = std::bind(
      std::normal_distribution<T>(0, 1), std::mt19937(std::random_device{}()));

public:
  CVector<T> mag;
  CVector<T> anis;
  T Ms;
  T thickness, surface, volume;
  std::vector<CVector<T>> demagTensor;
  T damping;
  T Tc;
  T susceptibility;
  T me;
  T alpha_perp_log, alpha_par_log;
  T K_log = 0;
  T T_log = 0;

  /// @brief
  /// @param id
  /// @param mag
  /// @param anis
  /// @param Ms
  /// @param thickness
  /// @param surface
  /// @param demagTensor
  /// @param damping
  /// @param Tc
  /// @param susceptibility
  /// @param me
  LLGBLayer(const std::string &id, const CVector<T> &mag,
            const CVector<T> &anis, T Ms, T thickness, T surface,
            const std::vector<CVector<T>> &demagTensor, T damping, T Tc,
            T susceptibility, T me)
      : AbstractLayer<T>(id), mag(mag), anis(anis), Ms(Ms),
        thickness(thickness), surface(surface), demagTensor(demagTensor),
        damping(damping), Tc(Tc), susceptibility(susceptibility), me(me) {
    this->volume = this->surface * this->thickness;
    if (this->volume == 0) {
      throw std::runtime_error("The volume of the LLGB layer cannot be 0!");
    }
    if (mag.length() == 0) {
      throw std::runtime_error(
          "Initial magnetisation was set to a zero vector!");
    }
    if (anis.length() == 0) {
      throw std::runtime_error("Anisotropy was set to a zero vector!");
    }
  }

  CVector<T> getMagnetisation() const override { return this->mag; }

  void setMagnetisation(const CVector<T> &newMag) override {
    if (newMag.length() == 0) {
      throw std::runtime_error(
          "Initial magnetisation was set to a zero vector!");
    }
    this->mag = newMag;
    this->mag.normalize();
  }

  T getAlphaParallel(const T &time) {
    const T temp = this->temperatureDriver->getCurrentScalarValue(time);
    this->alpha_par_log = this->damping * (temp / this->Tc) * (2. / 3.);
    return this->alpha_par_log;
  }

  T getAlphaPerpendicular(const T &time) {
    const T temp = this->temperatureDriver->getCurrentScalarValue(time);
    const T ratio = temp / this->Tc;
    if (temp >= this->Tc) {
      this->alpha_perp_log = this->damping * ratio * (2. / 3.);
    } else {
      this->alpha_perp_log = this->damping * (1. - ratio * (1. / 3.0));
    }
    return this->alpha_perp_log;
  }

  CVector<T> getLongitudinal(const T &time, const CVector<T> &m) {
    const T temp = this->temperatureDriver->getCurrentScalarValue(time);
    const T ratio_susc = 1. / (2. * this->susceptibility);
    const T m2 = pow(m.length(), 2);
    if (temp <= this->Tc) {
      const T ratio_m = m2 / pow(this->me, 2);
      return ratio_susc * (1. - ratio_m) * m;
    }
    const T ratio_T = (this->Tc / (this->Tc - temp));
    const T ratio_T_adj = (3. / 5.) * ratio_T * m2 - 1.;
    // this is given by some other paper
    const T ration_T_alt =
        (1. + (3. / 5.) * (this->Tc / (temp - this->Tc)) * m2);
    return -(1. / this->susceptibility) * ration_T_alt * m;
    // return ratio_susc * ratio_T_adj * m;
  }

  CVector<T> getAnisotropyField(const T &time, const CVector<T> &m) {
    return (-1. / this->anisotropyDriver->getCurrentScalarValue(time)) *
           CVector<T>(m[0], m[1], 0);
  }

  CVector<T> calculateAnisotropy(const CVector<T> &stepMag, const T &time) {
    this->K_log = this->anisotropyDriver->getCurrentScalarValue(time);
    const T nom = this->K_log * c_dot<T>(this->anis, stepMag);
    return this->anis * nom;
  }

  const CVector<T>
  calculateHeff(const T &time, const T &timeStep, const CVector<T> &stepMag,
                const CVector<T> &bottom, const CVector<T> &top,
                const CVector<T> &Hfluctuation = CVector<T>(),
                const CVector<T> &Hdipole = CVector<T>()) override {
    // this anisotropy is a bit different than in the LLG
    // const CVector<T> anisotropy = this->getAnisotropyField(time, m);
    const CVector<T> anisotropy = this->calculateAnisotropy(stepMag, time);
    const CVector<T> hext =
        this->externalFieldDriver.getCurrentAxialDrivers(time);
    const CVector<T> longField = this->getLongitudinal(time, stepMag);
    return anisotropy + hext + longField;
  }

  const CVector<T> calculateLLG(const T &time, const T &timeStep,
                                const CVector<T> &m,
                                const CVector<T> &bottom = CVector<T>(),
                                const CVector<T> &top = CVector<T>()) override {
    const CVector<T> heff = this->calculateHeff(time, timeStep, m, bottom, top);
    return solveLLG(time, timeStep, m, heff);
  }

  const CVector<T> solveLLG(const T &time, const T &timeStep,
                            const CVector<T> &m, const CVector<T> &heff) {
    T_log = this->temperatureDriver->getCurrentScalarValue(time);
    const CVector<T> mxh = c_cross<T>(m, heff);
    const CVector<T> mxmxh = c_cross<T>(m, mxh);
    const CVector<T> llbTerm = c_dot(m, heff) * m;
    const T inv_mlen = pow(1. / m.length(), 2);
    const T gamma_p = GYRO / (1 + pow(this->damping, 2)); // LLGS -> LL form
    const CVector<T> dmdt = -1 * mxh -
                            getAlphaPerpendicular(time) * mxmxh * inv_mlen +
                            llbTerm * getAlphaParallel(time) * inv_mlen;
    return gamma_p * dmdt;
  }

  CVector<T> nonadiabaticThermalField(const T &time, const T &timestep) {
    const T temp = this->temperatureDriver->getCurrentScalarValue(time);
    const T varianceDev =
        (2 * BOLTZMANN_CONST * temp *
         (this->getAlphaPerpendicular(time) - this->getAlphaParallel(time))) /
        (this->volume * this->Ms * GYRO *
         pow(this->getAlphaPerpendicular(time), 2));
    return 0 * sqrt(varianceDev) * CVector<T>(this->distributionA);
  }

  CVector<T> adiabaticThermalField(const T &time, const T &timestep) {
    const T temp = this->temperatureDriver->getCurrentScalarValue(time);
    // GYRO multiplies in the stochasticTorque for consistency
    const T varianceDev =
        (2 * BOLTZMANN_CONST * temp * this->getAlphaParallel(time)) /
        (GYRO * this->volume * this->Ms);
    return 0 * sqrt(varianceDev) * CVector<T>(this->distributionB);
  }

  CVector<T> stochasticTorque(const CVector<T> &m, T time,
                              const CVector<T> &nonAdiabatic,
                              const CVector<T> &adiabatic) {
    /*
        This formulation follows:
        Axitia, 2015, Fundamentals and applications of the Landau–Lifshitz–Bloch
       equation Evans, 2012, Stochastic form of the Landau-Lifshitz-Bloch
       equation Read Evans to understand the switch.

        This is more correct than stochasticTorqueOld, and used more recently
    */
    const T inv_mlen = pow(1. / m.length(), 2);
    const T gamma_p = GYRO / (1 + pow(this->damping, 2)); // LLGS -> LL form
    const CVector<T> nonAdiabaticTerm =
        c_cross<T>(m, c_cross<T>(m, nonAdiabatic));
    return -gamma_p * inv_mlen * getAlphaPerpendicular(time) *
               nonAdiabaticTerm +
           gamma_p * adiabatic;
  }

  CVector<T> stochasticTorqueOld(const CVector<T> &m, T time,
                                 const CVector<T> &nonAdiabatic,
                                 const CVector<T> &adiabatic) {
    /*
        This formulation follows:
        Atxitia, 2007, Micromagnetic modeling of laser-induced magnetization
       dynamics using the Landau-Lifshitz-Bloch equation And classical: Garanin,
       2004, Thermal fluctuations and longitudinal relaxation of single-domain
       magnetic particles at elevated temperatures
    */
    const T inv_mlen = pow(1. / m.length(), 2);
    const T gamma_p = GYRO / (1 + pow(this->damping, 2)); // LLGS -> LL form
    const CVector<T> nonAdiabaticTerm =
        c_cross<T>(m, c_cross<T>(m, nonAdiabatic));
    const CVector<T> adiabaticTerm = c_dot(m, adiabatic) * m;
    return gamma_p * inv_mlen *
           (adiabaticTerm * getAlphaParallel(time) -
            nonAdiabaticTerm * getAlphaPerpendicular(time));
  }

  void setEquilibriumMagnetisation(const T &me) { this->me = me; }

  void setSusceptibility(const T &susceptibility) {
    this->susceptibility = susceptibility;
  }
};

template <typename T = double>
class LLGBJunctionNew : public AbstractJunction<T> {
public:
  explicit LLGBJunctionNew(const std::vector<LLGBLayer<T>> &layers)
      : AbstractJunction<T>(layers) {}

  void setCouplingDriver(
      const std::string &bottomLayer, const std::string &topLayer,
      const ScalarDriver<T> &driver,
      typename AbstractLayer<T>::scalarDriverSetter setDriverFuncTop,
      typename AbstractLayer<T>::scalarDriverSetter setDriverFuncBottom)
      override {
    // LLGB implementation doesn't support coupling drivers directly
    throw std::runtime_error(
        "Coupling drivers not supported in LLGB implementation");
  }

  std::vector<T> getMagnetoresistance() override {
    // Implement magnetoresistance calculation for LLGB
    return std::vector<T>(); // Return empty vector for now
  }

  void heunSolverStep(const T &t, const T &timeStep,
                      const bool normalise = false) override {
    AbstractJunction<T>::heunSolverStep(t, timeStep, false);
  }

  void eulerHeunSolverStep(const T &t, const T &timeStep,
                           const bool normalise = false) override {
    AbstractJunction<T>::eulerHeunSolverStep(t, timeStep, false);
  }
};

#endif // CORE_JUNCTION_IMPL_HPP
