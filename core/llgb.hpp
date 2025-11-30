#include "constants.hpp"
#include "junction.hpp"
#include <map>
#include <tuple>

namespace LLGB {

template <typename T = double> T langevin(T x) {
  return (1.0 / tanh(x)) - (1.0 / x);
}

template <typename T = double> T langevinDerivative(T x) {
  return (-1.0 / pow(sinh(x), 2)) + (1. / pow(x, 2));
}

template <typename T = double>
std::tuple<T, T> MFAWeissCurie(T est, T temp, T J0, T relax = 0.2, T tol = 1e-6,
                               unsigned int maxIter = 1000) {
  /**
      This function solves the self-consistent Curie-Weiss equation in MFA
      The equation is given by:
      x = L(beta * J0 * x)
      where beta = 1/(k * T) and J0 is the exchange constant.
      The function returns the solution and the error.
      E.g. for FePt ~ 3.051823739e-20 J => Tc ~ 760 K

      @param est: initial guess
      @param temp: temperature
      @param J0: exchange constant
      @param relax: relaxation factor
      @param tol: tolerance
      @param maxIter: maximum number of iterations
  **/
  T beta = (1.0 / (BOLTZMANN_CONST * temp));
  T err = 0;
  for (unsigned int i = 0; i < maxIter; i++) {
    T xNext = langevin(beta * J0 * est);
    err = abs(xNext - est);
    if (err < tol) {
      return std::make_tuple(xNext, err);
    }
    est = relax * xNext + (1 - relax) * est;
  }
  return std::make_tuple(est, err);
}
} // namespace LLGB

template <typename T = double> class LLGBLayer {
protected:
  ScalarDriver<T> temperatureDriver;
  ScalarDriver<T> anisotropyDriver;
  AxialDriver<T> externalFieldDriver;
  // the distribution is binded (bound?) for faster generation
  // we need two distributions for the two types of noise in the LLB
  std::function<T()> distributionA = std::bind(
      std::normal_distribution<T>(0, 1), std::mt19937(std::random_device{}()));
  std::function<T()> distributionB = std::bind(
      std::normal_distribution<T>(0, 1), std::mt19937(std::random_device{}()));

public:
  std::string id;
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
      : id(id), mag(mag), anis(anis), Ms(Ms), thickness(thickness),
        surface(surface), demagTensor(demagTensor), damping(damping), Tc(Tc),
        susceptibility(susceptibility), me(me) {
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

  T getAlphaParallel(T &time) {
    const T temp = this->temperatureDriver.getCurrentScalarValue(time);
    this->alpha_par_log = this->damping * (temp / this->Tc) * (2. / 3.);
    return this->alpha_par_log;
  }

  T getAlphaPerpendicular(T &time) {
    const T temp = this->temperatureDriver.getCurrentScalarValue(time);
    const T ratio = temp / this->Tc;
    if (temp >= this->Tc) {
      this->alpha_perp_log = this->damping * ratio * (2. / 3.);
    } else {
      this->alpha_perp_log = this->damping * (1. - ratio * (1. / 3.0));
    }
    return this->alpha_perp_log;
  }

  CVector<T> getLongitudinal(T &time, const CVector<T> &m) {
    const T temp = this->temperatureDriver.getCurrentScalarValue(time);
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

  CVector<T> getAnisotropyField(T &time, const CVector<T> &m) {
    return (-1. / this->anisotropyDriver.getCurrentScalarValue(time)) *
           CVector<T>(m[0], m[1], 0);
  }

  CVector<T> calculateAnisotropy(const CVector<T> &stepMag, T &time) {
    this->K_log = this->anisotropyDriver.getCurrentScalarValue(time);
    const T nom = this->K_log * c_dot<T>(this->anis, stepMag);
    return this->anis * nom;
  }

  const CVector<T> calculateHeff(T time, const CVector<T> &m) {
    // this anisotropy is a bit different than in the LLG
    // const CVector<T> anisotropy = this->getAnisotropyField(time, m);
    const CVector<T> anisotropy = this->calculateAnisotropy(m, time);
    const CVector<T> hext =
        this->externalFieldDriver.getCurrentAxialDrivers(time);
    const CVector<T> longField = this->getLongitudinal(time, m);
    return anisotropy + hext + longField;
  }

  CVector<T> calculateLLG(const T &time, const T &timeStep,
                          const CVector<T> &m) {
    const CVector<T> heff = this->calculateHeff(time, m);
    return solveLLG(time, timeStep, m, heff);
  }

  const CVector<T> solveLLG(T time, T timeStep, const CVector<T> &m,
                            const CVector<T> &heff) {
    T_log = this->temperatureDriver.getCurrentScalarValue(time);
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

  CVector<T> nonadiabaticThermalField(T time, T timestep) {
    const T temp = this->temperatureDriver.getCurrentScalarValue(time);
    const T varianceDev =
        (2 * BOLTZMANN_CONST * temp *
         (this->getAlphaPerpendicular(time) - this->getAlphaParallel(time))) /
        (this->volume * this->Ms * GYRO *
         pow(this->getAlphaPerpendicular(time), 2));
    return 0 * sqrt(varianceDev) * CVector<T>(this->distributionA);
  }

  CVector<T> adiabaticThermalField(T time, T timestep) {
    const T temp = this->temperatureDriver.getCurrentScalarValue(time);
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

  // setters

  void setEquilibriumMagnetisation(const T &me) { this->me = me; }

  void setSusceptibility(const T &susceptibility) {
    this->susceptibility = susceptibility;
  }

  void setTemperatureDriver(const ScalarDriver<T> &driver) {
    this->temperatureDriver = driver;
  }

  void setExternalFieldDriver(const AxialDriver<T> &driver) {
    this->externalFieldDriver = driver;
  }

  void setAnisotropyDriver(const ScalarDriver<T> &driver) {
    this->anisotropyDriver = driver;
  }
};

template <typename T = double> class LLGBJunction {
private:
  // friend class LLGBLayer<T>;
  const std::vector<std::string> vectorNames = {"x", "y", "z"};
  std::vector<LLGBLayer<T>> layers;
  std::unordered_map<std::string, std::vector<T>> log;
  unsigned int logLength = 0;
  unsigned int layerNo = 0;
  T time = 0;

public:
  explicit LLGBJunction(const std::vector<LLGBLayer<T>> &layers) {
    this->layers = layers;
    this->layerNo = layers.size();
  }

  typedef void (LLGBLayer<T>::*scalarDriverSetter)(
      const ScalarDriver<T> &driver);
  typedef void (LLGBLayer<T>::*axialDriverSetter)(const AxialDriver<T> &driver);
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
  void setLayerTemperatureDriver(const std::string &layerID,
                                 const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &LLGBLayer<T>::setTemperatureDriver, driver);
  }
  void setLayerExternalFieldDriver(const std::string &layerID,
                                   const AxialDriver<T> &driver) {
    axiallayerSetter(layerID, &LLGBLayer<T>::setExternalFieldDriver, driver);
  }
  void setLayerAnisotropyDriver(const std::string &layerID,
                                const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerID, &LLGBLayer<T>::setAnisotropyDriver, driver);
  }
  void setLayerEquilibriumMagnetisation(const std::string &layerID,
                                        const T &me) {
    bool found = false;
    for (auto &l : this->layers) {
      if (l.id == layerID || layerID == "all") {
        l.setEquilibriumMagnetisation(me);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  void setLayerSusceptibility(const std::string &layerID,
                              const T &susceptibility) {
    bool found = false;
    for (auto &l : this->layers) {
      if (l.id == layerID || layerID == "all") {
        l.setSusceptibility(susceptibility);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  void heunSolverStep(const T &t, const T &timeStep) {
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
    std::vector<CVector<T>> fn(this->layerNo, CVector<T>());
    std::vector<CVector<T>> gn(this->layerNo, CVector<T>());
    std::vector<CVector<T>> nonAdiabatic(this->layerNo, CVector<T>());
    std::vector<CVector<T>> adiabatic(this->layerNo, CVector<T>());
    std::vector<CVector<T>> mNext(this->layerNo, CVector<T>());
    // first approximation

    // make sure that
    // 1. Thermal field is added if needed
    // 2. One/f noise is added if needed
    // 3. The timestep is correctly multiplied

    for (unsigned int i = 0; i < this->layerNo; i++) {
      fn[i] = this->layers[i].calculateLLG(t, timeStep, this->layers[i].mag);

      // draw the noise for each layer, dW
      nonAdiabatic[i] = this->layers[i].nonadiabaticThermalField(t, timeStep);
      adiabatic[i] = this->layers[i].adiabaticThermalField(t, timeStep);
      gn[i] = this->layers[i].stochasticTorque(this->layers[i].mag, t,
                                               nonAdiabatic[i], adiabatic[i]);

      mNext[i] =
          this->layers[i].mag + fn[i] * timeStep + gn[i] * sqrt(timeStep);
    }
    // second approximation
    for (unsigned int i = 0; i < this->layerNo; i++) {
      // first approximation is already multiplied by timeStep
      this->layers[i].mag =
          this->layers[i].mag +
          0.5 * timeStep *
              (fn[i] +
               this->layers[i].calculateLLG(t + timeStep, timeStep, mNext[i])) +
          0.5 *
              (gn[i] + this->layers[i].stochasticTorque(mNext[i], t + timeStep,
                                                        nonAdiabatic[i],
                                                        adiabatic[i])) *
              sqrt(timeStep);
      // normalise only in classical
      // this->layers[i].mag.normalize(); // LLB doesn't normalise
    }
  }
  void eulerHeunSolverStep(const T &t, const T &timeStep) {
    /*
        Euler Heun method (stochastic heun)

        y_np = y + g(y,t,dW)*dt
        g_sp = g(y_np,t+1,dW)
        y(t+1) = y + dt*f(y,t) + .5*(g(y,t,dW)+g_sp)*sqrt(dt)

        with f being the non-stochastic part and g the stochastic part
    */
    // draw the noise for each layer, dW
    std::vector<CVector<T>> mPrime(this->layerNo, CVector<T>());
    for (unsigned int i = 0; i < this->layerNo; i++) {
      // todo: after you're done, double check the thermal magnitude and dt
      // scaling there
      const CVector<T> nonAdiabaticTorque =
          this->layers[i].nonadiabaticThermalField(t, timeStep);
      const CVector<T> adiabaticTorque =
          this->layers[i].adiabaticThermalField(t, timeStep);

      const CVector<T> fnApprox =
          this->layers[i].calculateLLG(t, timeStep, this->layers[i].mag);
      const CVector<T> gnApprox = this->layers[i].stochasticTorque(
          this->layers[i].mag, t, nonAdiabaticTorque, adiabaticTorque);

      // theoretically we have 2 options
      // 1. calculate only the stochastic part with the second approximation
      // 2. calculate the second approximation of m with the stochastic and
      // non-stochastic
      //    part and then use if for torque est.
      const CVector<T> mNext = this->layers[i].mag + gnApprox * sqrt(timeStep);
      const CVector<T> gnPrimeApprox = this->layers[i].stochasticTorque(
          mNext, t + timeStep, nonAdiabaticTorque, adiabaticTorque);
      mPrime[i] = this->layers[i].mag + fnApprox * timeStep +
                  0.5 * (gnApprox + gnPrimeApprox) * sqrt(timeStep);
    }

    for (unsigned int i = 0; i < this->layerNo; i++) {
      this->layers[i].mag = mPrime[i];
      // this->layers[i].mag.normalize(); // LLB doesn't normalise
    }
  }

  typedef void (LLGBJunction<T>::*runnerFn)(const T &t, const T &timeStep);
  std::tuple<runnerFn, SolverMode> getSolver(SolverMode mode) {
    auto runner = &LLGBJunction<T>::heunSolverStep;
    if (mode == HEUN)
      runner = &LLGBJunction<T>::heunSolverStep;
    else if (mode == EULER_HEUN)
      runner = &LLGBJunction<T>::eulerHeunSolverStep;
    else
      throw std::runtime_error("The solver mode is not supported!");
    return std::make_tuple(runner, mode);
  }

  /**
   * @brief Log computed layer parameters.
   * This function logs all the necessayr parameters of the layers.
   * @param t: current time
   * @param timeStep: timeStep of the simulation (unsued for now)
   * @param calculateEnergies: if true, also include fields for energy
   * computation.
   */
  void logLayerParams(T t, const T timeStep) {
    for (const auto &layer : this->layers) {
      const std::string lId = layer.id;
      // always save magnetisation
      for (int i = 0; i < 3; i++) {
        this->log[lId + "_m" + vectorNames[i]].emplace_back(layer.mag[i]);
      }
      this->log[lId + "_alpha_parallel"].emplace_back(layer.alpha_par_log);
      this->log[lId + "_alpha_perpendicular"].emplace_back(
          layer.alpha_perp_log);
      this->log[lId + "_K"].emplace_back(layer.K_log);
      this->log[lId + "_T"].emplace_back(layer.T_log);
      this->log[lId + "_me"].emplace_back(layer.me);
      this->log[lId + "_Xpar"].emplace_back(layer.susceptibility);
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
   * Clears the simulation log.
   **/
  void clearLog() {
    this->log.clear();
    this->logLength = 0;
    this->time = 0;
  }

  std::unordered_map<std::string, std::vector<T>> &getLog() {
    return this->log;
  }

  /**
   * Main run simulation function. Use it to run the simulation.
   * @param totalTime: total time of a simulation, give it in seconds. Typical
   * length is in ~couple ns.
   * @param timeStep: the integration step of the RK45 method. Default is 1e-13
   * @param writeFrequency: how often is the log saved to? Must be no smaller
   * than `timeStep`. Default is 1e-11.
   * @param log: if you want some verbosity like timing the simulation. Default
   * is false
   * @param mode: Solver mode EULER_HEUN, RK4 or DORMAND_PRICE
   */
  void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-13,
                     bool log = false, SolverMode mode = HEUN)

  {
    if (timeStep > writeFrequency) {
      throw std::runtime_error(
          "The time step cannot be larger than write frequency!");
    }
    const unsigned int totalIterations = (int)(totalTime / timeStep);
    const unsigned int writeEvery = (int)(writeFrequency / timeStep);
    std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    // pick a solver based on drivers
    auto [runner, _] = getSolver(mode);

    for (unsigned int i = 0; i < totalIterations; i++) {
      this->time += timeStep;
      (*this.*runner)(this->time, timeStep);

      if (!(i % writeEvery)) {
        logLayerParams(this->time, timeStep);
      }
    }
    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    if (log) {
      std::cout << "Steps in simulation: " << totalIterations << std::endl;
      std::cout << "Write every: " << writeEvery << std::endl;
      std::cout << "Simulation time = "
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin)
                       .count()
                << "[s]" << std::endl;
    }
  }
};
