#ifndef CORE_ABSTRACT_HPP
#define CORE_ABSTRACT_HPP
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random> // for mt19937, normal_distribution
#include <string>
#include <utility>
#include <vector>

#include "cvector.hpp"
#include "driver.hpp"
#include "noise.hpp" // for OneFNoise

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 220880.0 // rad/Ts converted to m/As
#define TtoAm 795774.715459
#define HBAR 6.62607015e-34 / (2. * M_PI)
#define ELECTRON_CHARGE 1.60217662e-19
#define BOLTZMANN_CONST 1.380649e-23

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
  CVector<T> res(
      tensor[0][0] * m[0] + tensor[0][1] * m[1] + tensor[0][2] * m[2],
      tensor[1][0] * m[0] + tensor[1][1] * m[1] + tensor[1][2] * m[2],
      tensor[2][0] * m[0] + tensor[2][1] * m[1] + tensor[2][2] * m[2]);
  return res * (Ms / MAGNETIC_PERMEABILITY);
}

template <typename T>
inline CVector<T> calculate_tensor_interaction(
    const CVector<T> &m, const std::array<CVector<T>, 3> &tensor, const T &Ms) {
  CVector<T> res(
      tensor[0][0] * m[0] + tensor[0][1] * m[1] + tensor[0][2] * m[2],
      tensor[1][0] * m[0] + tensor[1][1] * m[1] + tensor[1][2] * m[2],
      tensor[2][0] * m[0] + tensor[2][1] * m[1] + tensor[2][2] * m[2]);
  return res * (Ms / MAGNETIC_PERMEABILITY);
}

template <typename T>
inline CVector<T> c_cross(const CVector<T> &a, const CVector<T> &b) {
  CVector<T> res(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                 a[0] * b[1] - a[1] * b[0]);

  return res;
}

template <typename T> inline T c_dot(const CVector<T> &a, const CVector<T> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
enum Reference { NONE = 0, FIXED, TOP, BOTTOM };
enum SolverMode { EULER_HEUN = 0, RK4 = 1, DORMAND_PRICE = 2, HEUN = 3 };
enum MRmode { MR_NONE = 0, CLASSIC = 1, STRIP = 2 };

template <typename T> class AbstractLayer {

protected:
  std::string id;

  Reference referenceType = Reference::NONE;
  std::shared_ptr<Driver<T>> anisotropyDriver =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> temperatureDriver =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> currentDriver =
      ScalarDriver<T>::getConstantDriver(0.0);
  AxialDriver<T> externalFieldDriver =
      AxialDriver<T>::getVectorAxialDriver(0, 0, 0);
  AxialDriver<T> HoeDriver = AxialDriver<T>::getVectorAxialDriver(0, 0, 0);
  AxialDriver<T> HdmiDriver = AxialDriver<T>::getVectorAxialDriver(0, 0, 0);
  AxialDriver<T> HreservedInteractionFieldDriver =
      AxialDriver<T>::getVectorAxialDriver(0, 0, 0);
  std::vector<CVector<T>> demagTensor = {
      CVector<T>(0, 0, 0), CVector<T>(0, 0, 0), CVector<T>(0, 0, 0)};

public:
  virtual ~AbstractLayer() = default;
  AbstractLayer(const std::string &id) : id(id) {};
  // Pure virtual functions that must be implemented
  virtual void setMagnetisation(const CVector<T> &newMag) = 0;
  virtual CVector<T> getMagnetisation() const = 0;
  void setReferenceLayer(Reference reference) {
    this->referenceType = reference;
  }
  Reference getReferenceType() const { return this->referenceType; }
  std::string getId() { return this->id; }

  // Common driver setters that must be implemented
  void setTemperatureDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->temperatureDriver = driver;
  }
  void setCurrentDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->currentDriver = driver;
  }
  void setAnisotropyDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->anisotropyDriver = driver;
  }
  void setExternalFieldDriver(const AxialDriver<T> &driver) {
    this->externalFieldDriver = driver;
  }
  void setOerstedFieldDriver(const AxialDriver<T> &driver) {
    this->HoeDriver = driver;
  }
  void setHdmiDriver(const AxialDriver<T> &driver) {
    this->HdmiDriver = driver;
  }
  void setReservedInteractionFieldDriver(const AxialDriver<T> &driver) {
    this->HreservedInteractionFieldDriver = driver;
  }
  // Solver methods that must be implemented
  virtual const CVector<T>
  calculateHeff(const T &time, const T &timeStep, const CVector<T> &stepMag,
                const CVector<T> &bottom, const CVector<T> &top,
                const CVector<T> &Hfluctuation = CVector<T>(),

                const CVector<T> &Hdipole = CVector<T>()) = 0;

  virtual const CVector<T>
  calculateLLG(const T &time, const T &timeStep, const CVector<T> &m,
               const CVector<T> &bottom = CVector<T>(),
               const CVector<T> &top = CVector<T>()) = 0;

  // Type definitions for driver setters
  typedef void (AbstractLayer<T>::*scalarDriverSetter)(const Driver<T> &);
  typedef void (AbstractLayer<T>::*axialDriverSetter)(const AxialDriver<T> &);
};

template <typename T> class AbstractJunction {
protected:
  // Common member variables
  std::vector<std::shared_ptr<AbstractLayer<T>>>
      layers; // Changed to store pointers to AbstractLayer
  MRmode MR_mode = MRmode::MR_NONE;
  unsigned int logLength = 0;
  unsigned int layerNo = 0;
  std::string Rtag = "R";
  std::unordered_map<std::string, std::vector<T>> log;

  const std::vector<std::string> vectorNames = {"x", "y", "z"};
  T time = 0;

public:
  AbstractJunction(const std::vector<AbstractLayer<T>> &layers) {
    this->layers = layers;
    this->layerNo = layers.size();
  }

  // Pure virtual functions for layer management
  const std::vector<std::string> getLayerIds() const {
    std::vector<std::string> ids;
    for (const auto &layer : layers) {
      ids.push_back(layer.id);
    }
    return ids;
  }
  unsigned int getLayerCount() const { return this->layerNo; }

  // Pure virtual functions for driver management
  virtual void
  scalarlayerSetter(const std::string &layerID,
                    typename AbstractLayer<T>::scalarDriverSetter functor,
                    Driver<T> driver) = 0;

  typedef void (AbstractLayer<T>::*scalarDriverSetter)(const Driver<T> &driver);
  typedef void (AbstractLayer<T>::*axialDriverSetter)(
      const AxialDriver<T> &driver);
  void scalarlayerSetter(const std::string &layerID, scalarDriverSetter functor,
                         Driver<T> driver) {
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

  virtual void setCouplingDriver(
      const std::string &bottomLayer, const std::string &topLayer,
      const Driver<T> &driver,
      typename AbstractLayer<T>::scalarDriverSetter setDriverFuncTop,
      typename AbstractLayer<T>::scalarDriverSetter setDriverFuncBottom) = 0;

  // Pure virtual functions for magnetoresistance
  virtual std::vector<T> getMagnetoresistance() = 0;

  virtual void calculateLLG(const T &t, CVector<T> &m, const CVector<T> &bottom,
                            const CVector<T> &top, const T timeStep) = 0;

  void heunSolverStep(const T &t, const T &timeStep,
                      const bool normalise = true) {
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

      const CVector<T> bottom =
          (i == 0) ? CVector<T>() : this->layers[i - 1].mag;
      const CVector<T> top =
          (i == this->layerNo - 1) ? CVector<T>() : this->layers[i + 1].mag;

      fn[i] = this->layers[i].calculateLLG(t, this->layers[i].mag, bottom, top,
                                           timeStep);

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
      if (normalise) {
        this->layers[i].mag.normalize(); // LLB doesn't normalise
      }
    }
  }

  void eulerHeunSolverStep(const T &t, const T &timeStep,
                           const bool normalise = true) {
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
      if (normalise) {
        this->layers[i].mag.normalize(); // LLB doesn't normalise
      }
    }
  }

  void commonHeunSolverStep(const T &t, const T &timeStep,
                            const bool normalise = true) {
    /*
    Euler Heun method (stochastic heun)

    y_np = y + g(y,t,dW)*dt
    g_sp = g(y_np,t+1,dW)
    y(t+1) = y + dt*f(y,t) + .5*(g(y,t,dW)+g_sp)*sqrt(dt)

    with f being the non-stochastic part and g the stochastic part
    */
    std::vector<CVector<T>> mPrime(this->layerNo, CVector<T>());
    for (unsigned int i = 0; i < this->layerNo; i++) {
    }
  }

  typedef void (AbstractJunction<T>::*runnerFn)(const T &t, const T &timeStep);
  std::tuple<runnerFn, SolverMode> getSolver(SolverMode mode) {
    auto runner = &AbstractJunction<T>::heunSolverStep;
    if (mode == HEUN)
      runner = &AbstractJunction<T>::heunSolverStep;
    else if (mode == EULER_HEUN)
      runner = &AbstractJunction<T>::eulerHeunSolverStep;
    else
      throw std::runtime_error("The solver mode is not supported!");
    return std::make_tuple(runner, mode);
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

  // log management
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
};

#endif // CORE_ABSTRACT_HPP
