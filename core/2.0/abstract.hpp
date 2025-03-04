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
#include <unordered_set>
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

template <typename T> struct BufferedNoiseParameters {
  T alphaNoise = 1.0;
  T scaleNoise = 0.0;
  T stdNoise = 0.0;
  Axis axis = Axis::all;
};

template <typename T> class AbstractLayer {

protected:
  std::string id;
  Reference referenceType = Reference::NONE;
  CVector<T> referenceLayer = CVector<T>(0, 0, 0);

  std::shared_ptr<Driver<T>> anisotropyDriver =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> temperatureDriver =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<Driver<T>> currentDriver =
      ScalarDriver<T>::getConstantDriver(0.0);
  std::shared_ptr<AxialDriver<T>> externalFieldDriver =
      std::make_shared<AxialDriver<T>>(
          AxialDriver<T>::getVectorAxialDriver(0, 0, 0));
  std::shared_ptr<AxialDriver<T>> HoeDriver = std::make_shared<AxialDriver<T>>(
      AxialDriver<T>::getVectorAxialDriver(0, 0, 0));
  std::shared_ptr<AxialDriver<T>> HdmiDriver = std::make_shared<AxialDriver<T>>(
      AxialDriver<T>::getVectorAxialDriver(0, 0, 0));
  std::shared_ptr<AxialDriver<T>> HreservedInteractionFieldDriver =
      std::make_shared<AxialDriver<T>>(
          AxialDriver<T>::getVectorAxialDriver(0, 0, 0));
  BufferedNoiseParameters<T> noiseParams;

public:
  std::vector<CVector<T>> demagTensor = {
      CVector<T>(0, 0, 0), CVector<T>(0, 0, 0), CVector<T>(0, 0, 0)};

  bool isStochastic = false;
  CVector<T> mag = CVector<T>(0, 0, 0);
  virtual ~AbstractLayer() = default;
  AbstractLayer(const std::string &id) : id(id) {};
  // Pure virtual functions that must be implemented
  CVector<T> getMagnetisation() const { return this->mag; }

  void setMagnetisation(const CVector<T> &newMag) {
    if (newMag.length() < 1e-10) {
      throw std::runtime_error("Magnetization vector cannot be zero!");
    }
    this->mag = newMag;
    this->mag.normalize();
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

  /**
   * @brief Set reference layer with enum
   * Can be used to refer to other layers in stack as reference
   * for this layer.
   * @param reference: an enum: FIXED, TOP, BOTTOM, or CUSTOM
   */
  void setReferenceType(Reference reference) {
    if ((reference == FIXED) && (!this->referenceLayer.length())) {
      throw std::runtime_error("Cannot set fixed polarisation layer to 0!"
                               " Set reference to NONE to disable reference.");
    }
    this->referenceType = reference;
  }

  CVector<T> getReferenceLayer() const { return this->referenceLayer; }

  Reference getReferenceType() const { return this->referenceType; }
  std::string getId() { return this->id; }

  // Common driver setters that must be implemented
  void setTemperatureDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->temperatureDriver = driver;
    this->isStochastic = true;
  }
  void setCurrentDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->currentDriver = driver;
  }
  void setAnisotropyDriver(const std::shared_ptr<Driver<T>> &driver) {
    this->anisotropyDriver = driver;
  }
  void setExternalFieldDriver(const std::shared_ptr<AxialDriver<T>> &driver) {
    this->externalFieldDriver = driver;
  }
  void setOerstedFieldDriver(const std::shared_ptr<AxialDriver<T>> &driver) {
    this->HoeDriver = driver;
  }
  void setHdmiDriver(const std::shared_ptr<AxialDriver<T>> &driver) {
    this->HdmiDriver = driver;
  }
  void setReservedInteractionFieldDriver(
      const std::shared_ptr<AxialDriver<T>> &driver) {
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
  typedef void (AbstractLayer<T>::*scalarDriverSetter)(
      const std::shared_ptr<Driver<T>> &);
  typedef void (AbstractLayer<T>::*axialDriverSetter)(
      const std::shared_ptr<AxialDriver<T>> &);

  virtual void rk4_step(const T &time, const T &timeStep,
                        const CVector<T> &bottom, const CVector<T> &top) = 0;

  BufferedNoiseParameters<T> getBufferedNoiseParameters() const {
    return this->noiseParams;
  }

  // void setBufferedNoiseParameters(const BufferedNoiseParameters<T>& params) {
  //     this->noiseParams = params;
  // }

  // void createBufferedAlphaNoise(const unsigned int totalIterations) = 0;
};

template <typename T> class AbstractJunction {
protected:
  // Common member variables
  std::vector<std::shared_ptr<AbstractLayer<T>>>
      layers; // Changed to store pointers to AbstractLayer
  unsigned int logLength = 0;
  unsigned int layerNo = 0;
  std::string Rtag = "R";
  std::unordered_map<std::string, std::vector<T>> log;

  const std::vector<std::string> vectorNames = {"x", "y", "z"};
  T time = 0;

public:
  AbstractJunction(
      const std::vector<std::shared_ptr<AbstractLayer<T>>> &layers) {
    this->layers = layers;
    this->layerNo = layers.size();
    if (this->layerNo == 0) {
      throw std::invalid_argument("Passed a zero length Layer vector!");
    }
    // Verify that all layers have unique ids
    std::unordered_set<std::string> _ids;
    for (const auto &layer : this->layers) {
      if (_ids.find(layer->getId()) != _ids.end()) {
        throw std::invalid_argument("Layers must have unique ids!");
      }
      _ids.insert(layer->getId());
    }
  }

  /**
   * @brief Get Ids of the layers in the junction.
   * @return vector of layer ids.
   */
  const std::vector<std::string> getLayerIds() const {
    std::vector<std::string> ids;
    std::transform(this->layers.begin(), this->layers.end(),
                   std::back_inserter(ids),
                   [](const std::shared_ptr<AbstractLayer<T>> &layer) {
                     return layer->getId();
                   });
    return ids;
  }

  /**
   * @brief Gets a specific layer from the junction
   *
   * @param index Layer index
   * @return Shared pointer to the layer
   */
  std::shared_ptr<AbstractLayer<T>> getLayer(size_t index) const {
    if (index >= layers.size()) {
      throw std::out_of_range("Layer index out of range");
    }
    return layers[index];
  }

  unsigned int getLayerCount() const { return this->layerNo; }

  typedef void (AbstractLayer<T>::*scalarDriverSetter)(
      const std::shared_ptr<Driver<T>> &);
  typedef void (AbstractLayer<T>::*axialDriverSetter)(
      const std::shared_ptr<AxialDriver<T>> &);
  void scalarlayerSetter(const std::string &layerID, scalarDriverSetter functor,
                         std::shared_ptr<Driver<T>> driver) {
    bool found = false;
    for (auto &l : this->layers) {
      if (l->getId() == layerID || layerID == "all") {
        ((*l).*functor)(driver);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  void axiallayerSetter(const std::string &layerID, axialDriverSetter functor,
                        std::shared_ptr<AxialDriver<T>> driver) {
    bool found = false;
    for (auto &l : this->layers) {
      if (l->getId() == layerID || layerID == "all") {
        ((*l).*functor)(driver);
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
  void setCouplingDriver(const std::string &bottomLayer,
                         const std::string &topLayer,
                         const std::shared_ptr<Driver<T>> driver,
                         void (AbstractLayer<T>::*setDriverFuncTop)(
                             const std::shared_ptr<Driver<T>> &),
                         void (AbstractLayer<T>::*setDriverFuncBottom)(
                             const std::shared_ptr<Driver<T>> &)) {
    bool found = false;
    for (unsigned int i = 0; i < this->layerNo - 1; i++) {
      // check if the layer above is actually top layer the user specified
      if ((this->layers[i]->getId() == bottomLayer) &&
          (this->layers[i + 1]->getId() == topLayer)) {
        (this->layers[i]->*setDriverFuncTop)(driver);
        (this->layers[i + 1]->*setDriverFuncBottom)(driver);
        found = true;
        break;
      } else if ((this->layers[i]->getId() == topLayer) &&
                 (this->layers[i + 1]->getId() == bottomLayer)) {
        (this->layers[i]->*setDriverFuncTop)(driver);
        (this->layers[i + 1]->*setDriverFuncBottom)(driver);
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
  void setCouplingDriverAxial(const std::string &bottomLayer,
                              const std::string &topLayer,
                              const std::shared_ptr<AxialDriver<T>> driver,
                              void (AbstractLayer<T>::*setDriverFuncTop)(
                                  const std::shared_ptr<AxialDriver<T>>),
                              void (AbstractLayer<T>::*setDriverFuncBottom)(
                                  const std::shared_ptr<AxialDriver<T>>)) {
    bool found = false;
    for (unsigned int i = 0; i < this->layerNo - 1; i++) {
      // check if the layer above is actually top layer the user specified
      if ((this->layers[i]->getId() == bottomLayer) &&
          (this->layers[i + 1]->getId() == topLayer)) {
        (this->layers[i]->*setDriverFuncTop)(driver);
        (this->layers[i + 1]->*setDriverFuncBottom)(driver);
        found = true;
        break;
      } else if ((this->layers[i]->getId() == topLayer) &&
                 (this->layers[i + 1]->getId() == bottomLayer)) {
        (this->layers[i]->*setDriverFuncTop)(driver);
        (this->layers[i + 1]->*setDriverFuncBottom)(driver);
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

  // Pure virtual functions for magnetoresistance
  virtual std::vector<T> getMagnetoresistance() = 0;
  typedef void (AbstractLayer<T>::*solverFn)(const T &, const T &,
                                             const CVector<T> &,
                                             const CVector<T> &);
  typedef void (AbstractJunction<T>::*runnerFn)(solverFn, T &, T &);
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
  void runSingleLayerSolver(solverFn functor, T &t, T &timeStep) {
    CVector<T> null;
    (this->layers[0].get()->*functor)(t, timeStep, null, null);
  }

  /**
   * @brief Select a solver based on the setup.
   *
   * Multilayer layer solver iteration.
   * @param functor: solver function.
   * @param t: current time
   * @param timeStep: integration step
   * */
  void runMultiLayerSolver(solverFn functor, T &t, T &timeStep) {
    // initialise with 0 CVectors
    std::vector<CVector<T>> magCopies(this->layerNo + 2, CVector<T>());
    // the first and the last layer get 0 vector coupled
    for (unsigned int i = 0; i < this->layerNo; i++) {
      magCopies[i + 1] = this->layers[i]->mag;
    }

    for (unsigned int i = 0; i < this->layerNo; i++) {
      (this->layers[i].get()->*functor)(t, timeStep, magCopies[i],
                                        magCopies[i + 2]);
    }
  }

  std::tuple<runnerFn, solverFn, SolverMode>
  getSolver(SolverMode mode, unsigned int totalIterations) {
    SolverMode localMode = mode;
    for (auto &l : this->layers) {
      if (l->isStochastic) {
        localMode = HEUN;
        break;
      }
    }

    // Define the correct solver function pointer
    solverFn solver = &AbstractLayer<T>::rk4_step;
    runnerFn runner;

    if (this->layerNo == 1) {
      runner = &AbstractJunction<T>::runSingleLayerSolver;
    } else {
      runner = &AbstractJunction<T>::runMultiLayerSolver;
    }

    return std::make_tuple(runner, solver, localMode);
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
                     SolverMode mode = HEUN) {
    if (timeStep > writeFrequency) {
      throw std::runtime_error(
          "The time step cannot be larger than write frequency!");
    }
    const unsigned int totalIterations =
        static_cast<unsigned int>(totalTime / timeStep);
    const unsigned int writeEvery =
        static_cast<unsigned int>(writeFrequency / timeStep);

    // pick a solver based on drivers
    auto [runner, solver, _] = getSolver(mode, totalIterations);

    for (unsigned int i = 0; i < totalIterations; i++) {
      this->time += timeStep;
      (*this.*runner)(solver, this->time, timeStep);

      if (!(i % writeEvery)) {
        logLayerParams(this->time);
      }
    }
    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
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
   * @brief Log junction parameters
   *
   * @param time Current simulation time
   */
  void logLayerParams(const T &time) {
    // Log each layer's parameters
    for (size_t i = 0; i < layers.size(); i++) {
      CVector<T> mag = layers[i]->getMagnetisation();

      // Log magnetization components
      std::string layerPrefix = "layer" + std::to_string(i) + "_";
      this->log[layerPrefix + "mx"].push_back(mag[0]);
      this->log[layerPrefix + "my"].push_back(mag[1]);
      this->log[layerPrefix + "mz"].push_back(mag[2]);
    }

    // Log junction parameters
    this->log["time"].push_back(time);

    this->logLength = this->log["time"].size();
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
};

#endif // CORE_ABSTRACT_HPP
