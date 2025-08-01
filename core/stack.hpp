#ifndef CORE_STACK_HPP_
#define CORE_STACK_HPP_

#include "cvector.hpp" // for CVector
#include "drivers.hpp" // for ScalarDriver, AxialDriver, NullDriver
#include "junction.hpp"
#include <algorithm>     // for for_each
#include <iostream>      // for operator<<, ofstream, string, basi...
#include <numeric>       // for accumulate
#include <stddef.h>      // for size_t
#include <stdexcept>     // for runtime_error
#include <string>        // for operator+, operator==, to_string
#include <unordered_map> // for unordered_map
#include <vector>        // for vector

template <typename T> class Stack {
  friend class Junction<T>;

private:
  ScalarDriver<T> currentDriver;
  std::unordered_map<std::string, std::vector<T>> stackLog;
  bool currentDriverSet = false;

  // Add function pointer for simulation method
  void (Stack<T>::*simulationMethod)(T, T, T);

  // Helper methods for common functionality
  struct SimulationSetup {
    unsigned int writeEvery;
    unsigned int totalIterations;
    std::vector<SolverMode> modes;
    RunnerFn<T> localRunner;
    SolverFn<T> solver;
  };

  SimulationSetup initializeSimulation(T totalTime, T timeStep,
                                       T writeFrequency) {
    SimulationSetup setup;
    setup.writeEvery = (int)(writeFrequency / timeStep);
    setup.totalIterations = (int)(totalTime / timeStep);

    if (timeStep > writeFrequency) {
      throw std::runtime_error(
          "The time step cannot be larger than write frequency!");
    }

    // pick a solver based on drivers
    setup.localRunner = &Junction<T>::runMultiLayerSolver;
    for (auto &j : this->junctionList) {
      auto [runner, solver, mode] = j.getSolver(RK4, setup.totalIterations);
      setup.modes.push_back(mode);
      setup.localRunner = runner;
    }
    setup.solver = &Layer<T>::rk4_step; // legacy, this actually doesn't matter
    if (!std::equal(setup.modes.begin() + 1, setup.modes.end(),
                    setup.modes.begin())) {
      throw std::runtime_error(
          "Junctions have different solver modes!"
          " Set the same solver mode for all junctions explicitly."
          " Do not mix stochastic and deterministic solvers!");
    }

    return setup;
  }

  void storeMagnetisations(std::vector<CVector<T>> &frozenMags,
                           std::vector<CVector<T>> &frozenPols,
                           bool isTwoLayerStack) {
    for (std::size_t j = 0; j < junctionList.size(); ++j) {
      frozenMags[j] = junctionList[j].getLayerMagnetisation(this->topId);
      if (isTwoLayerStack) {
        frozenPols[j] = junctionList[j].getLayerMagnetisation(this->bottomId);
      } else {
        frozenPols[j] = this->getPolarisationVector();
      }
    }
  }

  void executeJunctionStep(std::size_t junctionIndex,
                           ScalarDriver<T> &localDriver, SimulationSetup &setup,
                           T t, T timeStep, std::vector<T> &timeResistances,
                           std::vector<T> &timeCurrents) {
    junctionList[junctionIndex].setLayerCurrentDriver("all", localDriver);
    bool step_accepted = true;
    (junctionList[junctionIndex].*(setup.localRunner))(setup.solver, t,
                                                       timeStep, step_accepted);

    const auto resistance = junctionList[junctionIndex].getMagnetoresistance();
    timeResistances[junctionIndex] = resistance[0];
    timeCurrents[junctionIndex] = localDriver.getCurrentScalarValue(t);
  }

  void logSimulationStep(unsigned int iteration, unsigned int writeEvery, T t,
                         T timeStep, const std::vector<T> &timeResistances,
                         const std::vector<T> &timeCurrents,
                         const std::vector<T> &timeEffectiveCoupling) {
    if (!(iteration % writeEvery)) {
      const T magRes = this->calculateStackResistance(timeResistances);
      this->logStackData(t, magRes, timeCurrents, timeEffectiveCoupling);
      for (auto &jun : this->junctionList)
        jun.logLayerParams(t, timeStep, false);
    }
  }

protected:
  unsigned int stackSize;
  std::string topId, bottomId; // Ids of the top and bottom junctions
  std::vector<T> couplingStrength = {0};
  bool delayed = true;
  T phaseOffset = 0;
  virtual T calculateStackResistance(std::vector<T> resistances) = 0;
  virtual T getPhaseOffset(const unsigned int &order) const = 0;
  virtual T getEffectiveCouplingStrength(const unsigned int &order,
                                         const CVector<T> &m1,
                                         const CVector<T> &m2,
                                         const CVector<T> &p1,
                                         const CVector<T> &p2) = 0;

  T computeCouplingCurrentDensity(const unsigned int &order, T currentDensity,
                                  const CVector<T> &m1, const CVector<T> &m2,
                                  const CVector<T> &p1, const CVector<T> &p2) {
    return currentDensity *
           this->getEffectiveCouplingStrength(order, m1, m2, p1, p2);
  }

public:
  std::vector<Junction<T>> junctionList;

  Stack(std::vector<Junction<T>> inputStack, const std::string &topId,
        const std::string &bottomId, const T phaseOffset = 0,
        bool useKCL = true)
      : topId(topId), bottomId(bottomId), phaseOffset(phaseOffset) {

    // Set the simulation method based on constructor parameter
    if (useKCL) {
      simulationMethod = &Stack<T>::runSimulationKCL;
    } else {
      simulationMethod = &Stack<T>::runSimulationNonKCL;
    }

    if (inputStack.size() < 2) {
      throw std::runtime_error("Stack must have at least 2 junctions!");
    }
    this->junctionList = std::move(inputStack);
    if (std::any_of(this->junctionList.begin(), this->junctionList.end(),
                    [](const Junction<T> &j) {
                      return j.MR_mode != Junction<T>::MRmode::CLASSIC;
                    })) {
      throw std::runtime_error(
          "Junction has a non-classic magnetoresitance mode!"
          " Define the junction with Rp and Rap resistance values.");
    }
    stackSize = this->junctionList.size();
  }

  T getCoupling(const unsigned int &order) const {
    if (this->couplingStrength.empty()) {
      throw std::runtime_error("Coupling strength is not set!");
    }
    if (this->couplingStrength.size() == 1) {
      return this->couplingStrength[0];
    }
    return this->couplingStrength[order];
  }

  void setDelayed(bool delay) {
    if (!delay && !this->isTwoLayerMemberStack()) {
      throw std::runtime_error(
          "Non delayed coupling is only supported for 2 layer stacks!");
    }
    this->delayed = delay;
  }

  void setMagnetisation(unsigned int junctionId, const std::string &layerId,
                        CVector<T> mag) {
    this->junctionList[junctionId].setLayerMagnetisation(layerId, mag);
  }

  const CVector<T> getMagnetisation(unsigned int junctionId,
                                    const std::string &layerId) {
    return this->junctionList[junctionId].getLayerMagnetisation(layerId);
  }

  Junction<T> &getJunction(unsigned int junctionId) {
    return this->junctionList.at(junctionId);
  }

  void setJunctionAnisotropyDriver(unsigned int junctionId,
                                   const std::string &layerId,
                                   const ScalarDriver<T> &k) {
    this->junctionList[junctionId].setLayerAnisotropyDriver(layerId, k);
  }

  void setOerstedFieldDriver(const AxialDriver<T> &oDriver) {
    for (auto &j : this->junctionList) {
      j.setLayerOerstedFieldDriver("all", oDriver);
    }
  }

  void setExternalFieldDriver(const AxialDriver<T> &fDriver) {
    for (auto &j : this->junctionList) {
      j.setLayerExternalFieldDriver("all", fDriver);
    }
  }

  void resetCoupledCurrentDriver() {
    this->currentDriver = NullDriver<T>();
    for (auto &j : this->junctionList) {
      j.setLayerCurrentDriver("all", this->currentDriver);
    }
    this->currentDriverSet = false;
  }

  void setCoupledCurrentDriver(const ScalarDriver<T> &cDriver) {
    this->currentDriver = cDriver;
    for (auto &j : this->junctionList) {
      j.setLayerCurrentDriver("all", this->currentDriver);
    }
    this->currentDriverSet = true;
  }

  void saveLogs(const std::string &fileSave) {
    if (fileSave == "") {
      // if there's an empty fn, don't save
      std::cout << "Ignoring file save to an empty filename" << std::endl;
      return;
    }
    std::ofstream logFile;
    logFile.open(fileSave);
    for (const auto &keyPair : this->stackLog) {
      logFile << keyPair.first << ";";
    }
    logFile << "\n";
    for (unsigned int i = 0; i < this->stackLog["time"].size(); i++) {
      for (const auto &keyPair : this->stackLog) {
        logFile << keyPair.second[i] << ";";
      }
      logFile << "\n";
    }
    logFile.close();
  }

  void setCouplingStrength(const T &coupling) {
    this->couplingStrength = {coupling};
  }

  void setCouplingStrength(const std::vector<T> &coupling) {
    if (coupling.size() != this->stackSize - 1) {
      throw std::runtime_error(
          "Coupling strength vector must have size of stack size - 1!");
    }
    this->couplingStrength = coupling;
  }

  void logStackData(T t, T resistance, std::vector<T> timeCurrents,
                    std::vector<T> effectiveCoupling) {
    this->stackLog["Resistance"].push_back(resistance);
    for (std::size_t j = 0; j < timeCurrents.size(); ++j) {
      this->stackLog["I_" + std::to_string(j)].push_back(timeCurrents[j]);
      if (j < effectiveCoupling.size()) {
        this->stackLog["C_" + std::to_string(j) + "_" + std::to_string(j + 1)]
            .push_back(effectiveCoupling[j]);
      }
    }
    this->stackLog["time"].push_back(t);
  }

  void clearLogs() {
    for (auto &j : this->junctionList) {
      j.clearLog();
    }
    this->stackLog.clear();
  }

  std::unordered_map<std::string, std::vector<T>> &getLog() {
    return this->stackLog;
  }
  std::unordered_map<std::string, std::vector<T>> &getLog(unsigned int id) {
    if (id < this->junctionList.size()) {
      return this->junctionList[id].getLog();
    }
    throw std::runtime_error("Asking for id of a non-existing junction!");
  }

  const CVector<T> getPolarisationVector() {
    CVector<T> probe = junctionList[0].getLayer(this->topId).referenceLayer;
    for (std::size_t i = 1; i < junctionList.size(); ++i) {
      if (probe != junctionList[i].getLayer(this->topId).referenceLayer)
        throw std::runtime_error("Polarisation vectors are not equal in stack");
    }

    if (!probe.length()) {
      throw std::runtime_error("Polarisation is not set!");
    }
    return probe;
  }

  const bool isTwoLayerMemberStack() {
    for (const auto &j : this->junctionList) {
      if (j.layerNo >= 3) {
        throw std::runtime_error("At least one junction has more than 2 layers!"
                                 " It's not supported now.");
      } else if (j.layerNo != 2) {
        return false;
      }
    }
    return true;
  }

  // Public interface method that delegates to the chosen implementation
  void runSimulation(T totalTime, T timeStep = 1e-13,
                     T writeFrequency = 1e-11) {
    (this->*simulationMethod)(totalTime, timeStep, writeFrequency);
  }

private:
  void runSimulationNonKCL(T totalTime, T timeStep = 1e-13,
                           T writeFrequency = 1e-11) {

    auto setup = initializeSimulation(totalTime, timeStep, writeFrequency);

    std::vector<T> timeResistances(junctionList.size());
    std::vector<T> timeCurrents(junctionList.size());
    std::vector<T> timeEffectiveCoupling(junctionList.size() - 1);
    std::vector<CVector<T>> frozenMags(junctionList.size());
    std::vector<CVector<T>> frozenPols(junctionList.size());

    const bool isTwoLayerStack = this->isTwoLayerMemberStack();
    for (unsigned int i = 0; i < setup.totalIterations; i++) {
      T t = i * timeStep;

      // Store magnetisations first
      storeMagnetisations(frozenMags, frozenPols, isTwoLayerStack);

      T effectiveCoupling = 1;
      for (std::size_t j = 0; j < junctionList.size(); ++j) {

        /**
         * Coupling
         * Ik = Ik-1 + x* Ik-1 = (1+x)Ik-1
         * Ik+1 = Ik + x* Ik = (1+x)Ik = (1+x)(1+x)Ik-1
         * technically we could do (1+x)^n * I0 but
         * we want to expand to non-symmetric coupling x1, x2, ...
         */

        // skip first junction
        // modify the standing layer constant current
        if (j > 0) {
          if (this->delayed) {
            // accumulate coupling
            effectiveCoupling *=
                (1 + this->getEffectiveCouplingStrength(
                         j - 1, frozenMags[j - 1], frozenMags[j],
                         frozenPols[j - 1], frozenPols[j]));

          } else {
            effectiveCoupling *=
                (1 +
                 this->getEffectiveCouplingStrength(
                     j - 1,
                     junctionList[j - 1].getLayerMagnetisation(this->topId),
                     junctionList[j].getLayerMagnetisation(this->topId),
                     junctionList[j - 1].getLayerMagnetisation(this->bottomId),
                     junctionList[j].getLayerMagnetisation(this->bottomId)));
          }
          timeEffectiveCoupling[j - 1] = effectiveCoupling;
        }
        // set the current -- same for all layers
        // copy the driver and set the current value
        ScalarDriver<T> localDriver = this->currentDriver * effectiveCoupling;
        localDriver.phaseShift(this->getPhaseOffset(j));

        executeJunctionStep(j, localDriver, setup, t, timeStep, timeResistances,
                            timeCurrents);
      }
      logSimulationStep(i, setup.writeEvery, t, timeStep, timeResistances,
                        timeCurrents, timeEffectiveCoupling);
    }
  }

  void runSimulationKCL(T totalTime, T timeStep = 1e-13,
                        T writeFrequency = 1e-11) {

    auto setup = initializeSimulation(totalTime, timeStep, writeFrequency);

    std::vector<T> timeResistances(junctionList.size());
    std::vector<T> timeCurrents(junctionList.size());
    std::vector<T> timeEffectiveCoupling(junctionList.size() - 1);
    std::vector<CVector<T>> frozenMags(junctionList.size());
    std::vector<CVector<T>> frozenPols(junctionList.size());
    const bool isTwoLayerStack = this->isTwoLayerMemberStack();

    for (unsigned int i = 0; i < setup.totalIterations; i++) {
      T t = i * timeStep;

      // this is a base case
      T uncoupledCurrent = this->currentDriver.getCurrentScalarValue(t);

      // Store magnetisations first
      storeMagnetisations(frozenMags, frozenPols, isTwoLayerStack);

      T totalCurrent = uncoupledCurrent;
      for (std::size_t j = 0; j < junctionList.size() - 1; ++j) {
        timeEffectiveCoupling[j] = this->getEffectiveCouplingStrength(
            j, frozenMags[j], frozenMags[j + 1], frozenPols[j],
            frozenPols[j + 1]);
        totalCurrent += timeEffectiveCoupling[j] * uncoupledCurrent;
      }

      for (std::size_t j = 0; j < junctionList.size(); ++j) {
        // set the current -- same for all layers
        // copy the driver and set the current value
        ScalarDriver<T> localDriver =
            ScalarDriver<T>::getConstantDriver(totalCurrent);
        localDriver.phaseShift(this->getPhaseOffset(j));

        executeJunctionStep(j, localDriver, setup, t, timeStep, timeResistances,
                            timeCurrents);
      }
      logSimulationStep(i, setup.writeEvery, t, timeStep, timeResistances,
                        timeCurrents, timeEffectiveCoupling);
    }
  }
};
template <typename T> class SeriesStack : public Stack<T> {
  T calculateStackResistance(std::vector<T> resistances) override {
    const T resSum =
        std::accumulate(resistances.begin(), resistances.end(), 0.0);
    return resSum;
  }

  T getEffectiveCouplingStrength(const unsigned int &order,
                                 const CVector<T> &m1, const CVector<T> &m2,
                                 const CVector<T> &p1,
                                 const CVector<T> &p2) override {
    const T m1Comp = c_dot(m1, p1);
    const T m2Comp = c_dot(m2, p2);
    return this->getCoupling(order) * (m1Comp + m2Comp);
  }

  T getPhaseOffset(const unsigned int &order) const override {
    return this->phaseOffset * order;
  }

public:
  explicit SeriesStack(const std::vector<Junction<T>> &jL,
                       const std::string &topId = "free",
                       const std::string &bottomId = "bottom",
                       const T phaseOffset = 0, bool useKCL = true)
      : Stack<T>(jL, topId, bottomId, phaseOffset, useKCL) {}
};

template <typename T> class ParallelStack : public Stack<T> {
  T calculateStackResistance(std::vector<T> resistances) override {
    T invSum = 0.0;
    std::for_each(resistances.begin(), resistances.end(),
                  [&](T res) { invSum += 1.0 / res; });
    return 1. / invSum;
  }

  T getEffectiveCouplingStrength(const unsigned int &order,
                                 const CVector<T> &m1, const CVector<T> &m2,
                                 const CVector<T> &p1,
                                 const CVector<T> &p2) override {
    const T m1Comp = c_dot(m1, p1);
    const T m2Comp = c_dot(m2, p2);
    return this->getCoupling(order) * (m1Comp - m2Comp);
  }

  T getPhaseOffset(const unsigned int &order) const override {
    return this->phaseOffset;
  }

public:
  explicit ParallelStack(const std::vector<Junction<T>> &jL,
                         const std::string &topId = "free",
                         const std::string &bottomId = "bottom",
                         const T phaseOffset = 0, bool useKCL = true)
      : Stack<T>(jL, topId, bottomId, phaseOffset, useKCL) {}
};

#endif // CORE_STACK_HPP_
