#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "cvector.hpp"
#include "drivers.hpp"
#include "junction.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Computes the combinations
 * https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
 * @param N size of the set
 * @param K combination size
 */
void comb(int N, int K) {
  std::string bitmask(K, 1); // K leading 1's
  bitmask.resize(N, 0);      // N-K trailing 0's
  // print integers and permute bitmask
  do {
    for (int i = 0; i < N; ++i) // [0..N-1] integers
    {
      if (bitmask[i])
        std::cout << " " << i;
    }
    std::cout << std::endl;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

typedef std::array<CVector<double>, 3> tensor;
typedef std::vector<tensor> tensorList;

const tensor getDipoleTensorFromRelPositions(const CVector<double> &r1,
                                             const CVector<double> &r2) {
  const CVector<double> rij = r1 - r2; // 1-2 distance vector
  if (rij.length() < 1e-14) {
    throw std::runtime_error(
        "Points are too close for stable dipole calculation");
  }
  const double r_mag = pow(rij.length(), 2);
  const double mult = 3 / (4 * M_PI * pow(rij.length(), 5));
  const tensor dipoleTensor = {CVector<double>(pow(rij.x, 2) - (r_mag / 3),
                                               rij.x * rij.y, rij.x * rij.z) *
                                   mult,
                               CVector<double>(rij.x * rij.y,
                                               pow(rij.y, 2) - (r_mag / 3),
                                               rij.y * rij.z) *
                                   mult,
                               CVector<double>(rij.x * rij.z, rij.y * rij.z,
                                               pow(rij.z, 2) - (r_mag / 3)) *
                                   mult};
  return dipoleTensor;
}

typedef std::function<CVector<double>(
    const CVector<double> &, const CVector<double> &, const Layer<double> &,
    const Layer<double> &)>
    interactionFunction;

CVector<double> nullDipoleInteraction(const CVector<double> &,
                                      const CVector<double> &,
                                      const Layer<double> &,
                                      const Layer<double> &) {
  return CVector<double>(0, 0, 0);
}

/**
 * @brief Compute dipole interaction between two junctions.
 * From: Kanao et al, Reservoir Computing on Spin-Torque Oscillator Array (2019)
 * PRA
 * @param r1 1st junction position
 * @param r2 2nd junction position
 * @param layer1 1st junction layer
 * @param layer2 2nd junction layer
 */
CVector<double> computeDipoleInteraction(const CVector<double> &r1,
                                         const CVector<double> &r2,
                                         const Layer<double> &layer1,
                                         const Layer<double> &layer2) {
  const DVector rij = r1 - r2; // 1-2 distance vector
  if (rij.length() < 1e-10) {
    throw std::runtime_error(
        "Points are too close for stable dipole calculation");
  }
  const double r3 = pow(rij.length(), 3);
  const double r5 = pow(rij.length(), 5);
  const Layer<double> ref_magnetic_moment = layer2;
  const DVector m1 = ref_magnetic_moment.mag;
  const double V =
      ref_magnetic_moment.thickness * ref_magnetic_moment.cellSurface;
  const double prefactor = (ref_magnetic_moment.Ms / MAGNETIC_PERMEABILITY) * V;

  return prefactor * (3 * c_dot(m1, rij) * rij / r5 - m1 / r3);
}

/**
 * @brief Compute dipole interaction between two junctions.
 * From: Nomura et al, Reservoir computing with dipole-coupled nanomagnets
 * (2019) JJAP
 * @param r1 1st junction position
 * @param r2 2nd junction position
 * @param layer1 1st junction layer
 * @param layer2 2nd junction layer
 */
CVector<double> computeDipoleInteractionNoumra(const CVector<double> &r1,
                                               const CVector<double> &r2,
                                               const Layer<double> &layer1,
                                               const Layer<double> &layer2) {
  const tensor dipoleTensor = getDipoleTensorFromRelPositions(r1, r2);
  const double V = layer2.thickness * layer2.cellSurface;
  // remember that calculate_tensor_interaction normalises Ms by 1/mu0
  const CVector<double> dipoleVector =
      calculate_tensor_interaction(layer2.mag, dipoleTensor, layer2.Ms * V);
  // in the paper they don't multiply explicitly by V, but it's necessary for
  // units to match
  return dipoleVector;
}

class GroupInteraction {
  std::string topId; // Id of the top junction
  std::vector<CVector<double>> coordinateMatrix;
  std::vector<Junction<double>> junctionList;
  interactionFunction interactionFunc = computeDipoleInteraction;
  unsigned int noElements;

  void stepFunctionalSolver(double time, double &timeStep,
                            interactionFunction interaction,
                            RunnerFn<double> runner, SolverFn<double> solver,
                            bool &step_accepted) {
    // collect all frozen states
    // for each element, compute the extra field from all other elements
    for (unsigned int i = 0; i < this->noElements; i++) {
      CVector<double> H_extra;
      for (unsigned int j = 0; j < this->noElements; j++) {
        if (i == j)
          continue;
        H_extra +=
            interaction(this->coordinateMatrix[i], this->coordinateMatrix[j],
                        this->junctionList[i].getLayer(this->topId),
                        this->junctionList[j].getLayer(this->topId));
      }
      this->junctionList[i].setLayerReservedInteractionField(
          this->topId, AxialDriver<double>(H_extra));
    }
    // step the solver with the extra field
    for (unsigned int i = 0; i < this->noElements; i++) {
      (this->junctionList[i].*runner)(solver, time, timeStep, step_accepted);
    }
  }

public:
  GroupInteraction(std::vector<CVector<double>> coordinateMatrix,
                   std::vector<Junction<double>> junctionList,
                   const std::string &topId = "free")
      : topId(topId) {
    if (coordinateMatrix.size() != junctionList.size()) {
      throw std::runtime_error(
          "Coordinate matrix and junction list must have the same size!");
    }
    if (coordinateMatrix.empty() || junctionList.empty()) {
      throw std::runtime_error(
          "Coordinate matrix and junction list cannot be empty!");
    }
    this->noElements = junctionList.size();
    this->coordinateMatrix = std::move(coordinateMatrix);
    this->junctionList = std::move(junctionList);

    // check that all vectors in coordinateMatrix are unique
    for (size_t i = 0; i < this->coordinateMatrix.size(); i++) {
      for (size_t j = i + 1; j < this->coordinateMatrix.size(); j++) {
        if (this->coordinateMatrix[i] == this->coordinateMatrix[j]) {
          throw std::runtime_error("Coordinate vectors must be unique!");
        }
      }
    }
  }

  void setInteractionFunction(interactionFunction interactionFunc) {
    this->interactionFunc = interactionFunc;
  }

  void runSimulation(double totalTime, double timeStep = 1e-13,
                     double writeFrequency = 1e-13) {
    const unsigned int writeEvery = (int)(writeFrequency / timeStep);
    const unsigned int totalIterations = (int)(totalTime / timeStep);

    if (timeStep > writeFrequency) {
      throw std::runtime_error(
          "The time step cannot be larger than write frequency!");
    }

    // pick a solver based on drivers
    std::vector<SolverMode> modes;
    auto localRunner = &Junction<double>::runMultiLayerSolver;
    for (auto &j : this->junctionList) {
      // again, solver mode does not make a difference for legacy reasons
      auto [runner, solver, mode] = j.getSolver(RK4, totalIterations);
      modes.push_back(mode);
      localRunner = runner;
      // TODO: handle the rare case when the user mixes 1 layer with 2 layer
      // junction in the same stack -- i.e. runner is runSingleLayerSolver and
      // runMultiLayerSolver
    }
    auto solver = &Layer<double>::rk4_step;
    if (!std::equal(modes.begin() + 1, modes.end(), modes.begin())) {
      throw std::runtime_error(
          "Junctions have different solver modes!"
          " Set the same solver mode for all junctions explicitly."
          " Do not mix stochastic and deterministic solvers!");
    }

    bool step_accepted = true;
    for (unsigned int i = 0; i < totalIterations; i++) {
      double t = i * timeStep;
      stepFunctionalSolver(t, timeStep, this->interactionFunc, localRunner,
                           solver, step_accepted);

      if (!(i % writeEvery)) {
        for (auto &jun : this->junctionList)
          jun.logLayerParams(t, timeStep, false);
      }
    }
  }

  void clearLogs() {
    for (auto &j : this->junctionList) {
      j.clearLog();
    }
  }

  std::unordered_map<std::string, std::vector<double>> &
  getLog(unsigned int id) {
    if (id < this->junctionList.size()) {
      return this->junctionList[id].getLog();
    }
    throw std::runtime_error("Asking for id of a non-existing junction!");
  }

  std::unordered_map<std::string, std::vector<double>> &getLog() {
    throw std::runtime_error("Not implemented!");
  }
};

class Reservoir {
private:
  // log stuff
  std::unordered_map<std::string, std::vector<double>> reservoirLog;

  // reservoir matrices
  std::vector<std::vector<DVector>> coordinateMatrix;
  std::vector<DVector> frozenMMatrix;
  std::vector<double> MsMatrix, volumeMatrix;
  std::vector<std::vector<tensor>> reservoirDipoleTensor;
  std::vector<std::vector<Layer<double>>> layerMatrix;

  std::vector<std::vector<tensor>> computeReservoirDipoleMatrix(
      const std::vector<std::vector<CVector<double>>> &coordinateMatrix) {
    std::vector<std::vector<tensor>> localReservoirDipoleTensor;
    // reserve some place here
    localReservoirDipoleTensor.resize(this->noElements);

    // given a coordinate matrix, create a flat index of all indexed
    std::string bitmask(2, 1);           // K leading 1's
    bitmask.resize(this->noElements, 0); // N-K trailing 0's
    // print integers and permute bitmask
    do {
      std::array<int, 2> consideredPair;
      int asgn = 0;
      for (unsigned int i = 0; i < this->noElements; ++i) // [0..N-1] integers
      {
        if (bitmask[i])
          consideredPair[asgn++] =
              i; // currently selected index in the combination
      }
      // compute matrix position from index
      // this is row-first (row-major) ordering!
      const auto elIndx0 = getMatrixCoordinates(consideredPair[0]);
      const auto elIndx1 = getMatrixCoordinates(consideredPair[1]);
      // const unsigned int i0 = (int)consideredPair[0] / cols; // first
      // position of the first element -- row const unsigned int i1 =
      // consideredPair[0] % cols;      // second position of the first element
      // -- col
      const unsigned int i0 = std::get<0>(elIndx0);
      const unsigned int i1 = std::get<1>(elIndx0);
      const unsigned int j0 = std::get<0>(elIndx1);
      const unsigned int j1 = std::get<1>(elIndx1);
      std::cout << "i0: " << i0 << " i1: " << i1 << " j0: " << j0
                << " j1: " << j1 << std::endl;
      // // const unsigned int j0 = (int)consideredPair[1] / cols; // first
      // position of the second element
      // // const unsigned int j1 = consideredPair[1] % cols;      // second
      // position of the second element

      tensor dipTensor1 = this->getDipoleTensorFromRelPositions(
          coordinateMatrix[i0][i1], coordinateMatrix[j0][j1]);
      tensor dipTensor2 = this->getDipoleTensorFromRelPositions(
          coordinateMatrix[j0][j1], coordinateMatrix[i0][i1]);
      // dipole tensors should be symmetric (or anti-symmetric)
      // TODO: make sure this is the case
      // localReservoirDipoleTensor[i0][i1].push_back(dipTensor1);
      // localReservoirDipoleTensor[j0][j1].push_back(dipTensor2);
      localReservoirDipoleTensor[consideredPair[0]].push_back(dipTensor1);
      localReservoirDipoleTensor[consideredPair[1]].push_back(dipTensor2);

    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return localReservoirDipoleTensor;
  }

  std::tuple<int, int> getMatrixCoordinates(int elementIndx) {
    // this is row-major convention
    return std::make_tuple((int)(elementIndx / this->cols),
                           elementIndx % this->cols);
  }

  const tensor getDipoleTensorFromRelPositions(const CVector<double> &r1,
                                               const CVector<double> &r2) {
    const CVector<double> rij = r2 - r1; // 1-2 distance vector
    const double r_mag = pow(rij.length(), 2);
    const double mult = 3 / (4 * M_PI * pow(rij.length(), 5));
    const tensor dipoleTensor = {CVector<double>(pow(rij.x, 2) - (r_mag / 3),
                                                 rij.x * rij.y, rij.x * rij.z) *
                                     mult,
                                 CVector<double>(rij.x * rij.y,
                                                 pow(rij.y, 2) - (r_mag / 3),
                                                 rij.y * rij.z) *
                                     mult,
                                 CVector<double>(rij.x * rij.z, rij.y * rij.z,
                                                 pow(rij.z, 2) - (r_mag / 3)) *
                                     mult};
    return dipoleTensor;
  }

  CVector<double> computeDipoleInteraction(int currentIndx,
                                           double volumeNormaliser) {
    CVector<double> HdipoleEff;
    for (unsigned int i = 0;
         i < this->reservoirDipoleTensor[currentIndx].size(); i++) {
      HdipoleEff += calculate_tensor_interaction(
          this->frozenMMatrix[i], this->reservoirDipoleTensor[currentIndx][i],
          this->MsMatrix[i]);
    }
    return HdipoleEff * volumeNormaliser;
  }

public:
  unsigned int rows, cols;
  unsigned int noElements;

  Reservoir(std::vector<std::vector<DVector>> coordinateMatrix,
            std::vector<std::vector<Layer<double>>> layerMatrix)
      : coordinateMatrix(std::move(coordinateMatrix)),
        layerMatrix(std::move(layerMatrix)) {
    this->rows = this->coordinateMatrix.size();
    this->cols = this->coordinateMatrix[0].size();
    this->noElements = this->rows * this->cols;
    this->frozenMMatrix.resize(this->noElements);
    this->MsMatrix.reserve(this->noElements);
    this->volumeMatrix.reserve(this->noElements);
    for (unsigned int i = 0; i < this->rows; i++) {
      for (unsigned j = 0; j < this->cols; j++) {
        // must be multiplied by volume to avoid overflow
        this->volumeMatrix.push_back(this->layerMatrix[i][j].thickness *
                                     this->layerMatrix[i][j].cellSurface);
        this->MsMatrix.push_back(this->layerMatrix[i][j].Ms);
      }
    }
    this->reservoirDipoleTensor =
        this->computeReservoirDipoleMatrix(this->coordinateMatrix);
  }

  std::vector<CVector<double>> collectFrozenMMatrix() {
    for (unsigned int i = 0; i < this->noElements; i++) {
      const auto coords = this->getMatrixCoordinates(i);
      const unsigned int i0 = std::get<0>(coords);
      const unsigned int i1 = std::get<1>(coords);
      this->frozenMMatrix[i] = this->layerMatrix[i0][i1].mag;
    }
    return this->frozenMMatrix;
  }

  void runSolver(double time, double timeStep, bool parallel = false) {
    // collect all frozen states
    collectFrozenMMatrix();
    CVector<double> nullVec;
    for (unsigned int i = 0; i < this->noElements; i++) {
      const auto dipoleVector =
          computeDipoleInteraction(i, this->volumeMatrix[i]);
      const auto coords = this->getMatrixCoordinates(i);
      const unsigned int i0 = std::get<0>(coords);
      const unsigned int i1 = std::get<1>(coords);
      // std::cout << dipoleVector.x << " " << dipoleVector.y << " " <<
      // dipoleVector.z << std::endl;
      layerMatrix[i0][i1].rk4_stepDipoleInjection(time, timeStep, nullVec,
                                                  nullVec, dipoleVector);
    }
  }

  void logReservoirkData(double t) {
    this->reservoirLog["time"].push_back(t);
    for (unsigned int i = 0; i < this->noElements; i++) {
      const auto coords = this->getMatrixCoordinates(i);
      const unsigned int i0 = std::get<0>(coords);
      const unsigned int i1 = std::get<1>(coords);
      this->reservoirLog["m_" + std::to_string(i0) + "_" + std::to_string(i1) +
                         "_x"]
          .push_back(this->layerMatrix[i0][i1].mag.x);
      this->reservoirLog["m_" + std::to_string(i0) + "_" + std::to_string(i1) +
                         "_y"]
          .push_back(this->layerMatrix[i0][i1].mag.y);
      this->reservoirLog["m_" + std::to_string(i0) + "_" + std::to_string(i1) +
                         "_z"]
          .push_back(this->layerMatrix[i0][i1].mag.z);
    }
  }

  Layer<double> &getLayer(unsigned int index) {
    const auto coords = getMatrixCoordinates(index);
    const unsigned int i0 = std::get<0>(coords);
    const unsigned int i1 = std::get<1>(coords);
    return this->layerMatrix[i0][i1];
  }

  void setAllExternalField(const AxialDriver<double> &hdriver) {
    for (auto &r : this->layerMatrix) {
      for (auto &l : r) {
        l.setExternalFieldDriver(hdriver);
      }
    }
  }

  void setLayerExternalField(unsigned int index,
                             const AxialDriver<double> &hDriver) {
    this->getLayer(index).setExternalFieldDriver(hDriver);
  }

  void setLayerAnisotropy(unsigned int index,
                          const ScalarDriver<double> &anisotropyDriver) {
    this->getLayer(index).setAnisotropyDriver(anisotropyDriver);
  }

  CVector<double> getMagnetisation(unsigned int index) {
    const auto coords = getMatrixCoordinates(index);
    const unsigned int i0 = std::get<0>(coords);
    const unsigned int i1 = std::get<1>(coords);
    return this->layerMatrix[i0][i1].mag;
  }

  void saveLogs(const std::string &fileSave) {
    if (fileSave == "") {
      // if there's an empty fn, don't save
      std::cout << "Ignoring file save to an empty filename" << std::endl;
      return;
    }
    std::ofstream logFile;
    logFile.open(fileSave);
    for (const auto &keyPair : this->reservoirLog) {
      logFile << keyPair.first << ";";
    }
    logFile << "\n";
    for (unsigned int i = 0; i < this->reservoirLog["time"].size(); i++) {
      for (const auto &keyPair : this->reservoirLog) {
        logFile << keyPair.second[i] << ";";
      }
      logFile << "\n";
    }
    logFile.close();
  }
  void clearLogs() { this->reservoirLog.clear(); }

  void runSimulation(double totalTime, double timeStep) {
    const double totalIterations = (int)(totalTime / timeStep);
    // this->clearLogs();
    // this->prepareLog(totalIterations);
    for (unsigned int i = 0; i < totalIterations; i++) {
      double t = i * timeStep;
      runSolver(t, timeStep);
      logReservoirkData(t);
    }
  }
};

#endif
