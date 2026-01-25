#ifndef CORE_FDM_HPP_
#define CORE_FDM_HPP_

#include "drivers.hpp"
#include "junction.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

template <typename T> struct FDMGridSpec {
  T width = 0.0;
  T length = 0.0;
  T thickness = 0.0;
  T cellSizeXY = 0.0;
  T cellSizeZ = 0.0;
  unsigned int nx = 0;
  unsigned int ny = 0;
  unsigned int nz = 1;
  T dx = 0.0;
  T dy = 0.0;
  T dz = 0.0;

  FDMGridSpec() = default;
  FDMGridSpec(T width, T length, T thickness, T cellSizeXY, unsigned int nz = 1,
              T cellSizeZ = 0.0)
      : width(width), length(length), thickness(thickness),
        cellSizeXY(cellSizeXY), cellSizeZ(cellSizeZ), nz(nz) {
    if (width <= 0 || length <= 0 || thickness <= 0) {
      throw std::runtime_error("Grid dimensions must be positive");
    }
    if (cellSizeXY <= 0) {
      throw std::runtime_error("cellSizeXY must be positive");
    }
    if (nz < 1) {
      throw std::runtime_error("nz must be at least 1");
    }

    nx = static_cast<unsigned int>(std::ceil(width / cellSizeXY));
    ny = static_cast<unsigned int>(std::ceil(length / cellSizeXY));
    dx = width / static_cast<T>(nx);
    dy = length / static_cast<T>(ny);

    if (cellSizeZ <= 0.0) {
      dz = thickness / static_cast<T>(nz);
      cellSizeZ = dz;
    } else {
      dz = cellSizeZ;
      const T expectedThickness = dz * static_cast<T>(nz);
      const T tolerance = static_cast<T>(1e-12);
      if (std::abs(expectedThickness - thickness) > tolerance) {
        throw std::runtime_error("cellSizeZ * nz must match layer thickness");
      }
    }
  }
};

template <typename T> struct FDMGrid {
  unsigned int nx = 0;
  unsigned int ny = 0;
  unsigned int nz = 1;
  T dx = 0.0;
  T dy = 0.0;
  T dz = 0.0;

  FDMGrid() = default;
  explicit FDMGrid(const FDMGridSpec<T> &spec)
      : nx(spec.nx), ny(spec.ny), nz(spec.nz), dx(spec.dx), dy(spec.dy),
        dz(spec.dz) {}

  std::size_t cellCount() const {
    return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
           static_cast<std::size_t>(nz);
  }

  std::size_t index(unsigned int x, unsigned int y, unsigned int z) const {
    return (static_cast<std::size_t>(z) * ny + y) * nx + x;
  }

  std::array<unsigned int, 3> coordinates(std::size_t idx) const {
    const unsigned int x = static_cast<unsigned int>(idx % nx);
    const unsigned int y = static_cast<unsigned int>((idx / nx) % ny);
    const unsigned int z = static_cast<unsigned int>(idx / (nx * ny));
    return {x, y, z};
  }

  T cellVolume() const { return dx * dy * dz; }
  T cellSurface() const { return dx * dy; }
};

template <typename T> class FDMLayer {
public:
  using DemagTensorCell = std::array<CVector<T>, 3>;

private:
  Layer<T> baseLayer;
  FDMGrid<T> grid;
  std::vector<CVector<T>> magnetisation;
  std::vector<DemagTensorCell> demagTensor;
  T exchangeStiffness; // Exchange stiffness A [J/m]

public:
  FDMLayer() = default;
  FDMLayer(const Layer<T> &layerTemplate, const FDMGridSpec<T> &spec)
      : baseLayer(), grid(spec), exchangeStiffness(0) {
    baseLayer = layerTemplate;
    baseLayer.thickness = grid.dz;
    baseLayer.cellSurface = grid.cellSurface();
    baseLayer.cellVolume = grid.cellVolume();
    magnetisation.assign(grid.cellCount(), layerTemplate.mag);
    for (auto &m : magnetisation) {
      m.normalize();
    }
  }

  const std::string &getId() const { return baseLayer.id; }
  const FDMGrid<T> &getGrid() const { return grid; }
  const std::vector<CVector<T>> &getMagnetisationGrid() const {
    return magnetisation;
  }
  Layer<T> &getBaseLayer() { return baseLayer; }
  const Layer<T> &getBaseLayer() const { return baseLayer; }

  void setMagnetisation(const CVector<T> &mag) {
    CVector<T> norm = mag;
    norm.normalize();
    for (auto &m : magnetisation) {
      m = norm;
    }
  }

  void setMagnetisationGrid(const std::vector<CVector<T>> &mags) {
    if (mags.size() != magnetisation.size()) {
      throw std::runtime_error("Magnetisation grid size mismatch");
    }
    magnetisation = mags;
    for (auto &m : magnetisation) {
      m.normalize();
    }
  }

  void setDemagTensor(const std::vector<DemagTensorCell> &tensor) {
    if (!tensor.empty() && tensor.size() != grid.cellCount()) {
      throw std::runtime_error("Demag tensor size mismatch");
    }
    demagTensor = tensor;
  }

  const std::vector<DemagTensorCell> &getDemagTensor() const {
    return demagTensor;
  }

  bool hasDemagTensor() const { return !demagTensor.empty(); }

  void setExchangeStiffness(T A) { exchangeStiffness = A; }
  T getExchangeStiffness() const { return exchangeStiffness; }
};

template <typename T> class FDMJunction {
private:
  std::vector<FDMLayer<T>> layers;
  unsigned int layerNo = 0;

  using ScalarSetter = void (Layer<T>::*)(const ScalarDriver<T> &);
  using AxialSetter = void (Layer<T>::*)(const AxialDriver<T> &);

  void scalarlayerSetter(const std::string &layerID, ScalarSetter functor,
                         ScalarDriver<T> driver) {
    bool found = false;
    for (auto &layer : layers) {
      if (layer.getId() == layerID || layerID == "all") {
        (layer.getBaseLayer().*functor)(driver);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  void axiallayerSetter(const std::string &layerID, AxialSetter functor,
                        AxialDriver<T> driver) {
    bool found = false;
    for (auto &layer : layers) {
      if (layer.getId() == layerID || layerID == "all") {
        (layer.getBaseLayer().*functor)(driver);
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "Failed to find a layer with a given id: " + layerID + "!");
    }
  }

  void setCouplingDriver(
      const std::string &bottomLayer, const std::string &topLayer,
      const ScalarDriver<T> &driver,
      void (Layer<T>::*setDriverFuncTop)(const ScalarDriver<T> &),
      void (Layer<T>::*setDriverFuncBottom)(const ScalarDriver<T> &)) {
    bool found = false;
    for (unsigned int i = 0; i < layerNo - 1; i++) {
      const auto &current = layers[i].getId();
      const auto &next = layers[i + 1].getId();
      if ((current == bottomLayer && next == topLayer) ||
          (current == topLayer && next == bottomLayer)) {
        (layers[i].getBaseLayer().*setDriverFuncTop)(driver);
        (layers[i + 1].getBaseLayer().*setDriverFuncBottom)(driver);
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

  SolverMode resolveSolverMode(SolverMode mode, unsigned int totalIterations) {
    SolverMode localMode = mode;
    for (auto &layer : layers) {
      auto &baseLayer = layer.getBaseLayer();
      if (baseLayer.hasTemperature()) {
        if (localMode != HEUN && localMode != EULER_HEUN) {
          std::cout << "[WARNING] Solver automatically changed to Euler Heun "
                       "for stochastic calculation."
                    << std::endl;
          localMode = EULER_HEUN;
        }
      }
      if (baseLayer.noiseParams.scaleNoise != 0) {
        if (localMode != HEUN && localMode != EULER_HEUN) {
          std::cout << "[WARNING] Solver automatically changed to Euler Heun "
                       "for stochastic calculation."
                    << std::endl;
          localMode = EULER_HEUN;
        }
        baseLayer.createBufferedAlphaNoise(totalIterations);
      }
    }
    return localMode;
  }

  CVector<T> getCoupledMagnetisation(const std::vector<CVector<T>> &mags,
                                     const FDMGrid<T> &grid, unsigned int x,
                                     unsigned int y, unsigned int z) const {
    if (grid.cellCount() == 0) {
      return CVector<T>();
    }
    const unsigned int zIndex = std::min(z, grid.nz - 1);
    return mags[grid.index(x, y, zIndex)];
  }

  std::vector<CVector<T>>
  computeDemagFields(const FDMLayer<T> &layer,
                     const std::vector<CVector<T>> &mags) const {
    std::vector<CVector<T>> fields(mags.size(), CVector<T>());
    if (!layer.hasDemagTensor()) {
      return fields;
    }
    const auto &tensor = layer.getDemagTensor();
    const auto &grid = layer.getGrid();
    for (std::size_t idx = 0; idx < mags.size(); ++idx) {
      const auto coords = grid.coordinates(idx);
      CVector<T> field(0, 0, 0);
      for (std::size_t j = 0; j < mags.size(); ++j) {
        const auto jcoords = grid.coordinates(j);
        const unsigned int dx = (coords[0] > jcoords[0])
                                    ? (coords[0] - jcoords[0])
                                    : (jcoords[0] - coords[0]);
        const unsigned int dy = (coords[1] > jcoords[1])
                                    ? (coords[1] - jcoords[1])
                                    : (jcoords[1] - coords[1]);
        const unsigned int dz = (coords[2] > jcoords[2])
                                    ? (coords[2] - jcoords[2])
                                    : (jcoords[2] - coords[2]);
        const auto tensorIndex = grid.index(dx, dy, dz);
        field =
            field + calculate_tensor_interaction(mags[j], tensor[tensorIndex],
                                                 layer.getBaseLayer().Ms);
      }
      fields[idx] = field;
    }
    return fields;
  }

  std::vector<CVector<T>>
  computeExchangeFields(const FDMLayer<T> &layer,
                        const std::vector<CVector<T>> &mags) const {
    std::vector<CVector<T>> fields(mags.size(), CVector<T>());
    const T A = layer.getExchangeStiffness();
    if (A == 0) {
      return fields; // No exchange interaction
    }

    const auto &grid = layer.getGrid();
    const T Ms = layer.getBaseLayer().Ms;
    const T dx2 = grid.dx * grid.dx;
    const T dy2 = grid.dy * grid.dy;
    const T dz2 = grid.dz * grid.dz;

    // Exchange field: H_ex = (A/Ms) * ∇²m
    // Using finite differences: ∇²m_ij = (m_i+1,j + m_i-1,j - 2*m_ij)/dx² +
    //                                     (m_i,j+1 + m_i,j-1 - 2*m_ij)/dy² +
    //                                     (m_i,j+k+1 + m_i,j,k-1 - 2*m_ij)/dz²

    for (std::size_t idx = 0; idx < mags.size(); ++idx) {
      const auto coords = grid.coordinates(idx);
      const unsigned int ix = coords[0];
      const unsigned int iy = coords[1];
      const unsigned int iz = coords[2];

      CVector<T> laplacian(0, 0, 0);

      // X-direction (periodic boundary conditions)
      if (grid.nx > 1) {
        const unsigned int ix_plus = (ix + 1) % grid.nx;
        const unsigned int ix_minus = (ix + grid.nx - 1) % grid.nx;
        const auto idx_plus = grid.index(ix_plus, iy, iz);
        const auto idx_minus = grid.index(ix_minus, iy, iz);
        laplacian = laplacian +
                    (mags[idx_plus] + mags[idx_minus] - mags[idx] * 2.0) / dx2;
      }

      // Y-direction (periodic boundary conditions)
      if (grid.ny > 1) {
        const unsigned int iy_plus = (iy + 1) % grid.ny;
        const unsigned int iy_minus = (iy + grid.ny - 1) % grid.ny;
        const auto idx_plus = grid.index(ix, iy_plus, iz);
        const auto idx_minus = grid.index(ix, iy_minus, iz);
        laplacian = laplacian +
                    (mags[idx_plus] + mags[idx_minus] - mags[idx] * 2.0) / dy2;
      }

      // Z-direction (periodic boundary conditions)
      if (grid.nz > 1) {
        const unsigned int iz_plus = (iz + 1) % grid.nz;
        const unsigned int iz_minus = (iz + grid.nz - 1) % grid.nz;
        const auto idx_plus = grid.index(ix, iy, iz_plus);
        const auto idx_minus = grid.index(ix, iy, iz_minus);
        laplacian = laplacian +
                    (mags[idx_plus] + mags[idx_minus] - mags[idx] * 2.0) / dz2;
      }

      // H_ex = (A/Ms) * ∇²m
      fields[idx] = laplacian * (A / Ms);
    }
    return fields;
  }

  CVector<T> calculateEffectiveField(FDMLayer<T> &layer, const CVector<T> &m,
                                     const CVector<T> &bottom,
                                     const CVector<T> &top,
                                     const CVector<T> &Hdemag,
                                     const CVector<T> &Hexchange, T time) {
    auto &baseLayer = layer.getBaseLayer();
    const CVector<T> Hext = baseLayer.calculateExternalField(time);
    const CVector<T> Hoe = baseLayer.calculateHOeField(time);
    T currentTime = time;
    const CVector<T> HAnis = baseLayer.calculateAnisotropy(m, currentTime);
    const CVector<T> HAnis2 =
        baseLayer.calculateSecondOrderAnisotropy(m, currentTime);
    const CVector<T> HIEC = baseLayer.calculateIEC(time, m, bottom, top);
    const CVector<T> Hidmi = baseLayer.calculateIDMI(time, m, bottom, top);
    const CVector<T> Hdmi = baseLayer.calculateHdmiField(time);
    const CVector<T> Hreserved =
        baseLayer.calculateReservedInteractionField(time);
    return Hext + HAnis + HAnis2 + HIEC + Hidmi + Hoe + Hdmi + Hreserved +
           Hexchange - Hdemag;
  }

  CVector<T> calculateLLG(FDMLayer<T> &layer, const CVector<T> &m,
                          const CVector<T> &bottom, const CVector<T> &top,
                          const CVector<T> &Hdemag, const CVector<T> &Hexchange,
                          T time, T timeStep) {
    const CVector<T> heff =
        calculateEffectiveField(layer, m, bottom, top, Hdemag, Hexchange, time);
    return layer.getBaseLayer().solveLLG(time, m, timeStep, bottom, top, heff);
  }

  void runRK4Step(T t, T timeStep,
                  const std::vector<std::vector<CVector<T>>> &currentMags,
                  const std::vector<std::vector<CVector<T>>> &demagFields,
                  const std::vector<std::vector<CVector<T>>> &exchangeFields,
                  std::vector<std::vector<CVector<T>>> &updatedMags) {
    for (std::size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
      auto &layer = layers[layerIndex];
      const auto &grid = layer.getGrid();
      const auto &layerMags = currentMags[layerIndex];
      const auto &layerDemag = demagFields[layerIndex];
      const auto &layerExchange = exchangeFields[layerIndex];
      for (std::size_t idx = 0; idx < layerMags.size(); ++idx) {
        const auto coords = grid.coordinates(idx);
        const CVector<T> bottom =
            (layerIndex == 0)
                ? CVector<T>()
                : getCoupledMagnetisation(currentMags[layerIndex - 1],
                                          layers[layerIndex - 1].getGrid(),
                                          coords[0], coords[1], coords[2]);
        const CVector<T> top =
            (layerIndex + 1 >= layers.size())
                ? CVector<T>()
                : getCoupledMagnetisation(currentMags[layerIndex + 1],
                                          layers[layerIndex + 1].getGrid(),
                                          coords[0], coords[1], coords[2]);
        const CVector<T> m = layerMags[idx];
        const CVector<T> Hdemag = layerDemag[idx];
        const CVector<T> Hexchange = layerExchange[idx];
        const CVector<T> k1 = calculateLLG(layer, m, bottom, top, Hdemag,
                                           Hexchange, t, timeStep) *
                              timeStep;
        const CVector<T> k2 =
            calculateLLG(layer, m + k1 * 0.5, bottom, top, Hdemag, Hexchange,
                         t + 0.5 * timeStep, timeStep) *
            timeStep;
        const CVector<T> k3 =
            calculateLLG(layer, m + k2 * 0.5, bottom, top, Hdemag, Hexchange,
                         t + 0.5 * timeStep, timeStep) *
            timeStep;
        const CVector<T> k4 = calculateLLG(layer, m + k3, bottom, top, Hdemag,
                                           Hexchange, t + timeStep, timeStep) *
                              timeStep;
        CVector<T> updated = m + (k1 + k2 * 2.0 + k3 * 2.0 + k4) / 6.0;
        updated.normalize();
        updatedMags[layerIndex][idx] = updated;
      }
    }
  }

  void
  runEulerHeunStep(T t, T timeStep,
                   const std::vector<std::vector<CVector<T>>> &currentMags,
                   const std::vector<std::vector<CVector<T>>> &demagFields,
                   const std::vector<std::vector<CVector<T>>> &exchangeFields,
                   std::vector<std::vector<CVector<T>>> &updatedMags) {
    for (std::size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
      auto &layer = layers[layerIndex];
      const auto &grid = layer.getGrid();
      const auto &layerMags = currentMags[layerIndex];
      const auto &layerDemag = demagFields[layerIndex];
      const auto &layerExchange = exchangeFields[layerIndex];
      for (std::size_t idx = 0; idx < layerMags.size(); ++idx) {
        const auto coords = grid.coordinates(idx);
        const CVector<T> bottom =
            (layerIndex == 0)
                ? CVector<T>()
                : getCoupledMagnetisation(currentMags[layerIndex - 1],
                                          layers[layerIndex - 1].getGrid(),
                                          coords[0], coords[1], coords[2]);
        const CVector<T> top =
            (layerIndex + 1 >= layers.size())
                ? CVector<T>()
                : getCoupledMagnetisation(currentMags[layerIndex + 1],
                                          layers[layerIndex + 1].getGrid(),
                                          coords[0], coords[1], coords[2]);
        const CVector<T> m = layerMags[idx];
        const CVector<T> Hdemag = layerDemag[idx];
        const CVector<T> Hexchange = layerExchange[idx];
        const CVector<T> fn =
            calculateLLG(layer, m, bottom, top, Hdemag, Hexchange, t, timeStep);
        const CVector<T> dW =
            layer.getBaseLayer().getStochasticLangevinVector(t, timeStep) +
            layer.getBaseLayer().getOneFVector();
        const CVector<T> gn = layer.getBaseLayer().stochasticTorque(m, dW);
        const CVector<T> mNext = m + gn * std::sqrt(timeStep);
        const CVector<T> gnPrime =
            layer.getBaseLayer().stochasticTorque(mNext, dW);
        CVector<T> updated =
            m + fn * timeStep + 0.5 * (gn + gnPrime) * std::sqrt(timeStep);
        updated.normalize();
        updatedMags[layerIndex][idx] = updated;
      }
    }
  }

  void runHeunStep(T t, T timeStep,
                   const std::vector<std::vector<CVector<T>>> &currentMags,
                   const std::vector<std::vector<CVector<T>>> &demagFields,
                   const std::vector<std::vector<CVector<T>>> &exchangeFields,
                   std::vector<std::vector<CVector<T>>> &updatedMags) {
    for (std::size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
      auto &layer = layers[layerIndex];
      const auto &grid = layer.getGrid();
      const auto &layerMags = currentMags[layerIndex];
      const auto &layerDemag = demagFields[layerIndex];
      const auto &layerExchange = exchangeFields[layerIndex];
      for (std::size_t idx = 0; idx < layerMags.size(); ++idx) {
        const auto coords = grid.coordinates(idx);
        const CVector<T> bottom =
            (layerIndex == 0)
                ? CVector<T>()
                : getCoupledMagnetisation(currentMags[layerIndex - 1],
                                          layers[layerIndex - 1].getGrid(),
                                          coords[0], coords[1], coords[2]);
        const CVector<T> top =
            (layerIndex + 1 >= layers.size())
                ? CVector<T>()
                : getCoupledMagnetisation(currentMags[layerIndex + 1],
                                          layers[layerIndex + 1].getGrid(),
                                          coords[0], coords[1], coords[2]);
        const CVector<T> m = layerMags[idx];
        const CVector<T> Hdemag = layerDemag[idx];
        const CVector<T> Hexchange = layerExchange[idx];
        const CVector<T> fn =
            calculateLLG(layer, m, bottom, top, Hdemag, Hexchange, t, timeStep);
        const CVector<T> dW =
            layer.getBaseLayer().getStochasticLangevinVector(t, timeStep) +
            layer.getBaseLayer().getOneFVector();
        const CVector<T> gn = layer.getBaseLayer().stochasticTorque(m, dW);
        const CVector<T> mNext = m + fn * timeStep + gn * std::sqrt(timeStep);
        const CVector<T> fnPrime =
            calculateLLG(layer, mNext, bottom, top, Hdemag, Hexchange,
                         t + timeStep, timeStep);
        const CVector<T> gnPrime =
            layer.getBaseLayer().stochasticTorque(mNext, dW);
        CVector<T> updated = m + 0.5 * timeStep * (fn + fnPrime) +
                             0.5 * (gn + gnPrime) * std::sqrt(timeStep);
        updated.normalize();
        updatedMags[layerIndex][idx] = updated;
      }
    }
  }

public:
  FDMJunction() = default;
  explicit FDMJunction(const std::vector<FDMLayer<T>> &layersToSet)
      : layers(layersToSet) {
    layerNo = static_cast<unsigned int>(layers.size());
    if (layerNo == 0) {
      throw std::invalid_argument("Passed a zero length Layer vector!");
    }
    const auto &grid = layers[0].getGrid();
    for (const auto &layer : layers) {
      if (layer.getGrid().nx != grid.nx || layer.getGrid().ny != grid.ny) {
        throw std::runtime_error("All FDM layers must share nx and ny");
      }
    }
  }

  const std::vector<std::string> getLayerIds() const {
    std::vector<std::string> ids;
    std::transform(layers.begin(), layers.end(), std::back_inserter(ids),
                   [](const FDMLayer<T> &layer) { return layer.getId(); });
    return ids;
  }

  FDMLayer<T> &getLayer(const std::string &layerId) {
    for (auto &layer : layers) {
      if (layer.getId() == layerId) {
        return layer;
      }
    }
    throw std::runtime_error("Failed to find a layer with id: " + layerId);
  }

  const std::vector<CVector<T>> &
  getLayerMagnetisationGrid(const std::string &layerId) {
    return getLayer(layerId).getMagnetisationGrid();
  }

  void setLayerMagnetisation(const std::string &layerId,
                             const CVector<T> &mag) {
    getLayer(layerId).setMagnetisation(mag);
  }

  void setLayerMagnetisationGrid(const std::string &layerId,
                                 const std::vector<CVector<T>> &mags) {
    getLayer(layerId).setMagnetisationGrid(mags);
  }

  void setLayerDemagTensor(
      const std::string &layerId,
      const std::vector<typename FDMLayer<T>::DemagTensorCell> &tensor) {
    getLayer(layerId).setDemagTensor(tensor);
  }

  void setLayerExchangeStiffness(const std::string &layerId, T A) {
    getLayer(layerId).setExchangeStiffness(A);
  }

  void setLayerExternalFieldDriver(const std::string &layerId,
                                   const AxialDriver<T> &driver) {
    axiallayerSetter(layerId, &Layer<T>::setExternalFieldDriver, driver);
  }

  void setLayerOerstedFieldDriver(const std::string &layerId,
                                  const AxialDriver<T> &driver) {
    axiallayerSetter(layerId, &Layer<T>::setOerstedFieldDriver, driver);
  }

  void setLayerHdmiDriver(const std::string &layerId,
                          const AxialDriver<T> &driver) {
    axiallayerSetter(layerId, &Layer<T>::setHdmiDriver, driver);
  }

  void setLayerCurrentDriver(const std::string &layerId,
                             const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerId, &Layer<T>::setCurrentDriver, driver);
  }

  void setLayerAnisotropyDriver(const std::string &layerId,
                                const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerId, &Layer<T>::setAnisotropyDriver, driver);
  }

  void setLayerSecondOrderAnisotropyDriver(const std::string &layerId,
                                           const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerId, &Layer<T>::setSecondOrderAnisotropyDriver,
                      driver);
  }

  void setLayerTemperatureDriver(const std::string &layerId,
                                 const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerId, &Layer<T>::setTemperatureDriver, driver);
  }

  void setLayerFieldLikeTorqueDriver(const std::string &layerId,
                                     const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerId, &Layer<T>::setFieldLikeTorqueDriver, driver);
  }

  void setLayerDampingLikeTorqueDriver(const std::string &layerId,
                                       const ScalarDriver<T> &driver) {
    scalarlayerSetter(layerId, &Layer<T>::setDampingLikeTorqueDriver, driver);
  }

  void setLayerReferenceLayer(const std::string &layerId,
                              const CVector<T> &reference) {
    for (auto &layer : layers) {
      if (layer.getId() == layerId || layerId == "all") {
        layer.getBaseLayer().setReferenceLayer(reference);
      }
    }
  }

  void setLayerReferenceType(const std::string &layerId, Reference reference) {
    for (auto &layer : layers) {
      if (layer.getId() == layerId || layerId == "all") {
        layer.getBaseLayer().setReferenceLayer(reference);
      }
    }
  }

  void setIECDriver(const std::string &bottomLayer, const std::string &topLayer,
                    const ScalarDriver<T> &driver) {
    setCouplingDriver(bottomLayer, topLayer, driver, &Layer<T>::setIECDriverTop,
                      &Layer<T>::setIECDriverBottom);
  }

  void setQuadIECDriver(const std::string &bottomLayer,
                        const std::string &topLayer,
                        const ScalarDriver<T> &driver) {
    setCouplingDriver(bottomLayer, topLayer, driver,
                      &Layer<T>::setQuadIECDriverTop,
                      &Layer<T>::setQuadIECDriverBottom);
  }

  void setIDMIDriver(const std::string &bottomLayer,
                     const std::string &topLayer,
                     const AxialDriver<T> &driver) {
    bool found = false;
    for (unsigned int i = 0; i < layerNo - 1; i++) {
      const auto &current = layers[i].getId();
      const auto &next = layers[i + 1].getId();
      if ((current == bottomLayer && next == topLayer) ||
          (current == topLayer && next == bottomLayer)) {
        layers[i].getBaseLayer().setIDMIDriverTop(driver);
        layers[i + 1].getBaseLayer().setIDMIDriverBottom(driver);
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

  void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11,
                     bool verbose = false, SolverMode mode = RK4) {
    if (timeStep > writeFrequency) {
      throw std::runtime_error(
          "The time step cannot be larger than write frequency!");
    }
    const unsigned int totalIterations =
        static_cast<unsigned int>(totalTime / timeStep);
    SolverMode solverMode = resolveSolverMode(mode, totalIterations);
    if (solverMode == DORMAND_PRINCE) {
      throw std::runtime_error("Dormand-Prince is not supported for FDM yet");
    }
    std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();

    std::vector<std::vector<CVector<T>>> currentMags(layers.size());
    std::vector<std::vector<CVector<T>>> updatedMags(layers.size());
    for (std::size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
      currentMags[layerIndex] = layers[layerIndex].getMagnetisationGrid();
      updatedMags[layerIndex] = currentMags[layerIndex];
    }

    for (unsigned int i = 0; i < totalIterations; i++) {
      const T t = static_cast<T>(i) * timeStep;
      std::vector<std::vector<CVector<T>>> demagFields;
      std::vector<std::vector<CVector<T>>> exchangeFields;
      demagFields.reserve(layers.size());
      exchangeFields.reserve(layers.size());
      for (std::size_t layerIndex = 0; layerIndex < layers.size();
           ++layerIndex) {
        demagFields.push_back(
            computeDemagFields(layers[layerIndex], currentMags[layerIndex]));
        exchangeFields.push_back(
            computeExchangeFields(layers[layerIndex], currentMags[layerIndex]));
      }

      if (solverMode == EULER_HEUN) {
        runEulerHeunStep(t, timeStep, currentMags, demagFields, exchangeFields,
                         updatedMags);
      } else if (solverMode == HEUN) {
        runHeunStep(t, timeStep, currentMags, demagFields, exchangeFields,
                    updatedMags);
      } else {
        runRK4Step(t, timeStep, currentMags, demagFields, exchangeFields,
                   updatedMags);
      }

      currentMags.swap(updatedMags);
    }

    for (std::size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
      layers[layerIndex].setMagnetisationGrid(currentMags[layerIndex]);
    }

    if (verbose) {
      std::chrono::steady_clock::time_point end =
          std::chrono::steady_clock::now();
      std::cout << "Steps in simulation: " << totalIterations << std::endl;
      std::cout << "Simulation time = "
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin)
                       .count()
                << "[s]" << std::endl;
    }
  }
};

#endif // CORE_FDM_HPP_
