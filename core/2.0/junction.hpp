#ifndef CORE_JUNCTION_HPP_
#define CORE_JUNCTION_HPP_

#include "abstract.hpp"
#include "fm.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

template <typename T> class FMJunction : public AbstractJunction<T> {
protected:
  T Rap;
  T Rp;
  MRmode mode = MRmode::CLASSIC;

public:
  void setMRMode(MRmode newMode) { mode = newMode; }
  void setRp(T value) { Rp = value; }
  void setRap(T value) { Rap = value; }

  /**
   * @brief Creates a ferromagnetic junction
   *
   * @param layers Vector of layers
   */
  explicit FMJunction(
      const std::vector<std::shared_ptr<AbstractLayer<T>>> &layers)
      : AbstractJunction<T>(layers) {}

  FMJunction(const std::vector<std::shared_ptr<AbstractLayer<T>>> &layers, T Rp,
             T Rap)
      : AbstractJunction<T>(layers), Rp(Rp), Rap(Rap) {}

  T calculateMagnetoresistance(T cosTheta) {
    return this->Rp + (((this->Rap - this->Rp) / 2.0) * (1.0 - cosTheta));
  }

  std::vector<T> getMagnetoresistance() override {
    // this is classical bilayer case
    std::cout << "MR mode: " << this->mode << ", layerNo: " << this->layerNo
              << std::endl;
    if (this->mode == MRmode::CLASSIC && this->layerNo == 2) {
      return {calculateMagnetoresistance(
          c_dot<T>(this->layers[0]->getMagnetisation(),
                   this->layers[1]->getMagnetisation()))};
    }
    // this is the case when we use the pinning layer
    else if (this->mode == MRmode::CLASSIC && this->layerNo == 1) {
      return {calculateMagnetoresistance(
          c_dot<T>(this->layers[0]->getMagnetisation(),
                   this->layers[0]->getReferenceLayer()))};
    } else {
      throw std::runtime_error(
          "Magnetisation calculation is not supported for this structure!");
    }
  }

  void setLayerTemperatureDriver(const std::string &layerID,
                                 std::shared_ptr<Driver<T>> driver) {
    this->scalarlayerSetter(layerID, &Layer<T>::setTemperatureDriver, driver);
  }

  void setLayerCurrentDriver(const std::string &layerID,
                             std::shared_ptr<Driver<T>> driver) {
    this->scalarlayerSetter(layerID, &Layer<T>::setCurrentDriver, driver);
  }

  void setLayerAnisotropyDriver(const std::string &layerID,
                                std::shared_ptr<Driver<T>> driver) {
    this->scalarlayerSetter(layerID, &Layer<T>::setAnisotropyDriver, driver);
  }

  void setLayerOerstedFieldDriver(const std::string &layerID,
                                  std::shared_ptr<AxialDriver<T>> driver) {
    this->axiallayerSetter(layerID, &Layer<T>::setOerstedFieldDriver, driver);
  }

  void setLayerExternalFieldDriver(const std::string &layerID,
                                   std::shared_ptr<AxialDriver<T>> driver) {
    this->axiallayerSetter(layerID, &Layer<T>::setExternalFieldDriver, driver);
  }
};

#endif // CORE_JUNCTION_HPP_
