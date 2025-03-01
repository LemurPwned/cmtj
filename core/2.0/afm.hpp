#ifndef CORE_AFM_HPP_
#define CORE_AFM_HPP_

#include "abstract.hpp"
#include "drivers.hpp"
#include "fm_layer.hpp"

template <typename T> class AFMLayer : public AbstractLayer<T> {
private:
  Layer<T> sublatticeA;
  Layer<T> sublatticeB;
  T exchangeCoupling; // Exchange coupling between sublattices

public:
  AFMLayer(const std::string &id, const CVector<T> &magA,
           const CVector<T> &magB, const CVector<T> &anis, T Ms, T thickness,
           T cellSurface, const std::vector<CVector<T>> &demagTensor, T damping,
           T exchangeCoupling)
      : sublatticeA(id + "_A", magA, anis, Ms, thickness, cellSurface,
                    demagTensor, damping),
        sublatticeB(id + "_B", magB, anis, Ms, thickness, cellSurface,
                    demagTensor, damping),
        exchangeCoupling(exchangeCoupling) {
    // Initialize sublattices with antiparallel alignment
    if (magA.length() == 0 || magB.length() == 0) {
      throw std::runtime_error(
          "Initial magnetisations cannot be zero vectors!");
    }
  }

  // Implementation of pure virtual functions
  void setMagnetisation(CVector<T> &mag) override {
    // Set opposite magnetizations for AFM sublattices
    sublatticeA.setMagnetisation(mag);
    sublatticeB.setMagnetisation(mag * -1.0);
  }

  CVector<T> getMagnetisation() const override {
    // Return net magnetization (usually near zero for AFM)
    return (sublatticeA.mag + sublatticeB.mag) * 0.5;
  }

  void setReferenceLayer(const CVector<T> &reference) override {
    sublatticeA.setReferenceLayer(reference);
    sublatticeB.setReferenceLayer(reference);
  }

  void setReferenceLayer(Reference reference) override {
    sublatticeA.setReferenceLayer(reference);
    sublatticeB.setReferenceLayer(reference);
  }

  CVector<T> getReferenceLayer() const override {
    return sublatticeA.getReferenceLayer();
  }

  Reference getReferenceType() const override {
    return sublatticeA.getReferenceType();
  }

  // Driver setters implementation
  void setTemperatureDriver(const ScalarDriver<T> &driver) override {
    sublatticeA.setTemperatureDriver(driver);
    sublatticeB.setTemperatureDriver(driver);
  }

  void setCurrentDriver(const ScalarDriver<T> &driver) override {
    sublatticeA.setCurrentDriver(driver);
    sublatticeB.setCurrentDriver(driver);
  }

  void setAnisotropyDriver(const ScalarDriver<T> &driver) override {
    sublatticeA.setAnisotropyDriver(driver);
    sublatticeB.setAnisotropyDriver(driver);
  }

  void setExternalFieldDriver(const AxialDriver<T> &driver) override {
    sublatticeA.setExternalFieldDriver(driver);
    sublatticeB.setExternalFieldDriver(driver);
  }

  // Solver methods implementation
  CVector<T>
  calculateHeff(T time, T timeStep, const CVector<T> &stepMag,
                const CVector<T> &bottom, const CVector<T> &top,
                const CVector<T> &Hfluctuation = CVector<T>()) override {
    // Calculate effective field including AFM exchange interaction
    CVector<T> HeffA = sublatticeA.calculateHeff(
        time, timeStep, sublatticeA.mag, bottom, top, Hfluctuation);
    CVector<T> HeffB = sublatticeB.calculateHeff(
        time, timeStep, sublatticeB.mag, bottom, top, Hfluctuation);

    // Add exchange coupling between sublattices
    HeffA = HeffA - exchangeCoupling * sublatticeB.mag;
    HeffB = HeffB - exchangeCoupling * sublatticeA.mag;

    // Return average effective field
    return (HeffA + HeffB) * 0.5;
  }

  void rk4_step(T time, T timeStep, const CVector<T> &bottom,
                const CVector<T> &top) override {
    // Perform RK4 integration for both sublattices
    sublatticeA.rk4_step(time, timeStep, bottom, top);
    sublatticeB.rk4_step(time, timeStep, bottom, top);
  }
};

#endif // CORE_AFM_HPP_
