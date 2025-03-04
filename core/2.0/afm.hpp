#ifndef CORE_AFM_HPP_
#define CORE_AFM_HPP_

#include "abstract.hpp"
#include "drivers.hpp"
#include "fm.hpp"

template <typename T> class LayerAFM : public Layer<T> {
public:
  // Second sublattice magnetization
  CVector<T> mag1, mag2;

  // AFM exchange coupling constant between sublattices
  Driver<T> J_AFM = ScalarDriver<T>::getConstantDriver(0.0);

  // Damping for the second sublattice (can be different)

  /**
   * @brief Creates a layer with Antiferromagnetic properties
   *
   * @param id Layer identifier
   * @param mag1 Initial magnetization vector for first sublattice
   * @param mag2 Initial magnetization vector for second sublattice
   * @param anis Anisotropy vector
   * @param Ms Saturation magnetization
   * @param thickness Layer thickness
   * @param cellSurface Cell surface area
   * @param demagTensor Demagnetization tensor
   * @param damping Damping parameter for first sublattice
   * @param J_AFM Exchange coupling constant between sublattices (negative for
   * AFM)
   */
  explicit LayerAFM(const std::string &id, const CVector<T> &mag1,
                    const CVector<T> &mag2, const CVector<T> &anis, T Ms,
                    T thickness, T cellSurface,
                    const std::vector<CVector<T>> &demagTensor, T damping,
                    T J_AFM)
      : Layer<T>(id, CVector<T>(0, 0, 0), anis, Ms, thickness, cellSurface,
                 demagTensor, damping) {
    // Initialize second sublattice
    this->mag1 = mag1;
    this->mag2 = mag2;
    this->mag1.normalize();
    this->mag2.normalize();
    // Set AFM parameters
  }

  CVector<T> getMagnetisation() const override { return getNetMagnetisation(); }

  void setAFMExchangeCoupling(std::shared_ptr<Driver<T>> J_AFM) {
    this->J_AFM = J_AFM;
  }

  // Get magnetization of the second sublattice
  CVector<T> getMagnetisation(unsigned int sublattice) const {
    if (sublattice == 0) {
      return this->mag1;
    } else if (sublattice == 1) {
      return this->mag2;
    } else {
      throw std::runtime_error("Invalid sublattice index");
    }
  }

  void setMagnetisation(const CVector<T> &newMag) override {
    throw std::runtime_error(
        "AFM setMagnetisation makes no sense. "
        "Use setMagnetisation(sublattice, newMag) instead.");
  }

  void setMagnetisation(unsigned int sublattice, const CVector<T> &newMag) {
    if (sublattice == 1) {
      this->mag1 = newMag;
    } else if (sublattice == 2) {
      this->mag2 = newMag;
    } else {
      throw std::runtime_error("Invalid sublattice index");
    }
  }

  // Calculate the exchange field between sublattices - Common method
  CVector<T> calculateAFMExchangeField(const T &time,
                                       const CVector<T> &targetMag) {
    // J_AFM should be negative for antiferromagnetic coupling
    const T J_AFM_value = this->J_AFM.getCurrentScalarValue(time);
    return targetMag * (J_AFM_value / (this->Ms * this->thickness));
  }

  // Calculate the effective field for any sublattice
  const CVector<T>
  calculateHeffCommon(const T &time, const T &timeStep, const CVector<T> &m1,
                      const CVector<T> &m2, const CVector<T> &bottom,
                      const CVector<T> &top,
                      const CVector<T> &Hfluctuation = CVector<T>(),
                      const CVector<T> &Hdipole = CVector<T>()) {
    // Calculate the standard effective field from parent class
    // passed by const reference to avoid modifying the original magnetization
    CVector<T> Heff = Layer<T>::calculateHeff(time, timeStep, m1, bottom, top,
                                              Hfluctuation, Hdipole);
    // Add the AFM exchange field contribution
    CVector<T> Hafm = calculateAFMExchangeField(time, m2);

    return Heff + Hafm;
  }

  // Calculate LLG for second sublattice
  const CVector<T>
  calculateLLGSublattice(const T &time, const T &timeStep, const CVector<T> &m1,
                         const CVector<T> &m2,
                         const CVector<T> &bottom = CVector<T>(),
                         const CVector<T> &top = CVector<T>()) {
    // Get effective field for second sublattice
    CVector<T> Heff = calculateHeffCommon(time, timeStep, m1, m2, bottom, top);

    // Use common method with second sublattice damping
    return Layer<T>::calculateLLG(time, timeStep, m1, bottom, top);
  }

  // Override the main LLG calculation method
  const CVector<T> calculateLLG(const T &time, const T &timeStep,
                                const CVector<T> &m,
                                const CVector<T> &bottom = CVector<T>(),
                                const CVector<T> &top = CVector<T>()) override {
    // Calculate dynamics for first sublattice
    throw std::runtime_error("AFM LLG calculation not implemented");
    return CVector<T>(0, 0, 0);
  }

  // RK4 step implementation that handles both sublattices
  void rk4_step(const T &time, const T &timeStep, const CVector<T> &bottom,
                const CVector<T> &top) {
    // For the first sublattice (same as parent class)
    CVector<T> m1_t = this->mag1;
    CVector<T> m2_t = this->mag2;

    CVector<T> k1_1 =
        calculateLLGSublattice(time, timeStep, m1_t, m2_t, bottom, top) *
        timeStep;
    CVector<T> k2_1 =
        calculateLLGSublattice(time + timeStep / 2, timeStep, m1_t + k1_1 / 2,
                               m2_t, bottom, top) *
        timeStep;
    CVector<T> k3_1 =
        calculateLLGSublattice(time + timeStep / 2, timeStep, m1_t + k2_1 / 2,
                               m2_t, bottom, top) *
        timeStep;
    CVector<T> k4_1 = calculateLLGSublattice(time + timeStep, timeStep,
                                             m1_t + k3_1, m2_t, bottom, top) *
                      timeStep;

    m1_t += (k1_1 + k2_1 * 2 + k3_1 * 2 + k4_1) / 6;
    m1_t.normalize();

    // For the second sublattice
    CVector<T> k1_2 =
        calculateLLGSublattice(time, timeStep, m2_t, m1_t, bottom, top) *
        timeStep;
    CVector<T> k2_2 =
        calculateLLGSublattice(time + timeStep / 2, timeStep, m2_t + k1_2 / 2,
                               m1_t, bottom, top) *
        timeStep;
    CVector<T> k3_2 =
        calculateLLGSublattice(time + timeStep / 2, timeStep, m2_t + k2_2 / 2,
                               m1_t, bottom, top) *
        timeStep;
    CVector<T> k4_2 = calculateLLGSublattice(time + timeStep, timeStep,
                                             m2_t + k3_2, m1_t, bottom, top) *
                      timeStep;

    m2_t += (k1_2 + k2_2 * 2 + k3_2 * 2 + k4_2) / 6;
    m2_t.normalize();

    // Update both sublattices
    this->mag1 = m1_t;
    this->mag2 = m2_t;

    if (isnan(this->mag.x) || isnan(this->mag2.x)) {
      throw std::runtime_error("NAN magnetisation in AFM layer");
    }
  }

  // Getter for net magnetization (difference between sublattices)
  CVector<T> getNetMagnetisation() const {
    return (this->mag1 + this->mag2) / 2;
  }

  // Getter for Néel vector (sum of sublattices)
  CVector<T> getNeelVector() const { return (this->mag1 - this->mag2) / 2; }
};

#endif // CORE_AFM_HPP_
