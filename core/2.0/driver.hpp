#ifndef DRIVERS_H
#define DRIVERS_H

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <cassert>  // for assert
#include <stdlib.h> // for abs
#define _USE_MATH_DEFINES
#include "cvector.hpp" // for CVector
#include <cmath>       // for M_PI
#include <memory>      // for shared_ptr
#include <stdexcept>   // for runtime_error
#include <utility>     // for move
#include <vector>      // for vector

enum UpdateType {
  constant,
  pulse,
  sine,
  step,
  posine,
  halfsine,
  trapezoid,
  gaussimpulse,
  gaussstep
};

// Abstract base class
template <typename T> class Driver {
protected:
  T constantValue;

public:
  Driver(T constantValue = 0.0) : constantValue(constantValue) {}
  virtual ~Driver() = default;

  // Pure virtual method - must be implemented by derived classes
  virtual T getCurrentScalarValue(const T &time) const = 0;

  // Common methods
  void setConstantValue(const T &val) { this->constantValue = val; }
  T getConstantValue() const { return constantValue; }
};

// Default implementation - returns constant value only
template <typename T> class ConstantDriver : public Driver<T> {
public:
  explicit ConstantDriver(T constantValue = 0.0) : Driver<T>(constantValue) {}

  T getCurrentScalarValue(const T &time) const override {
    return this->constantValue;
  }
};

// Alias for ConstantDriver
template <typename T> class NullDriver : public ConstantDriver<T> {
public:
  explicit NullDriver(T constantValue = 0.0)
      : ConstantDriver<T>(constantValue) {}
};

// Sine wave driver
template <typename T> class SineDriver : public Driver<T> {
private:
  T amplitude;
  T frequency;
  T phase;

public:
  explicit SineDriver(T constantValue = 0.0, T amplitude = 0.0,
                      T frequency = 0.0, T phase = 0.0)
      : Driver<T>(constantValue), amplitude(amplitude), frequency(frequency),
        phase(phase) {
    if (frequency < 0) {
      throw std::runtime_error("Sine driver requires positive frequency");
    }
  }

  T getCurrentScalarValue(const T &time) const override {
    return this->constantValue +
           amplitude * sin(2 * M_PI * time * frequency + phase);
  }

  // Getters and setters
  T getAmplitude() const { return amplitude; }
  void setAmplitude(T value) { amplitude = value; }

  T getFrequency() const { return frequency; }
  void setFrequency(T value) {
    if (value < 0) {
      throw std::runtime_error("Frequency must be positive");
    }
    frequency = value;
  }

  T getPhase() const { return phase; }
  void setPhase(T value) { phase = value; }
};

// Pulse train driver
template <typename T> class PulseDriver : public Driver<T> {
private:
  T amplitude;
  T period;
  T cycle;

  T pulseTrain(T time) const {
    const int n = static_cast<int>(time / period);
    const T dT = cycle * period;
    const T nT = n * period;
    if (nT <= time && time <= (nT + dT)) {
      return amplitude;
    } else {
      return 0;
    }
  }

public:
  explicit PulseDriver(T constantValue = 0.0, T amplitude = 0.0, T period = 0.0,
                       T cycle = 0.5)
      : Driver<T>(constantValue), amplitude(amplitude), period(period),
        cycle(cycle) {
    if (period <= 0 || cycle < 0 || cycle > 1) {
      throw std::runtime_error(
          "Pulse driver requires positive period and cycle in [0,1]");
    }
  }

  T getCurrentScalarValue(const T &time) const override {
    return this->constantValue + pulseTrain(time);
  }

  // Getters and setters
  T getAmplitude() const { return amplitude; }
  void setAmplitude(T value) { amplitude = value; }

  T getPeriod() const { return period; }
  void setPeriod(T value) {
    if (value <= 0) {
      throw std::runtime_error("Period must be positive");
    }
    period = value;
  }

  T getCycle() const { return cycle; }
  void setCycle(T value) {
    if (value < 0 || value > 1) {
      throw std::runtime_error("Cycle must be in range [0,1]");
    }
    cycle = value;
  }
};

// Step driver
template <typename T> class StepDriver : public Driver<T> {
private:
  T amplitude;
  T timeStart;
  T timeStop;

  T stepUpdate(T time) const {
    if (time >= timeStart && time <= timeStop) {
      return amplitude;
    } else {
      return 0.0;
    }
  }

public:
  StepDriver(T constantValue = 0.0, T amplitude = 0.0, T timeStart = 0.0,
             T timeStop = 0.0)
      : Driver<T>(constantValue), amplitude(amplitude), timeStart(timeStart),
        timeStop(timeStop) {
    if (timeStop <= timeStart) {
      throw std::runtime_error("Step driver requires timeStop > timeStart");
    }
  }

  T getCurrentScalarValue(const T &time) const override {
    return this->constantValue + stepUpdate(time);
  }

  // Getters and setters
  T getAmplitude() const { return amplitude; }
  void setAmplitude(T value) { amplitude = value; }

  T getTimeStart() const { return timeStart; }
  void setTimeStart(T value) {
    if (value >= timeStop) {
      throw std::runtime_error("timeStart must be less than timeStop");
    }
    timeStart = value;
  }

  T getTimeStop() const { return timeStop; }
  void setTimeStop(T value) {
    if (value <= timeStart) {
      throw std::runtime_error("timeStop must be greater than timeStart");
    }
    timeStop = value;
  }
};

// Scalar driver factory - creates appropriate driver types
template <typename T> class ScalarDriver : public Driver<T> {
public:
  static std::shared_ptr<Driver<T>> getConstantDriver(T constantValue) {
    return std::make_shared<ConstantDriver<T>>(constantValue);
  }

  static std::shared_ptr<Driver<T>> getSineDriver(T constantValue, T amplitude,
                                                  T frequency, T phase) {
    return std::make_shared<SineDriver<T>>(constantValue, amplitude, frequency,
                                           phase);
  }

  static std::shared_ptr<Driver<T>> getPulseDriver(T constantValue, T amplitude,
                                                   T period, T cycle) {
    return std::make_shared<PulseDriver<T>>(constantValue, amplitude, period,
                                            cycle);
  }

  static std::shared_ptr<Driver<T>> getStepDriver(T constantValue, T amplitude,
                                                  T timeStart, T timeStop) {
    return std::make_shared<StepDriver<T>>(constantValue, amplitude, timeStart,
                                           timeStop);
  }

  // Implement the pure virtual function from Driver
  T getCurrentScalarValue(const T &time) const override {
    return this->constantValue;
  }
};

// Axial driver - manages three scalar drivers for x, y, z components
template <typename T> class AxialDriver : public Driver<T> {
private:
  std::vector<std::shared_ptr<Driver<T>>> drivers;

public:
  AxialDriver() {
    this->drivers = {std::make_shared<ConstantDriver<T>>(),
                     std::make_shared<ConstantDriver<T>>(),
                     std::make_shared<ConstantDriver<T>>()};
  }

  AxialDriver(std::shared_ptr<Driver<T>> x, std::shared_ptr<Driver<T>> y,
              std::shared_ptr<Driver<T>> z) {
    this->drivers = {x, y, z};
  }

  explicit AxialDriver(const CVector<T> &xyz) {
    this->drivers = {ScalarDriver<T>::getConstantDriver(xyz.x),
                     ScalarDriver<T>::getConstantDriver(xyz.y),
                     ScalarDriver<T>::getConstantDriver(xyz.z)};
  }

  static AxialDriver getVectorAxialDriver(T x, T y, T z) {
    return AxialDriver(CVector<T>(x, y, z));
  }

  CVector<T> getCurrentAxialDrivers(const T &time) const {
    return CVector<T>(this->drivers[0]->getCurrentScalarValue(time),
                      this->drivers[1]->getCurrentScalarValue(time),
                      this->drivers[2]->getCurrentScalarValue(time));
  }

  CVector<T> getConstantValues() const {
    return CVector<T>(this->drivers[0]->getConstantValue(),
                      this->drivers[1]->getConstantValue(),
                      this->drivers[2]->getConstantValue());
  }

  CVector<T> getUnitAxis() const {
    CVector<T> values = getConstantValues();
    return CVector<T>(values.x != 0.0 ? values.x / std::abs(values.x) : 0.0,
                      values.y != 0.0 ? values.y / std::abs(values.y) : 0.0,
                      values.z != 0.0 ? values.z / std::abs(values.z) : 0.0);
  }

  // Implement the pure virtual function from Driver
  T getCurrentScalarValue(const T &time) const override {
    throw std::runtime_error("Cannot get scalar value from AxialDriver");
  }
};

#endif // DRIVERS_H
