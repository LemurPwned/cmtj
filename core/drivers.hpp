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
#include <pybind11/pybind11.h>
#include <stdexcept> // for runtime_error
#include <utility>   // for move
#include <vector>    // for vector

enum UpdateType {
  constant,
  pulse,
  sine,
  step,
  posine,
  halfsine,
  trapezoid,
  gaussimpulse,
  gaussstep,
  custom = 100
};

template <typename T> class Driver {
protected:
  // if the user wants to update, let them do that
  T constantValue, amplitude, frequency, phase, period, cycle, timeStart,
      timeStop;
  UpdateType update;

public:
  Driver() {
    this->constantValue = 0.0;
    this->amplitude = 0.0;
    this->frequency = 0.0;
    this->phase = 0.0;
    this->period = 0.0;
    this->timeStart = 0.0;
    this->cycle = 0.0;
    this->timeStop = 0.0;
    this->update = constant;
  };
  Driver(UpdateType update, T constantValue, T amplitude, T frequency, T phase,
         T period, T cycle, T timeStart, T timeStop)
      : constantValue(constantValue), amplitude(amplitude),
        frequency(frequency), phase(phase), period(period), cycle(cycle),
        timeStart(timeStart), timeStop(timeStop), update(update)

  {}
  virtual T getCurrentScalarValue(T &time) { return 0; };
  virtual ~Driver() = default;

  void setConstantValue(const T &val) { this->constantValue = val; }

  void phaseShift(const T &phase) { this->phase += phase; }
};

template <typename T> class ScalarDriver : public Driver<T> {

private:
  T edgeTime = 0;
  T steadyTime = 0;
  pybind11::function m_callback;

protected:
  T stepUpdate(T amplitude, T time, T timeStart, T timeStop) {
    if (time >= timeStart && time <= timeStop) {
      return amplitude;
    }
    return 0.0;
  }
  T pulseTrain(T amplitude, T time, T period, T cycle) {
    const int n = static_cast<int>(time / period);
    const T dT = cycle * period;
    const T nT = n * period;
    if (nT <= time && time <= (nT + dT)) {
      return amplitude;
    }
    return 0.0;
  }

  T trapezoidalUpdate(T amplitude, T time, T timeStart, T edgeTime,
                      T steadyTime) {
    if (time < timeStart) {
      return 0;
    }
    // growth
    else if (time <= timeStart + edgeTime) {
      return (amplitude / edgeTime) * (time - timeStart);
    }
    // steady
    else if (time <= timeStart + edgeTime + steadyTime) {
      return amplitude;
    }
    // decay
    else if (time <= timeStart + 2 * edgeTime + steadyTime) {
      return amplitude - (amplitude / edgeTime) *
                             (time - (timeStart + edgeTime + steadyTime));
    }
    return 0;
  }

public:
  explicit ScalarDriver(UpdateType update = constant, T constantValue = 0,
                        T amplitude = 0, T frequency = -1, T phase = 0,
                        T period = -1, T cycle = -1, T timeStart = -1,
                        T timeStop = -1, T edgeTime = -1, T steadyTime = -1,
                        pybind11::function m_callback = pybind11::function())
      : Driver<T>(update, constantValue, amplitude, frequency, phase, period,
                  cycle, timeStart, timeStop) {
    this->edgeTime = edgeTime;
    this->steadyTime = steadyTime;
    if (update == pulse && ((period == -1) || (cycle == -1))) {
      throw std::runtime_error("Selected pulse train driver type but either "
                               "period or cycle were not set");
    } else if (update == sine && (frequency == -1)) {
      throw std::runtime_error(
          "Selected sine driver type but frequency was not set");
    }
    this->m_callback = m_callback;
  }

  /**
   * Constant driver produces a constant signal of a fixed amplitude.
   * @param constantValue: constant value of the driver (constant
   * offset/amplitude)
   */
  static ScalarDriver getConstantDriver(T constantValue) {
    return ScalarDriver(constant, constantValue);
  }

  /**
   * Produces a square pulse of certain period and cycle
   * @param constantValue: offset (vertical) of the pulse. The pulse amplitude
   * will be added to this.
   * @param amplitude: amplitude of the pulse signal
   * @param period: period of the signal in seconds
   * @param cycle: duty cycle of the signal -- a fraction between [0 and 1].
   */
  static ScalarDriver getPulseDriver(T constantValue, T amplitude, T period,
                                     T cycle) {
    return ScalarDriver(pulse, constantValue, amplitude, -1, -1, period, cycle);
  }

  /**
   * Produces a sinusoidal signal with some offset (constantValue), amplitude
   * frequency and phase offset.
   * @param constantValue: vertical offset. The sine will oscillate around this
   * value.
   * @param amplitude: amplitude of the sine wave
   * @param frequency: frequency of the sine
   * @param phase: phase of the sine in radians.
   */
  static ScalarDriver getSineDriver(T constantValue, T amplitude, T frequency,
                                    T phase) {
    return ScalarDriver(sine, constantValue, amplitude, frequency, phase);
  }

  /**
   * Produces a positive sine signal with some offset (constantValue), amplitude
   * frequency and phase offset.
   * @param constantValue: vertical offset. The sine will oscillate around this
   * value.
   * @param amplitude: amplitude of the sine wave
   * @param frequency: frequency of the sine
   * @param phase: phase of the sine in radians.
   */
  static ScalarDriver getPosSineDriver(T constantValue, T amplitude,
                                       T frequency, T phase) {
    return ScalarDriver(posine, constantValue, amplitude, frequency, phase);
  }

  static ScalarDriver getHalfSineDriver(T constantValue, T amplitude,
                                        T frequency, T phase) {
    return ScalarDriver(halfsine, constantValue, amplitude, frequency, phase);
  }
  /**
   * Get a step driver. It has amplitude between timeStart and timeStop and 0
   * elsewhere
   * @param constantValue: offset of the pulse (vertical)
   * @param amplitude: amplitude that is added on top of the constantValue
   * @param timeStart: start of the pulse
   * @param timeStop: when the pulse ends
   */
  static ScalarDriver getStepDriver(T constantValue, T amplitude, T timeStart,
                                    T timeStop) {
    if (timeStop <= timeStart) {
      throw std::runtime_error("Start time cannot be later than stop time!");
    }
    return ScalarDriver(step, constantValue, amplitude, -1, -1, -1, -1,
                        timeStart, timeStop);
  }

  /**
   * Get a trapezoidal driver. It has amplitude between timeStart and timeStop
   * and 0 elsewhere
   * @param constantValue: offset of the pulse (vertical)
   * @param amplitude: amplitude that is added on top of the constantValue
   * @param timeStart: start of the pulse
   * @param edgeTime: time it takes to reach the maximum amplitude
   * @param steadyTime: time it spends in a steady state
   */
  static ScalarDriver getTrapezoidDriver(T constantValue, T amplitude,
                                         T timeStart, T edgeTime,
                                         T steadyTime) {
    return ScalarDriver(trapezoid, constantValue, amplitude, -1, -1, -1, -1,
                        timeStart, -1, edgeTime, steadyTime);
  }

  /**
   * @brief Get the Gaussian Impulse Driver object
   *
   * @param constantValue
   * @param amplitude
   * @param t0 center of the pulse
   * @param sigma sigma of the gaussian
   * @return ScalarDriver
   */
  static ScalarDriver getGaussianImpulseDriver(T constantValue, T amplitude,
                                               T t0, T sigma) {
    return ScalarDriver(gaussimpulse, constantValue, amplitude, -1, -1, -1, -1,
                        t0, -1, sigma);
  }

  /**
   * @brief Get the Gaussian Impulse Driver object
   *
   * @param constantValue
   * @param amplitude
   * @param t0 center of the growth
   * @param sigma sigma of the gaussian
   * @return ScalarDriver
   */
  static ScalarDriver getGaussianStepDriver(T constantValue, T amplitude, T t0,
                                            T sigma) {
    return ScalarDriver(gaussimpulse, constantValue, amplitude, -1, -1, -1, -1,
                        t0, -1, sigma);
  }

  static ScalarDriver getCustomDriver(pybind11::function callback) {
    if (!callback) {
      throw std::runtime_error("Callback function is not set");
    }
    // check if the callback is callable has one argument
    // Check if the callback function has exactly one argument
    // Using cast to int to avoid type mismatch error
    if (pybind11::cast<int>(callback.attr("__code__").attr("co_argcount")) !=
        1) {
      throw std::runtime_error("Callback function must have one argument");
    }

    return ScalarDriver(custom, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, callback);
  }

  T getCurrentScalarValue(T &time) override {
    T returnValue = this->constantValue;
    if (this->update == pulse) {
      returnValue +=
          pulseTrain(this->amplitude, time, this->period, this->cycle);
    } else if (this->update == sine) {
      returnValue += this->amplitude *
                     sin(2 * M_PI * time * this->frequency + this->phase);
    } else if (this->update == posine) {
      returnValue += abs(this->amplitude *
                         sin(2 * M_PI * time * this->frequency + this->phase));
    } else if (this->update == halfsine) {
      const T tamp = this->amplitude *
                     sin(2 * M_PI * time * this->frequency + this->phase);
      if (tamp <= 0) {
        returnValue += tamp; // ? tamp >= 0. : 0.;
      }
    } else if (this->update == step) {
      returnValue +=
          stepUpdate(this->amplitude, time, this->timeStart, this->timeStop);
    } else if (this->update == trapezoid) {
      returnValue += trapezoidalUpdate(this->amplitude, time, this->timeStart,
                                       this->edgeTime, this->steadyTime);
    } else if (this->update == gaussimpulse) {
      const T gaussImp = this->amplitude * exp(-pow(time - this->timeStart, 2) /
                                               (2 * pow(this->edgeTime, 2)));
      returnValue += gaussImp;
    } else if (this->update == gaussstep) {
      const T gaussStep =
          0.5 * this->amplitude *
          (1 + std::erf((time - this->timeStart) / (sqrt(2) * this->edgeTime)));
      returnValue += gaussStep;
    } else if (this->update == custom) {
      // If it is, call the Python function
      pybind11::gil_scoped_acquire gil;
      try {
        return pybind11::cast<double>(m_callback(time));
      } catch (pybind11::error_already_set &e) {
        std::cerr << "Error in Python callback: " << e.what() << std::endl;
        throw std::runtime_error("Error in Python callback");
      }
    }
    return returnValue;
  }

  CVector<T> getUnitAxis() {
    return CVector<T>(1 ? this->constantValue : 0, 1 ? this->constantValue : 0,
                      1 ? this->constantValue : 0);
  }

  // override multiplication operator
  ScalarDriver<T> operator*(const T &val) {
    return ScalarDriver<T>(this->update, this->constantValue * val,
                           this->amplitude * val, this->frequency, this->phase,
                           this->period, this->cycle, this->timeStart,
                           this->timeStop, this->edgeTime, this->steadyTime);
  }

  ScalarDriver<T> operator*(const T &val) const {
    return ScalarDriver<T>(this->update, this->constantValue * val,
                           this->amplitude * val, this->frequency, this->phase,
                           this->period, this->cycle, this->timeStart,
                           this->timeStop, this->edgeTime, this->steadyTime);
  }

  // override *= operator
  ScalarDriver<T> operator*=(const T &val) {
    this->constantValue *= val;
    this->amplitude *= val;
    return *this;
  }

  // override addition operator
  ScalarDriver<T> operator+(const T &val) {
    return ScalarDriver<T>(this->update, this->constantValue + val,
                           this->amplitude + val, this->frequency, this->phase,
                           this->period, this->cycle, this->timeStart,
                           this->timeStop, this->edgeTime, this->steadyTime);
  }

  ScalarDriver operator+(const T &v) const {
    // Use non-const operator+ here
    return ScalarDriver<T>(this->update, this->constantValue + v,
                           this->amplitude + v, this->frequency, this->phase,
                           this->period, this->cycle, this->timeStart,
                           this->timeStop, this->edgeTime, this->steadyTime);
  };
  ScalarDriver<T> operator+=(const T &val) {
    this->constantValue += val;
    this->amplitude += val;
    return *this;
  }
};

template <typename T> class NullDriver : public ScalarDriver<T> {
public:
  NullDriver() = default;
  T getCurrentScalarValue(T &time) override { return 0.0; }
};

template <typename T> class AxialDriver : public Driver<T> {
private:
  std::vector<ScalarDriver<T>> drivers;

public:
  static AxialDriver getVectorAxialDriver(T x, T y, T z) {
    return AxialDriver(CVector<T>(x, y, z));
  }

  void applyMask(const std::vector<unsigned int> &mask) {
    assert(mask.size() == 3);
    for (int i = 0; i < 3; i++) {
      if (mask[i] == 0) {
        // Mask asks to nullify the driver
        this->drivers[i] = NullDriver<T>();
      } else if (mask[i] != 1) {
        throw std::runtime_error("Invalid mask value, mask must be binary!");
      }
    }
  }

  void applyMask(const CVector<T> &mask) {
    this->applyMask(std::vector<unsigned int>{(unsigned int)(mask[0]),
                                              (unsigned int)(mask[1]),
                                              (unsigned int)(mask[2])});
  }

  AxialDriver() {
    this->drivers = {NullDriver<T>(), NullDriver<T>(), NullDriver<T>()};
  }

  AxialDriver(ScalarDriver<T> x, ScalarDriver<T> y, ScalarDriver<T> z) {
    this->drivers = {x, y, z};
  }

  explicit AxialDriver(const CVector<T> &xyz)
      : AxialDriver(ScalarDriver<T>::getConstantDriver(xyz.x),
                    ScalarDriver<T>::getConstantDriver(xyz.y),
                    ScalarDriver<T>::getConstantDriver(xyz.z)) {}

  explicit AxialDriver(const T x, const T y, const T z)
      : AxialDriver(ScalarDriver<T>::getConstantDriver(x),
                    ScalarDriver<T>::getConstantDriver(y),
                    ScalarDriver<T>::getConstantDriver(z)) {}

  explicit AxialDriver(std::vector<ScalarDriver<T>> axialDrivers) {
    if (axialDrivers.size() != 3) {
      throw std::runtime_error("The axial driver can only have 3 axes!");
    }
    this->drivers = std::move(axialDrivers);
  }

  static AxialDriver getUniAxialDriver(const ScalarDriver<T> &in, Axis axis) {
    switch (axis) {
    case xaxis:
      return AxialDriver(in, NullDriver<T>(), NullDriver<T>());
    case yaxis:
      return AxialDriver(NullDriver<T>(), in, NullDriver<T>());
    case zaxis:
      return AxialDriver(NullDriver<T>(), NullDriver<T>(), in);
    case all:
      return AxialDriver(in, in, in);
    case none:
      return AxialDriver(NullDriver<T>(), NullDriver<T>(), NullDriver<T>());
    }
    return AxialDriver(NullDriver<T>(), NullDriver<T>(), NullDriver<T>());
  }
  CVector<T> getCurrentAxialDrivers(T time) {
    return CVector<T>(this->drivers[0].getCurrentScalarValue(time),
                      this->drivers[1].getCurrentScalarValue(time),
                      this->drivers[2].getCurrentScalarValue(time));
  }

  CVector<T> getConstantValues() {
    return CVector<T>(this->drivers[0].constantValue,
                      this->drivers[1].constantValue,
                      this->drivers[2].constantValue);
  }

  /**
   * Returns the mask for the Axial Driver.
   * For instance: a vector (1213, 123, 0) returns (1, 1, 0)
   * Note: This is not normalised
   * @return CVector<T>: mask for the driver
   */
  CVector<T> getUnitAxis() {
    return CVector<T>(this->drivers[0].constantValue != 0.0
                          ? this->drivers[0].constantValue /
                                std::abs(this->drivers[0].constantValue)
                          : 0.0,
                      this->drivers[1].constantValue != 0.0
                          ? this->drivers[1].constantValue /
                                std::abs(this->drivers[1].constantValue)
                          : 0.0,
                      this->drivers[2].constantValue != 0.0
                          ? this->drivers[2].constantValue /
                                std::abs(this->drivers[2].constantValue)
                          : 0.0);
  }
};

template <typename T> class NullAxialDriver : public AxialDriver<T> {
public:
  NullAxialDriver() = default;
};

#endif
