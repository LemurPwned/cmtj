#ifndef DRIVERS_H
#define DRIVERS_H

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <stdlib.h>     // for abs
#include <cassert>      // for assert
#define _USE_MATH_DEFINES
#include <cmath>        // for M_PI
#include <stdexcept>    // for runtime_error
#include <vector>       // for vector
#include "cvector.hpp"  // for CVector

enum UpdateType
{
    constant,
    pulse,
    sine,
    step,
    posine,
    halfsine,
    trapezoid
};

template <typename T>
class Driver
{
protected:
    // if the user wants to update, let them do that
    T constantValue, amplitude, frequency, phase,
        period, cycle, timeStart, timeStop;
    UpdateType update;
public:
    Driver()
    {
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
    Driver(UpdateType update,
        T constantValue,
        T amplitude,
        T frequency,
        T phase,
        T period,
        T cycle,
        T timeStart,
        T timeStop) : constantValue(constantValue),
        amplitude(amplitude),
        frequency(frequency),
        phase(phase),
        period(period),
        cycle(cycle),
        timeStart(timeStart),
        timeStop(timeStop),
        update(update)

    {
    }
    virtual T getCurrentScalarValue(T& time)
    {
        return 0;
    };
    virtual ~Driver() = default;
};

template <typename T>
class ScalarDriver : public Driver<T>
{
private:
    T edgeTime = 0;
    T steadyTime = 0;
protected:
    T stepUpdate(T amplitude, T time, T timeStart, T timeStop)
    {
        if (time >= timeStart && time <= timeStop)
        {
            return amplitude;
        }
        else
        {
            return 0.0;
        }
    }
    T pulseTrain(T amplitude, T time, T period, T cycle)
    {
        const int n = (int)(time / period);
        const T dT = cycle * period;
        const T nT = n * period;
        if (nT <= time && time <= (nT + dT))
        {
            return amplitude;
        }
        else
        {
            return 0;
        }
    }

    T trapezoidalUpdate(T amplitude, T time, T timeStart, T edgeTime, T steadyTime) {
        if (time < timeStart) {
            return 0;
        }
        // growth
        else if (time <= timeStart + edgeTime) {
            return  (amplitude / edgeTime) * (time - timeStart);
        }
        // steady
        else if (time <= timeStart + edgeTime + steadyTime) {
            return amplitude;
        }
        // decay
        else if (time <= timeStart + 2 * edgeTime + steadyTime) {
            return  amplitude - (amplitude / edgeTime) * (time - (timeStart + edgeTime + steadyTime));
        }
        return 0;
    }


public:
    ScalarDriver(
        UpdateType update = constant,
        T constantValue = 0,
        T amplitude = 0,
        T frequency = -1,
        T phase = 0,
        T period = -1,
        T cycle = -1,
        T timeStart = -1,
        T timeStop = -1,
        T edgeTime = -1,
        T steadyTime = -1)
        : Driver<T>(update,
            constantValue,
            amplitude,
            frequency,
            phase,
            period,
            cycle,
            timeStart,
            timeStop)
    {
        this->edgeTime = edgeTime;
        this->steadyTime = steadyTime;
        if (update == pulse && ((period == -1) || (cycle == -1)))
        {
            throw std::runtime_error("Selected pulse train driver type but either period or cycle were not set");
        }
        else if (update == sine && (frequency == -1))
        {
            throw std::runtime_error("Selected sine driver type but frequency was not set");
        }
    }

    /**
     * Constant driver produces a constant signal of a fixed amplitude.
     * @param constantValue: constant value of the driver (constant offset/amplitude)
     */
    static ScalarDriver getConstantDriver(T constantValue)
    {
        return ScalarDriver(
            constant,
            constantValue);
    }

    /**
     * Produces a square pulse of certain period and cycle
     * @param constantValue: offset (vertical) of the pulse. The pulse amplitude will be added to this.
     * @param amplitude: amplitude of the pulse signal
     * @param period: period of the signal in seconds
     * @param cycle: duty cycle of the signal -- a fraction between [0 and 1].
     */
    static ScalarDriver getPulseDriver(T constantValue, T amplitude, T period, T cycle)
    {
        return ScalarDriver(
            pulse,
            constantValue,
            amplitude,
            -1, -1, period, cycle);
    }

    /**
     * Produces a sinusoidal signal with some offset (constantValue), amplitude frequency and phase offset.
     * @param constantValue: vertical offset. The sine will oscillate around this value.
     * @param amplitude: amplitude of the sine wave
     * @param frequency: frequency of the sine
     * @param phase: phase of the sine in radians.
     */
    static ScalarDriver getSineDriver(T constantValue, T amplitude, T frequency, T phase)
    {
        return ScalarDriver(
            sine,
            constantValue,
            amplitude,
            frequency, phase);
    }

    /**
     * Produces a positive sine signal with some offset (constantValue), amplitude frequency and phase offset.
     * @param constantValue: vertical offset. The sine will oscillate around this value.
     * @param amplitude: amplitude of the sine wave
     * @param frequency: frequency of the sine
     * @param phase: phase of the sine in radians.
     */
    static ScalarDriver getPosSineDriver(T constantValue, T amplitude, T frequency, T phase)
    {
        return ScalarDriver(
            posine,
            constantValue,
            amplitude,
            frequency, phase);
    }

    static ScalarDriver getHalfSineDriver(T constantValue, T amplitude, T frequency, T phase)
    {
        return ScalarDriver(
            halfsine,
            constantValue,
            amplitude,
            frequency, phase);
    }
    /**
     * Get a step driver. It has amplitude between timeStart and timeStop and 0 elsewhere
     * @param constantValue: offset of the pulse (vertical)
     * @param amplitude: amplitude that is added on top of the constantValue
     * @param timeStart: start of the pulse
     * @param timeStop: when the pulse ends
     */
    static ScalarDriver getStepDriver(T constantValue, T amplitude, T timeStart, T timeStop)
    {
        if (timeStop <= timeStart)
        {
            throw std::runtime_error("Start time cannot be later than stop time!");
        }
        return ScalarDriver(
            step,
            constantValue,
            amplitude,
            -1, -1, -1, -1, timeStart, timeStop);
    }

    /**
     * Get a trapezoidal driver. It has amplitude between timeStart and timeStop and 0 elsewhere
     * @param constantValue: offset of the pulse (vertical)
     * @param amplitude: amplitude that is added on top of the constantValue
     * @param timeStart: start of the pulse
     * @param edgeTime: time it takes to reach the maximum amplitude
     * @param steadyTime: time it spends in a steady state
     */
    static ScalarDriver getTrapezoidDriver(T constantValue, T amplitude, T timeStart, T edgeTime, T steadyTime) {
        return ScalarDriver(
            trapezoid,
            constantValue,
            amplitude,
            -1, -1, -1, -1, timeStart, -1, edgeTime, steadyTime);
    }

    T getCurrentScalarValue(T& time) override
    {
        T returnValue = this->constantValue;
        if (this->update == pulse)
        {
            returnValue += pulseTrain(this->amplitude, time, this->period, this->cycle);
        }
        else if (this->update == sine)
        {
            returnValue += this->amplitude * sin(2 * M_PI * time * this->frequency + this->phase);
        }
        else if (this->update == posine)
        {
            returnValue += abs(this->amplitude * sin(2 * M_PI * time * this->frequency + this->phase));
        }
        else if (this->update == halfsine)
        {
            const T tamp = this->amplitude * sin(2 * M_PI * time * this->frequency + this->phase);
            if (tamp <= 0)
            {
                returnValue += tamp; // ? tamp >= 0. : 0.;
            }
        }
        else if (this->update == step)
        {
            returnValue += stepUpdate(this->amplitude, time, this->timeStart, this->timeStop);
        }
        else if (this->update == trapezoid) {
            returnValue += trapezoidalUpdate(this->amplitude, time, this->timeStart, this->edgeTime, this->steadyTime);
        }

        return returnValue;
    }
    void setConstantValue(const T& val)
    {
        this->constantValue = val;
    }
};



template <typename T>
class NullDriver : public ScalarDriver<T>
{
public:
    NullDriver() = default;
    T getCurrentScalarValue(T& time) override
    {
        return 0.0;
    }
};

enum Axis
{
    xaxis,
    yaxis,
    zaxis
};

template <typename T>
class AxialDriver : public Driver<T>
{
private:
    std::vector<ScalarDriver<T>> drivers;

public:
    static AxialDriver getVectorAxialDriver(T x, T y, T z)
    {
        return AxialDriver(CVector<T>(x, y, z));
    }

    void applyMask(std::vector<unsigned int> mask)
    {
        assert(mask.size() == 3);
        for (int i = 0; i < 3; i++)
        {
            if (mask[i] == 0)
            {
                // Mask asks to nullify the driver
                this->drivers[i] = NullDriver<T>();
            }
            else if (mask[i] != 1)
            {
                throw std::runtime_error("Invalid mask value, mask must be binary!");
            }
        }
    }

    void applyMask(CVector<T> mask)
    {
        this->applyMask(std::vector<unsigned int>{(unsigned int)(mask[0]),
            (unsigned int)(mask[1]),
            (unsigned int)(mask[2])});
    }

    AxialDriver()
    {
        this->drivers = {
            NullDriver<T>(),
            NullDriver<T>(),
            NullDriver<T>() };
    }

    AxialDriver(ScalarDriver<T> x,
        ScalarDriver<T> y,
        ScalarDriver<T> z)
    {
        this->drivers = { x, y, z };
    }

    explicit AxialDriver(const CVector<T>& xyz) : AxialDriver(
        ScalarDriver<T>::getConstantDriver(xyz.x),
        ScalarDriver<T>::getConstantDriver(xyz.y),
        ScalarDriver<T>::getConstantDriver(xyz.z))
    {
    }

    explicit AxialDriver(std::vector<ScalarDriver<T>> axialDrivers)
    {
        if (axialDrivers.size() != 3)
        {
            throw std::runtime_error("The axial driver can only have 3 axes!");
        }
        this->drivers = std::move(axialDrivers);
    }

    static AxialDriver getUniAxialDriver(const ScalarDriver<T>& in, Axis axis)
    {
        switch (axis)
        {
        case xaxis:
            return AxialDriver(in, NullDriver<T>(), NullDriver<T>());
        case yaxis:
            return AxialDriver(NullDriver<T>(), in, NullDriver<T>());
        case zaxis:
            return AxialDriver(NullDriver<T>(), NullDriver<T>(), in);
        }
        return AxialDriver(NullDriver<T>(), NullDriver<T>(), NullDriver<T>());
    }
    CVector<T>
        getCurrentAxialDrivers(T time)
    {
        return CVector<T>(
            this->drivers[0].getCurrentScalarValue(time),
            this->drivers[1].getCurrentScalarValue(time),
            this->drivers[2].getCurrentScalarValue(time));
    }

    CVector<T> getConstantValues()
    {
        return CVector<T>(
            this->drivers[0].constantValue,
            this->drivers[1].constantValue,
            this->drivers[2].constantValue);
    }
    /**
     * Returns the mask for the Axial Driver.
     * For instance: a vector (1213, 123, 0) returns (1, 1, 0)
     * Note: This is not normalised
     * @return CVector<T>: mask for the driver
     */
    CVector<T> getUnitAxis()
    {
        return CVector<T>(
            this->drivers[0].constantValue != 0.0 ? this->drivers[0].constantValue / std::abs(this->drivers[0].constantValue) : 0.0,
            this->drivers[1].constantValue != 0.0 ? this->drivers[1].constantValue / std::abs(this->drivers[1].constantValue) : 0.0,
            this->drivers[2].constantValue != 0.0 ? this->drivers[2].constantValue / std::abs(this->drivers[2].constantValue) : 0.0);
    }
};

template <typename T>
class NullAxialDriver : public AxialDriver<T>
{
public:
    NullAxialDriver() = default;
    CVector<T> getCurrentAxialDrivers([[maybe_unused]] T time)
    {
        return CVector<T>(0., 0., 0.);
    }
    CVector<T> getConstantValues()
    {
        return CVector<T>(0., 0., 0.);
    }
};

#endif
