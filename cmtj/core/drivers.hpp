#ifndef DRIVERS_H
#define DRIVERS_H

#include "cvector.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

enum UpdateType
{
    constant,
    pulse,
    sine,
    step
};

template <typename T>
class Driver
{
public:
    // if the user wants to update, let them do that
    T constantValue, amplitude, frequency, phase,
        cycle, period, timeStart, timeStop;
    UpdateType update;
    Driver()
    {
        this->constantValue = 0.0;
        this->amplitude = 0.0;
        this->frequency = 0.0;
        this->phase = 0.0;
        this->period = 0.0;
        this->timeStart = 0.0;
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
};

template <typename T>
class ScalarDriver : public Driver<T>
{
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
        T timeStop = -1)
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
        if (update == pulse && ((period == -1) || (cycle == -1)))
        {
            throw std::runtime_error("Selected pulse train driver type but either period or cycle were not set");
        }
        else if (update == sine && (frequency == -1))
        {
            throw std::runtime_error("Selected sine driver type but frequency was not set");
        }
    }

    static ScalarDriver getConstantDriver(T constantValue)
    {
        return ScalarDriver(
            constant,
            constantValue);
    }

    static ScalarDriver getPulseDriver(T constantValue, T amplitude, T period, T cycle)
    {
        return ScalarDriver(
            pulse,
            constantValue,
            amplitude,
            -1, -1, period, cycle);
    }

    static ScalarDriver getSineDriver(T constantValue, T amplitude, T frequency, T phase)
    {
        return ScalarDriver(
            sine,
            constantValue,
            amplitude,
            frequency, phase);
    }

    static ScalarDriver getStepDriver(T constantValue, T amplitude, T timeStart, T timeStop)
    {
        return ScalarDriver(
            step,
            constantValue,
            amplitude,
            -1, -1, -1, -1, timeStart, timeStop);
    }

    T getCurrentScalarValue(T &time)
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
        else if (this->update == step)
        {
            returnValue += stepUpdate(this->amplitude, time, this->timeStart, this->timeStop);
        }

        return returnValue;
    }
};

template <typename T>
class NullDriver : public ScalarDriver<T>
{
public:
    NullDriver() = default;
    T getCurrentScalarValue(T time)
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
        this->applyMask((std::vector<unsigned int>){(unsigned int)(mask[0]),
                                                    (unsigned int)(mask[1]),
                                                    (unsigned int)(mask[2])});
    }

    AxialDriver()
    {
        this->drivers = {
            NullDriver<T>(),
            NullDriver<T>(),
            NullDriver<T>()};
    }

    AxialDriver(ScalarDriver<T> x,
                ScalarDriver<T> y,
                ScalarDriver<T> z)
    {
        this->drivers = {x, y, z};
    }

    AxialDriver(CVector<T> xyz) : AxialDriver(
                                      ScalarDriver<T>::getConstantDriver(xyz.x),
                                      ScalarDriver<T>::getConstantDriver(xyz.y),
                                      ScalarDriver<T>::getConstantDriver(xyz.z))
    {
    }

    AxialDriver(std::vector<ScalarDriver<T>> axialDrivers)
    {
        if (axialDrivers.size() != 3)
        {
            throw std::runtime_error("The axial driver can only have 3 axes!");
        }
        this->drivers = std::move(axialDrivers);
    }

    static AxialDriver getUniAxialDriver(ScalarDriver<T> in, Axis axis)
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
    CVector<T> getCurrentAxialDrivers(T time)
    {
        return CVector<T>(0., 0., 0.);
    }
    CVector<T> getConstantValues()
    {
        return CVector<T>(0., 0., 0.);
    }
};

#endif