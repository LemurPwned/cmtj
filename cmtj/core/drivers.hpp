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

class Driver
{
public:
    // if the user wants to update, let them do that
    double constantValue, amplitude, frequency, phase,
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
           double constantValue,
           double amplitude,
           double frequency,
           double phase,
           double period,
           double cycle,
           double timeStart,
           double timeStop) : update(update), constantValue(constantValue),
                              amplitude(amplitude),
                              frequency(frequency),
                              phase(phase),
                              period(period),
                              cycle(cycle),
                              timeStart(timeStart),
                              timeStop(timeStop)
    {
    }
};

class ScalarDriver : public Driver
{
protected:
    double stepUpdate(double amplitude, double time, double timeStart, double timeStop)
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
    double pulseTrain(double amplitude, double time, double T, double cycle)
    {
        const int n = (int)(time / T);
        const double dT = cycle * T;
        const double nT = n * T;
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
        double constantValue = 0,
        double amplitude = 0,
        double frequency = -1,
        double phase = 0,
        double period = -1,
        double cycle = -1,
        double timeStart = -1,
        double timeStop = -1)
        : Driver(update, constantValue,
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

    static ScalarDriver getConstantDriver(double constantValue)
    {
        return ScalarDriver(
            constant,
            constantValue);
    }

    static ScalarDriver getPulseDriver(double constantValue, double amplitude, double period, double cycle)
    {
        return ScalarDriver(
            pulse,
            constantValue,
            amplitude,
            -1, -1, period, cycle);
    }

    static ScalarDriver getSineDriver(double constantValue, double amplitude, double frequency, double phase)
    {
        return ScalarDriver(
            sine,
            constantValue,
            amplitude,
            frequency, phase);
    }

    static ScalarDriver getStepDriver(double constantValue, double amplitude, double timeStart, double timeStop)
    {
        return ScalarDriver(
            step,
            constantValue,
            amplitude,
            -1, -1, -1, -1, timeStart, timeStop);
    }

    double getCurrentScalarValue(double time)
    {
        double returnValue = constantValue;
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

class NullDriver : public ScalarDriver
{
public:
    NullDriver() = default;
    double getCurrentScalarValue(double time)
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

class AxialDriver : public Driver
{
private:
    std::vector<ScalarDriver> drivers;

public:
    void applyMask(std::vector<unsigned int> mask)
    {
        assert(mask.size() == 3);
        for (int i = 0; i < 3; i++)
        {
            if (mask[i] == 0)
            {
                // Mask asks to nullify the driver
                this->drivers[i] = NullDriver();
            }
            else if (mask[i] != 1)
            {
                throw std::runtime_error("Invalid mask value, mask must be binary!");
            }
        }
    }

    void applyMask(CVector mask)
    {
        this->applyMask((std::vector<unsigned int>){(unsigned int)(mask[0]),
                                                    (unsigned int)(mask[1]),
                                                    (unsigned int)(mask[2])});
    }

    AxialDriver()
    {
        this->drivers = {
            NullDriver(),
            NullDriver(),
            NullDriver()};
    }

    AxialDriver(ScalarDriver x,
                ScalarDriver y,
                ScalarDriver z)
    {
        this->drivers = {x, y, z};
    }

    AxialDriver(CVector xyz) : AxialDriver(
                                   ScalarDriver::getConstantDriver(xyz.x),
                                   ScalarDriver::getConstantDriver(xyz.y),
                                   ScalarDriver::getConstantDriver(xyz.z))
    {
    }

    AxialDriver(std::vector<ScalarDriver> axialDrivers)
    {
        if (axialDrivers.size() != 3)
        {
            throw std::runtime_error("The axial driver can only have 3 axes!");
        }
        this->drivers = std::move(axialDrivers);
    }

    static AxialDriver getUniAxialDriver(ScalarDriver in, Axis axis)
    {
        switch (axis)
        {
        case xaxis:
            return AxialDriver(in, NullDriver(), NullDriver());
        case yaxis:
            return AxialDriver(NullDriver(), in, NullDriver());
        case zaxis:
            return AxialDriver(NullDriver(), NullDriver(), in);
        }
    }
    CVector
    getCurrentAxialDrivers(double time)
    {
        return CVector(
            this->drivers[0].getCurrentScalarValue(time),
            this->drivers[1].getCurrentScalarValue(time),
            this->drivers[2].getCurrentScalarValue(time));
    }
};

#endif