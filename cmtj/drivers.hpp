#ifndef DRIVERS_H
#define DRIVERS_H

#include "cvector.hpp"
#include <iostream>
#include <cmath>

enum UpdateType
{
    constant,
    pulse,
    sine
};

class Driver
{
public:
    // if the user wants to update, let them do that
    double constantValue, amplitude, frequency, phase, cycle, period;
    UpdateType update;
    Driver(double constantValue,
           double amplitude,
           double frequency,
           double phase,
           double period,
           double cycle,
           UpdateType update) : constantValue(constantValue),
                                amplitude(amplitude),
                                frequency(frequency),
                                phase(phase),
                                period(period),
                                cycle(cycle),
                                update(update) {}
};

class AxialDriver
{
private:
    std::vector<ScalarDriver> drivers;

public:
    AxialDriver(std::vector<ScalarDriver> axialDrivers)
    {
        this->drivers = std::move(axialDrivers);
    }

    CVector getCurrentAxialDrivers(double time)
    {
        return CVector(
            this->drivers[0].getCurrentScalarValue(time),
            this->drivers[1].getCurrentScalarValue(time),
            this->drivers[2].getCurrentScalarValue(time));
    }
};

class ScalarDriver : public Driver
{

public:
    ScalarDriver(
        double constantValue = 0,
        double amplitude = 0,
        double frequency = -1,
        double phase = 0,
        double period = -1,
        double cycle = -1,
        UpdateType update = constant) : Driver(constantValue,
                                               amplitude,
                                               frequency,
                                               phase,
                                               period,
                                               cycle,
                                               update)
    {
        if (update == pulse && (period == -1 || cycle == -1))
        {
            std::runtime_error("Selected pulse train driver type but either period or cycle were not set");
        }
        else if (update == sine && (frequency == -1))
        {
            std::runtime_error("Selected sine driver type but frequency was not set");
        }
    }

    double pulseTrain(double amplitude, double time, double T, double cycle)
    {
        const int n = (int)(time / T);
        const double dT = cycle * T;
        const double nT = n * T;
        if (nT <= time <= (nT + dT))
        {
            return amplitude;
        }
        else
        {
            return 0;
        }
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

        return returnValue;
    }
};

class NullDriver : public ScalarDriver
{
    double getCurrentScalarValue(double time)
    {
        return 0.0;
    }
};

#endif