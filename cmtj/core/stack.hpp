#ifndef STACK_H
#define STACK_H

#include "junction.hpp"
#include "parallel.hpp"

class Stack
{
    friend class Junction;

protected:
    std::vector<Junction> junctionList;
    std::map<std::string, std::vector<double>> stackLog;
    virtual double calculateStackResistance(std::vector<double> resistances) = 0;

public:
    Stack(std::vector<Junction> junctionList)
    {
        this->junctionList = std::move(junctionList);
        for (auto &j : this->junctionList)
        {
            if (j.MR_mode != Junction::CLASSIC)
            {
                throw std::runtime_error("Junction has a non-classic magnetoresitance mode!");
            }
        }
    }

    void logStackData(double t, double resistance)
    {
        this->stackLog["Resistance"].push_back(resistance);
        this->stackLog["time"].push_back(t);
    }

    void dumpStackData(std::string filename)
    {
        ComputeUtil::customResultMap(this->stackLog, filename);
    }

    void runSimulation(double totalTime, double timeStep = 1e-13, double writeFrequency = 1e-11)
    {
        const unsigned int writeEvery = (int)(writeFrequency / timeStep) - 1;
        const unsigned int totalIterations = (int)(totalTime / timeStep);
        double t;

        std::vector<double> timeResistances(junctionList.size());
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (unsigned int i = 0; i < totalIterations; i++)
        {
            t = i * timeStep;
            for (Junction &junction : junctionList)
            {
                junction.runSingleRK4Iteration(t, timeStep);
                const auto resistance = junction.getMagnetoresistance();
                timeResistances.push_back(resistance[0]);
            }
            if (!(i % writeEvery))
            {
                const double magRes = calculateStackResistance(timeResistances);
                this->logStackData(t, magRes);
            }
            timeResistances.clear();
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    }
};
class SeriesStack : public Stack
{
    double calculateStackResistance(std::vector<double> resistances)
    {
        const double resSum = std::accumulate(resistances.begin(),
                                              resistances.end(),
                                              0.0);
        return resSum;
    }

public:
    SeriesStack(std::vector<Junction> jL) : Stack(jL) {}
};
class ParallelStack : public Stack
{
    double calculateStackResistance(std::vector<double> resistances)
    {
        double invSum = 0.0;
        std::for_each(resistances.begin(), resistances.end(), [&](double res) {
            invSum += 1.0 / res;
        });
        return 1 / invSum;
    }

public:
    ParallelStack(std::vector<Junction> jL) : Stack(jL) {}
};
#endif