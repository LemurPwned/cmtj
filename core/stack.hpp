#ifndef CORE_STACK_HPP_
#define CORE_STACK_HPP_

#include <algorithm>               // for for_each
#include <stddef.h>                // for size_t
#include <iostream>                // for operator<<, ofstream, string, basi...
#include <numeric>                 // for accumulate
#include <stdexcept>               // for runtime_error
#include <string>                  // for operator+, operator==, to_string
#include <unordered_map>           // for unordered_map
#include <vector>                  // for vector
#include "cvector.hpp"             // for CVector
#include "drivers.hpp"             // for ScalarDriver, AxialDriver, NullDriver
#include "junction.hpp"

template <typename T>
class Stack
{
    friend class Junction<T>;

private:
    ScalarDriver<T> currentDriver;
    std::unordered_map<std::string, std::vector<T>> stackLog;
    bool currentDriverSet = false;

protected:
    std::string topId, bottomId; // Ids of the top and bottom junctions
    T couplingStrength = 0;
    bool delayed = false;
    virtual T calculateStackResistance(std::vector<T> resistances) = 0;
    virtual T computeCouplingCurrentDensity(T currentDensity,
        CVector<T> m1, CVector<T> m2, CVector<T> p) = 0;

public:
    std::vector<Junction<T>> junctionList;

    void setDelayed(bool delay)
    {
        this->delayed = delay;
    }

    void setMagnetisation(unsigned int junctionId, const std::string& layerId, CVector<T> mag)
    {
        this->junctionList[junctionId].setLayerMagnetisation(layerId, mag);
    }

    const CVector<T> getMagnetisation(unsigned int junctionId, const std::string& layerId)
    {
        return this->junctionList[junctionId].getLayerMagnetisation(layerId);
    }

    void setOerstedFieldDriver(const AxialDriver<T>& oDriver)
    {
        for (auto& j : this->junctionList)
        {
            j.setLayerOerstedFieldDriver("all", oDriver);
        }
    }

    void setExternalFieldDriver(const AxialDriver<T>& fDriver)
    {
        for (auto& j : this->junctionList)
        {
            j.setLayerExternalFieldDriver("all", fDriver);
        }
    }

    void resetCoupledCurrentDriver()
    {
        this->currentDriver = NullDriver<T>();
        for (auto& j : this->junctionList)
        {
            j.setLayerCurrentDriver("all", this->currentDriver);
        }
        this->currentDriverSet = false;
    }

    void setCoupledCurrentDriver(const ScalarDriver<T>& cDriver)
    {
        this->currentDriver = cDriver;
        for (auto& j : this->junctionList)
        {
            j.setLayerCurrentDriver("all", this->currentDriver);
        }
        this->currentDriverSet = true;
    }

    Stack(std::vector<Junction<T>> inputStack,
        const std::string& topId,
        const std::string& bottomId) : topId(topId), bottomId(bottomId)
    {
        if (inputStack.size() < 2)
        {
            throw std::runtime_error("Stack must have at least 2 junctions!");
        }
        this->junctionList = std::move(inputStack);
        if (std::any_of(this->junctionList.begin(),
            this->junctionList.end(),
            [](const Junction<T>& j) { return j.MR_mode != Junction<T>::MRmode::CLASSIC; }))
        {
            throw std::runtime_error("Junction has a non-classic magnetoresitance mode!");
        }
    }
    void
        saveLogs(std::string fileSave)
    {
        if (fileSave == "")
        {
            // if there's an empty fn, don't save
            std::cout << "Ignoring file save to an empty filename" << std::endl;
            return;
        }
        std::ofstream logFile;
        logFile.open(fileSave);
        for (const auto& keyPair : this->stackLog)
        {
            logFile << keyPair.first << ";";
        }
        logFile << "\n";
        for (unsigned int i = 0; i < this->stackLog["time"].size(); i++)
        {
            for (const auto& keyPair : this->stackLog)
            {
                logFile << keyPair.second[i] << ";";
            }
            logFile << "\n";
        }
        logFile.close();
    }

    void setCouplingStrength(T coupling)
    {
        this->couplingStrength = coupling;
    }

    void logStackData(T t, T resistance, std::vector<T> timeCurrents)
    {
        this->stackLog["Resistance"].push_back(resistance);
        for (std::size_t j = 0; j < timeCurrents.size(); ++j)
        {
            this->stackLog["I_" + std::to_string(j)].push_back(timeCurrents[j]);
        }
        this->stackLog["time"].push_back(t);
    }

    void clearLogs()
    {
        for (auto& j : this->junctionList)
        {
            j.clearLog();
        }
        this->stackLog.clear();
    }

    std::unordered_map<std::string, std::vector<T>>& getLog()
    {
        return this->stackLog;
    }
    std::unordered_map<std::string, std::vector<T>>& getLog(unsigned int id)
    {
        if (id <= this->junctionList.size())
        {
            return this->junctionList[id].getLog();
        }
        throw std::runtime_error("Asking for id of a non-existing junction!");
    }

    const CVector<T> getPolarisationVector()
    {
        CVector<T> probe = junctionList[0].getLayer(this->topId).referenceLayer;
        for (std::size_t i = 1; i < junctionList.size(); ++i)
        {
            if (probe != junctionList[i].getLayer(this->topId).referenceLayer)
                throw std::runtime_error("Polarisation vectors are not equal in stack");
        }

        if (!probe.length())
        {
            throw std::runtime_error("Polarisation is not set!");
        }
        return probe;
    }

    const bool isTwoLayerMemberStack() {
        for (const auto& j : this->junctionList)
        {
            if (j.layerNo >= 3)
            {
                throw std::runtime_error("At least one junction has more than 2 layers!"
                    " It's not supported now.");
            }
            else if (j.layerNo != 2)
            {
                return false;
            }
        }
        return true;
    }

    void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11)
    {
        const unsigned int writeEvery = (int)(writeFrequency / timeStep);
        const unsigned int totalIterations = (int)(totalTime / timeStep);

        if (timeStep > writeFrequency)
        {
            std::runtime_error("The time step cannot be larger than write frequency!");
        }

        // pick a solver based on drivers
        std::vector<SolverMode> modes;
        auto localRunner = &Junction<T>::runMultiLayerSolver;
        for (auto& j : this->junctionList)
        {
            auto [runner, solver, mode] = j.getSolver(RK4, totalIterations);
            modes.push_back(mode);
            localRunner = runner;
            // TODO: handle the rare case when the user mixes 1 layer with 2 layer junction
            // in the same stack -- i.e. runner is runSingleLayerSolver and runMultiLayerSolver
        }
        auto solver = &Layer<T>::rk4_step;
        if (!std::equal(modes.begin() + 1, modes.end(), modes.begin()))
        {
            throw std::runtime_error("Junctions have different solver modes!"
                " Set the same solver mode for all junctions explicitly."
                " Do not mix stochastic and deterministic solvers!");
        }


        T tCurrent;
        std::vector<T> timeResistances(junctionList.size());
        std::vector<T> timeCurrents(junctionList.size());
        std::vector<CVector<T>> frozenMags(junctionList.size());
        std::vector<CVector<T>> frozenPols(junctionList.size());

        const bool isTwoLayerStack = this->isTwoLayerMemberStack();
        for (unsigned int i = 0; i < totalIterations; i++)
        {
            T t = i * timeStep;

            // stash the magnetisations first
            for (std::size_t j = 0; j < junctionList.size(); ++j) {
                frozenMags[j] = junctionList[j].getLayerMagnetisation(this->topId);
                if (isTwoLayerStack) {
                    frozenPols[j] = junctionList[j].getLayerMagnetisation(this->bottomId);
                }
                else {
                    frozenPols[j] = this->getPolarisationVector();
                }
            }
            const T plainCurrent = this->currentDriver.getCurrentScalarValue(t);
            T coupledCurrent = plainCurrent;
            for (std::size_t j = 0; j < junctionList.size(); ++j)
            {
                // skip first junction
                // modify the standing layer constant current
                if (j > 0) {
                    if (this->delayed) {
                        // accumulate coupling
                        coupledCurrent = coupledCurrent + this->computeCouplingCurrentDensity(
                            // j -> k, j-1 -> k'
                            coupledCurrent, frozenMags[j], frozenMags[j - 1], frozenPols[j]);
                    }
                    else {
                        coupledCurrent = coupledCurrent + this->computeCouplingCurrentDensity(
                            // j -> k, j-1 -> k'
                            coupledCurrent, junctionList[j].getLayerMagnetisation(this->topId),
                            junctionList[j - 1].getLayerMagnetisation(this->topId),
                            junctionList[j].getLayerMagnetisation(this->bottomId));
                    }
                    tCurrent = coupledCurrent;
                }
                else {
                    tCurrent = plainCurrent;
                }

                (junctionList[j].*localRunner)(solver, t, timeStep);
                // change the instant value of the current before the
                // the resistance is calculated
                // compute the next j+1 input to the current.
                const auto resistance = junctionList[j].getMagnetoresistance();
                timeResistances[j] = resistance[0];
                timeCurrents[j] = tCurrent;
            }
            if (!(i % writeEvery))
            {
                const T magRes = this->calculateStackResistance(timeResistances);
                this->logStackData(t, magRes, timeCurrents);
                for (auto& jun : this->junctionList)
                    jun.logLayerParams(t, timeStep, false);
            }
        }
    }
};
template <typename T>
class SeriesStack : public Stack<T>
{
    T calculateStackResistance(std::vector<T> resistances) override
    {
        const T resSum = std::accumulate(resistances.begin(),
            resistances.end(),
            0.0);
        return resSum;
    }

    T computeCouplingCurrentDensity(T currentDensity, CVector<T> m1, CVector<T> m2, CVector<T> p) override
    {
        const T m1Comp = c_dot(m1, p);
        const T m2Comp = c_dot(m2, p);
        const T coupledI = currentDensity * this->couplingStrength * (m1Comp + m2Comp);
        return coupledI;
    }

public:
    explicit SeriesStack(const std::vector<Junction<T>>& jL,
        const std::string& topId = "free",
        const std::string& bottomId = "bottom") : Stack<T>(jL, topId, bottomId) {}
};
template <typename T>
class ParallelStack : public Stack<T>
{
    T calculateStackResistance(std::vector<T> resistances) override
    {
        T invSum = 0.0;
        std::for_each(resistances.begin(), resistances.end(), [&](T res)
            { invSum += 1.0 / res; });
        return 1. / invSum;
    }

    T computeCouplingCurrentDensity(T currentDensity, CVector<T> m1, CVector<T> m2, CVector<T> p) override
    {
        const T m1Comp = c_dot(m1, p);
        const T m2Comp = c_dot(m2, p);
        const T coupledI = currentDensity * this->couplingStrength * (m1Comp - m2Comp);
        return coupledI;
    }

public:
    explicit ParallelStack(const std::vector<Junction<T>>& jL,
        const std::string& topId = "free",
        const std::string& bottomId = "bottom") : Stack<T>(jL, topId, bottomId) {}
};
#endif // CORE_STACK_HPP_
