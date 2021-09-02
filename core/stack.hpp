#ifndef STACK_H
#define STACK_H

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
    T couplingStrength = 0;
    virtual T calculateStackResistance(std::vector<T> resistances) = 0;
    virtual T computeCouplingCurrentDensity(T currentDensity, CVector<T> m1, CVector<T> m2, CVector<T> p) = 0;

public:
    std::vector<Junction<T>> junctionList;

    void setMagnetisation(unsigned int junctionId, std::string layerId, CVector<T> mag)
    {
        this->junctionList[junctionId].setLayerMagnetisation(layerId, mag);
    }

    void setOerstedFieldDriver(AxialDriver<T> oDriver)
    {
        for (auto &j : this->junctionList)
        {
            j.setLayerOerstedFieldDriver("all", oDriver);
        }
    }

    void setExternalFieldDriver(AxialDriver<T> fDriver)
    {
        for (auto &j : this->junctionList)
        {
            j.setLayerExternalFieldDriver("all", fDriver);
        }
    }

    void resetCoupledCurrentDriver()
    {
        this->currentDriver = NullDriver<T>();
        for (auto &j : this->junctionList)
        {
            j.setLayerCurrentDriver("all", this->currentDriver);
        }
        this->currentDriverSet = false;
    }

    void setCoupledCurrentDriver(ScalarDriver<T> cDriver)
    {
        this->currentDriver = cDriver;
        for (auto &j : this->junctionList)
        {
            j.setLayerCurrentDriver("all", this->currentDriver);
        }
        this->currentDriverSet = true;
    }

    Stack(std::vector<Junction<T>> inputStack)
    {
        this->junctionList = std::move(inputStack);

        for (auto &j : this->junctionList)
        {
            if (j.MR_mode != Junction<T>::MRmode::CLASSIC)
            {
                throw std::runtime_error("Junction has a non-classic magnetoresitance mode!");
            }
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
        for (const auto &keyPair : this->stackLog)
        {
            logFile << keyPair.first << ";";
        }
        logFile << "\n";
        for (unsigned int i = 0; i < this->stackLog["time"].size(); i++)
        {
            for (const auto &keyPair : this->stackLog)
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
            this->stackLog["Current_" + std::to_string(j)].push_back(timeCurrents[j]);
        }
        this->stackLog["time"].push_back(t);
    }

    void clearLogs()
    {
        for (auto &j : this->junctionList)
        {
            j.clearLog();
        }
        this->stackLog.clear();
    }

    std::unordered_map<std::string, std::vector<T>> &getLog()
    {
        return this->stackLog;
    }
    std::unordered_map<std::string, std::vector<T>> &getLog(unsigned int id)
    {
        if (id <= this->junctionList.size())
        {
            return this->junctionList[id].getLog();
        }
        throw std::runtime_error("Asking for id of non-existing junction!");
    }

    const CVector<T> getPolarisationVector()
    {
        std::vector<CVector<T>> polarisation(junctionList.size());
        for (std::size_t i = 0; i < junctionList.size(); ++i)
            polarisation[i] = junctionList[i].getLayer("free").referenceLayer;

        if (!(std::adjacent_find(polarisation.begin(), polarisation.end(), std::not_equal_to<>()) == polarisation.end()))
        {
            throw std::runtime_error("Polarisation vectors are not equal in stack");
        }
        if (!polarisation[0].length())
        {
            throw std::runtime_error("Polarisation is not set!");
        }
        return polarisation[0];
    }

    void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11)
    {
        const unsigned int writeEvery = (int)(writeFrequency / timeStep);
        const unsigned int totalIterations = (int)(totalTime / timeStep);

        if (timeStep > writeFrequency)
        {
            std::runtime_error("The time step cannot be larger than write frequency!");
        }
        T t;

        // pick a solver based on drivers
        auto solv = &Layer<T>::rk4_step;
        for (auto &j : this->junctionList)
        {
            for (auto &l : j.layers)
            {
                if (l.hasTemperature())
                {
                    // if at least one temp. driver is set
                    // then use euler_heun for consistency
                    std::cout << "Warning: using Euler-Heun in stack computation" << std::endl;
                    solv = &Layer<T>::euler_heun;
                    goto labelEndLoop;
                }
            }
        }
    labelEndLoop:
        T coupledCurrent = 0;
        T tCurrent = 0;
        std::vector<T> timeResistances(junctionList.size());
        std::vector<T> timeCurrents(junctionList.size());
        std::vector<CVector<T>> frozenMags(junctionList.size());

        const CVector<T> pol = this->getPolarisationVector();

        for (unsigned int i = 0; i < totalIterations; i++)
        {

            t = i * timeStep;
            coupledCurrent = this->currentDriver.getCurrentScalarValue(t);

            // stash the magnetisations first
            for (std::size_t j = 0; j < junctionList.size(); ++j)
                frozenMags[j] = junctionList[j].getLayerMagnetisation("free");

            for (std::size_t j = 0; j < junctionList.size(); ++j)
            {
                // skip first junction
                // modify the standing layer constant current
                if (j > 0)
                    tCurrent = coupledCurrent + this->computeCouplingCurrentDensity(
                                                    coupledCurrent, frozenMags[j], frozenMags[j - 1], pol);
                else
                    tCurrent = coupledCurrent;

                junctionList[j].setLayerCurrentDriver("all", ScalarDriver<T>::getConstantDriver(
                                                                 tCurrent));

                // solve the equation
                if (this->junctionList[j].layerNo == 1)
                {
                    junctionList[j].runSingleLayerSolver(solv, t, timeStep);
                }
                else
                {
                    junctionList[j].runMultiLayerSolver(solv, t, timeStep);
                }

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
                for (auto &jun : this->junctionList)
                    jun.logLayerParams(t, false);
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
    SeriesStack(std::vector<Junction<T>> jL) : Stack<T>(jL) {}
};
template <typename T>
class ParallelStack : public Stack<T>
{
    T calculateStackResistance(std::vector<T> resistances) override
    {
        T invSum = 0.0;
        std::for_each(resistances.begin(), resistances.end(), [&](T res)
                      { invSum += 1.0 / res; });
        return 1 / invSum;
    }

    T computeCouplingCurrentDensity(T currentDensity, CVector<T> m1, CVector<T> m2, CVector<T> p) override
    {
        const T m1Comp = c_dot(m1, p);
        const T m2Comp = c_dot(m2, p);
        const T coupledI = currentDensity * this->couplingStrength * (m1Comp - m2Comp);
        return coupledI;
    }

public:
    ParallelStack(std::vector<Junction<T>> jL) : Stack<T>(jL) {}
};
#endif