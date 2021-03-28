#ifndef JUNCTION_H
#define JUNCTION_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdio.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "cvector.hpp"
#include "drivers.hpp"

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 220880.0
#define TtoAm 795774.715459
#define HBAR 6.62607015e-34 / (2. * M_PI)
#define ELECTRON_CHARGE 1.60217662e-19
#define BOLTZMANN_CONST 1.380649e-23

typedef CVector<double> DVector;
typedef CVector<float> FVector;

double operator"" _ns(unsigned long long timeUnit)
{
    return ((double)timeUnit) / 1e9;
}
double operator"" _ns(long double timeUnit)
{
    return ((double)timeUnit) / 1e9;
}

double operator"" _mT(unsigned long long tesla)
{
    return ((double)tesla) / 1000.0;
}

double operator"" _mT(long double tesla)
{
    return ((double)tesla) / 1000.0;
}

template <typename T>
inline CVector<T> calculate_tensor_interaction(CVector<T> &m,
                                               std::vector<CVector<T>> &tensor,
                                               T &Ms)
{
    CVector<T> res(
        tensor[0][0] * m[0] + tensor[0][1] * m[1] + tensor[0][2] * m[2],
        tensor[1][0] * m[0] + tensor[1][1] * m[1] + tensor[1][2] * m[2],
        tensor[2][0] * m[0] + tensor[2][1] * m[1] + tensor[2][2] * m[2]);
    return res * (Ms / MAGNETIC_PERMEABILITY);
}

template <typename T>
inline CVector<T> c_cross(const CVector<T> &a, const CVector<T> &b)
{
    CVector<T> res(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);

    return res;
}

template <typename T>
inline T c_dot(CVector<T> &a, CVector<T> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
class EnergyDriver
{
public:
    static T calculateZeemanEnergy(CVector<T> mag, CVector<T> Hext, T cellVolume, T Ms)
    {
        return -MAGNETIC_PERMEABILITY * Ms * c_dot<T>(mag, Hext) * cellVolume;
    }

    static T calculateAnisotropyEnergy(CVector<T> mag, CVector<T> anis, CVector<T> K, T cellVolume)
    {
        const T sinSq = 1 - pow(c_dot<T>(mag, anis) / (anis.length() * mag.length()), 2);
        return K * sinSq * cellVolume;
    }

    static T calculateIECEnergy(CVector<T> mag, CVector<T> other, T J, T cellSurface)
    {
        return -c_dot<T>(mag, other) * J * cellSurface;
    }

    static T calculateDemagEnergy(CVector<T> mag, CVector<T> Hdemag, T Ms, T cellVolume)
    {
        return -0.5 * MAGNETIC_PERMEABILITY * Ms * c_dot<T>(mag, Hdemag) * cellVolume;
    }
};
static std::default_random_engine generator;
template <typename T>
class Layer
{
private:
    ScalarDriver<T> currentDriver = NullDriver<T>();
    ScalarDriver<T> IECDriver = NullDriver<T>();
    ScalarDriver<T> anisotropyDriver = NullDriver<T>();
    AxialDriver<T> externalFieldDriver = NullAxialDriver<T>();
    AxialDriver<T> HoeDriver = NullAxialDriver<T>();

public:
    T cellVolume = 0.0, cellSurface = 0.0;

    std::string id;

    CVector<T> H_log, Hoe_log, Hconst, mag, anisAxis, anis;
    CVector<T> Hext, Hdipole, Hdemag, HIEC, Hoe, HAnis, Hthermal, Hfl;
    T K_log = 0.0;
    T Ms = 0.0;
    T J_log = 0.0, I_log = 0.0;
    T thickness = 0.0;

    std::vector<CVector<T>>
        demagTensor,
        dipoleTensor;

    // resting temperature in Kelvin
    T temperature;

    T Hstart = 0.0, Hstop = 0.0, Hstep = 0.0;

    bool includeSTT = false;

    // LLG params
    T damping;
    T SlonczewskiSpacerLayerParameter;
    T beta; // usually either set to 0 or to damping
    T spinPolarisation;

    std::normal_distribution<T> distribution;
    const T stochasticTorqueMean = 0.0;
    Layer() {}
    Layer(std::string id,
          CVector<T> mag,
          CVector<T> anis,
          T Ms,
          T thickness,
          T cellSurface,
          std::vector<CVector<T>> demagTensor,
          std::vector<CVector<T>> dipoleTensor,
          T temperature = 0.0,
          bool includeSTT = false,
          T damping = 0.011,
          T SlonczewskiSpacerLayerParameter = 1.0,
          T beta = 0.0,
          T spinPolarisation = 0.8,
          bool silent = true) : id(id),
                                mag(mag),
                                anis(anis),
                                Ms(Ms),
                                thickness(thickness),
                                cellSurface(cellSurface),
                                demagTensor(demagTensor),
                                dipoleTensor(dipoleTensor),
                                temperature(temperature),
                                includeSTT(includeSTT),
                                damping(damping),
                                SlonczewskiSpacerLayerParameter(SlonczewskiSpacerLayerParameter),
                                beta(beta),
                                spinPolarisation(spinPolarisation)
    {
        this->cellVolume = this->cellSurface * this->thickness;
        // this is Langevin fluctuation field from the Kaiser paper
        const T torqueStd = this->getLangevinStochasticStandardDeviation();
        if (!silent)
        {
            std::cout << "Langevin torque std: " << torqueStd << std::endl;
            std::cout << "Cell Volume: " << cellVolume << std::endl;
        }
        this->distribution = std::normal_distribution<T>(stochasticTorqueMean, torqueStd);
    }

    T getLangevinStochasticStandardDeviation()
    {
        if (this->cellVolume == 0.0)
            return 0.0;
        const T scaledDamping = 2 * (this->damping) / (1 + pow(this->damping, 2));
        const T mainFactor = BOLTZMANN_CONST * this->temperature / (GYRO * this->Ms * this->cellVolume);
        return sqrt(scaledDamping * mainFactor);
    }

    void setCurrentDriver(ScalarDriver<T> &driver)
    {
        this->currentDriver = driver;
    }

    void setAnisotropyDriver(ScalarDriver<T> &driver)
    {
        this->anisotropyDriver = driver;
    }

    void setIECDriver(ScalarDriver<T> &driver)
    {
        this->IECDriver = driver;
    }

    void setExternalFieldDriver(AxialDriver<T> &driver)
    {
        this->externalFieldDriver = driver;
    }
    void setLayerOerstedFieldDriver(AxialDriver<T> &driver)
    {
        this->HoeDriver = driver;
    }

    void setMagnetisation(CVector<T> &mag)
    {
        this->mag = mag;
    }

    const CVector<T> calculateHeff(T time, T timeStep, CVector<T> &stepMag, CVector<T> &otherMag)
    {
        this->Hext = calculateExternalField(time);
        this->Hoe = calculateHOeField(time);
        this->Hdipole = calculate_tensor_interaction(otherMag, this->dipoleTensor, this->Ms);
        this->Hdemag = calculate_tensor_interaction(stepMag, this->demagTensor, this->Ms);
        this->HIEC = calculateIEC(time, stepMag, otherMag);
        this->HAnis = calculateAnisotropy(stepMag, time);
        this->Hfl = calculateLangevinStochasticField(timeStep);
        const CVector<T> Heff = this->Hext    // external
                                + this->HAnis // anistotropy
                                + this->HIEC  // IEC
                                + this->Hoe   // Oersted field
                                // demag -- negative contribution
                                - this->Hdemag
                                // dipole -- negative contribution
                                - this->Hdipole
                                // stochastic field dependent on the temperature
                                + this->Hfl;

        return Heff;
    }

    CVector<T> calculateHOeField(T &time)
    {
        this->Hoe_log = this->HoeDriver.getCurrentAxialDrivers(time);
        return this->Hoe_log;
    }

    CVector<T> calculateLangevinStochasticField(T &timeStep)
    {
        if (this->cellVolume == 0.0)
            return CVector<T>();
        // becomes zero if the temperature is 0
        CVector<T> res(this->distribution, generator);
        // either of those expressions may be correct -- undecided for now
        const T nom = sqrt((2 * this->damping * BOLTZMANN_CONST * this->temperature) /
                           (MAGNETIC_PERMEABILITY * GYRO * this->cellVolume * this->Ms * timeStep));
        return res * nom;
        // return res;
    }

    CVector<T> calculateExternalField(T &time)
    {
        this->H_log =
            this->externalFieldDriver.getCurrentAxialDrivers(time);
        return this->H_log;
    }

    CVector<T> calculateAnisotropy(CVector<T> &stepMag, T &time)
    {
        // this->K_log = this->anisotropyDriver.getCurrentAxialDrivers(time);
        // return this->K_log * c_dot<T>(this->anisAxis, this->mag) * 2 / this->Ms;
        this->K_log = this->anisotropyDriver.getCurrentScalarValue(time);
        const T nom = (2 * this->K_log) * c_dot<T>(this->anis, stepMag) / (this->Ms);
        return this->anis * nom;
    }

    CVector<T> calculateIEC(T &time, CVector<T> &stepMag, CVector<T> &coupledMag)
    {
        this->J_log = this->IECDriver.getCurrentScalarValue(time);
        const T nom = this->J_log / (this->Ms * this->thickness);
        return (coupledMag - stepMag) * nom;
        // either of those may be correct -- undecided for now
        // return coupledMag * nom;
    }

    const CVector<T> llg(T time, CVector<T> m, CVector<T> &coupledMag, T &timeStep)
    {
        const CVector<T> heff = calculateHeff(time, timeStep, m, coupledMag);
        const CVector<T> prod = c_cross<T>(m, heff);
        const CVector<T> prod2 = c_cross<T>(m, prod);
        if (this->includeSTT)
        {
            // we will use coupledMag as the reference layer
            this->I_log = this->currentDriver.getCurrentScalarValue(time);
            // damping-like torque factor
            const T aJ = HBAR * this->I_log /
                         (ELECTRON_CHARGE * this->Ms * this->thickness);

            const T slonSq = pow(this->SlonczewskiSpacerLayerParameter, 2);
            // field like
            const T eta = (this->spinPolarisation * slonSq) / (slonSq + 1 + (slonSq - 1) * c_dot<T>(m, coupledMag));
            const T sttTerm = GYRO * aJ * eta;

            const CVector<T> prod3 = c_cross<T>(m, coupledMag);
            // first term is "damping-like torque"
            // second term is "field-like torque"
            CVector<T> dmdt = (prod * -GYRO) - (prod2 * GYRO * this->damping) + c_cross<T>(m, prod3) * -sttTerm + prod3 * sttTerm * this->beta;
            return dmdt;
        }
        CVector<T> dmdt = (prod * -GYRO) - (prod2 * GYRO * this->damping);
        return dmdt;
    }

    void rk4_step(T &time, T &timeStep, CVector<T> &coupledMag)
    {
        CVector<T> m_t = this->mag;
        const CVector<T> k1 = llg(time, m_t, coupledMag, timeStep) * timeStep;
        const CVector<T> k2 = llg(time + 0.5 * timeStep, m_t + k1 * 0.5, coupledMag, timeStep) * timeStep;
        const CVector<T> k3 = llg(time + 0.5 * timeStep, m_t + k2 * 0.5, coupledMag, timeStep) * timeStep;
        const CVector<T> k4 = llg(time + timeStep, m_t + k3, coupledMag, timeStep) * timeStep;
        m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        this->mag = m_t;
    }
};

template <typename T>
class Junction
{
    friend class Layer<T>;
    const std::vector<const std::string> vectorNames = {"x", "y", "z"};

public:
    enum MRmode
    {
        NONE = 0,
        CLASSIC = 1,
        STRIP = 2
    };
    MRmode MR_mode;
    std::vector<Layer<T>> layers;
    T Rp, Rap = 0.0;

    std::vector<T> Rx0, Ry0, AMR_X, AMR_Y, SMR_X, SMR_Y, AHE;
    std::unordered_map<std::string, std::vector<T>> log;
    std::string fileSave;
    unsigned int logLength = 0;
    int layerNo;
    Junction() {}
    Junction(std::vector<Layer<T>> layersToSet, std::string filename = "")
    {
        this->MR_mode = NONE;
        this->layers = std::move(layersToSet);
        this->layerNo = this->layers.size();
        if (this->layerNo == 0)
        {
            throw std::invalid_argument("Passed a zero length Layer vector!");
        }
        this->fileSave = std::move(filename);
    }
    explicit Junction(std::vector<Layer<T>> layersToSet, std::string filename, T Rp, T Rap) : Junction(
                                                                                                  layersToSet, filename)
    {
        if (this->layerNo != 2)
        {
            throw std::invalid_argument("This constructor supports only bilayers! Choose the other one with the strip resistance!");
        }
        this->Rp = Rp;
        this->Rap = Rap;
        this->MR_mode = CLASSIC;
    }

    /**
     * Creates a junction with a STRIP magnetoresistance.  
     * Each of the Rx0, Ry, AMR, AMR and SMR is list matching the 
     * length of the layers passed (they directly correspond to each layer).
     * @param Rx0
     * @param Ry0
     * @param AMR_X
     * @param AMR_Y
     * @param SMR_X
     * @param SMR_Y
     * @param AHE
     */
    explicit Junction(std::vector<Layer<T>> layersToSet,
                      std::string filename,
                      std::vector<T> Rx0,
                      std::vector<T> Ry0,
                      std::vector<T> AMR_X,
                      std::vector<T> AMR_Y,
                      std::vector<T> SMR_X,
                      std::vector<T> SMR_Y,
                      std::vector<T> AHE) : Rx0(std::move(Rx0)),
                                            Ry0(std::move(Ry0)),
                                            AMR_X(std::move(AMR_X)),
                                            AMR_Y(std::move(AMR_Y)),
                                            SMR_X(std::move(SMR_X)),
                                            SMR_Y(std::move(SMR_Y)),
                                            AHE(std::move(AHE))

    {
        this->layers = std::move(layersToSet);
        this->layerNo = this->layers.size();
        if (this->layerNo == 0)
        {
            throw std::invalid_argument("Passed a zero length Layer vector!");
        }
        if ((this->layerNo != this->Rx0.size()) ||
            (this->layerNo != this->Ry0.size()) ||
            (this->layerNo != this->AMR_X.size()) ||
            (this->layerNo != this->AMR_Y.size()) ||
            (this->layerNo != this->AHE.size()) ||
            (this->layerNo != this->SMR_X.size()) ||
            (this->layerNo != this->SMR_Y.size()))
        {
            throw std::invalid_argument("Layers and Rx0, Ry, AMR, AMR and SMR must be of the same size!");
        }
        this->fileSave = std::move(filename);
        this->MR_mode = STRIP;
    }

    void clearLog()
    {
        this->log.clear();
        this->logLength = 0;
    }

    std::unordered_map<std::string, std::vector<T>> &getLog()
    {
        return this->log;
    }

    Layer<T> &findLayerByID(std::string lID)
    {
        for (auto &l : layers)
        {
            if (l.id == lID)
            {
                return l;
            }
        }
        throw std::runtime_error("Failed to find the specified layer!");
    }

    typedef void (Layer<T>::*scalarDriverSetter)(ScalarDriver<T> &driver);
    typedef void (Layer<T>::*axialDriverSetter)(AxialDriver<T> &driver);
    void scalarlayerSetter(std::string &layerID, scalarDriverSetter functor, ScalarDriver<T> driver)
    {
        bool found = false;
        for (auto &l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                (l.*functor)(driver);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }
    void axiallayerSetter(std::string &layerID, axialDriverSetter functor, AxialDriver<T> driver)
    {
        bool found = false;
        for (auto &l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                (l.*functor)(driver);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }

    void setLayerAnisotropyDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setAnisotropyDriver, driver);
    }
    void setLayerExternalFieldDriver(std::string layerID, AxialDriver<T> driver)
    {
        axiallayerSetter(layerID, &Layer<T>::setExternalFieldDriver, driver);
    }
    void setLayerOerstedFieldDriver(std::string layerID, AxialDriver<T> driver)
    {
        axiallayerSetter(layerID, &Layer<T>::setLayerOerstedFieldDriver, driver);
    }
    void setLayerCurrentDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setCurrentDriver, driver);
    }
    void setLayerIECDriver(std::string layerID, ScalarDriver<T> driver)
    {
        scalarlayerSetter(layerID, &Layer<T>::setIECDriver, driver);
    }

    void setLayerMagnetisation(std::string layerID, CVector<T> &mag)
    {
        bool found = false;
        for (auto &l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                l.setMagnetisation(mag);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }

    void logLayerParams(T &t, bool calculateEnergies = false)
    {
        for (const auto &layer : this->layers)
        {
            const std::string lId = layer.id;
            this->log[lId + "_K"].emplace_back(layer.K_log);
            for (int i = 0; i < 3; i++)
            {
                this->log[lId + "_m" + vectorNames[i]].emplace_back(layer.mag[i]);
                this->log[lId + "_Hext" + vectorNames[i]].emplace_back(layer.H_log[i]);
                // this->log[lId + "_K" + vectorNames[i]].emplace_back(layer.K_log[i]);
            }

            if (layer.includeSTT)
                this->log[lId + "_I"].emplace_back(layer.I_log);

            if (calculateEnergies)
            {
                this->log[lId + "_EZeeman"].emplace_back(EnergyDriver<T>::calculateZeemanEnergy(layer.mag,
                                                                                             layer.Hext,
                                                                                             layer.cellVolume,
                                                                                             layer.Ms));
                // this->log[lId + "_EAnis"].push_back(EnergyDriver::calculateAnisotropyEnergy(layer.mag,
                //                                                                                  layer.anisAxis,
                //                                                                                  layer.K_log,
                //                                                                                  layer.cellVolume));
                // this->log[lId + "_EIEC"] = EnergyDriver::calculateDemagEnergy(layer.mag,
                //                                                                    layer.other,
                //                                                                    layer.J_log,
                //                                                                    layer.cellSurface);
                this->log[lId + "_EDemag"].emplace_back(EnergyDriver<T>::calculateDemagEnergy(layer.mag,
                                                                                           layer.Hdemag,
                                                                                           layer.Ms,
                                                                                           layer.cellVolume));
                this->log[lId + "_EDipole"].emplace_back(EnergyDriver<T>::calculateDemagEnergy(layer.mag,
                                                                                            layer.Hdipole,
                                                                                            layer.Ms,
                                                                                            layer.cellVolume));
            }
        }
        if (MR_mode == CLASSIC)
        {
            const auto magnetoresistance = calculateMagnetoresistance(c_dot<T>(this->layers[0].mag, this->layers[1].mag));
            this->log["R_free_bottom"].emplace_back(magnetoresistance);
        }
        else if (MR_mode == STRIP)
        {
            const auto magnetoresistance = stripMagnetoResistance(this->Rx0, this->Ry0,
                                                                  this->AMR_X,
                                                                  this->SMR_X,
                                                                  this->AMR_Y,
                                                                  this->SMR_Y,
                                                                  this->AHE);
            this->log["Rx"].emplace_back(magnetoresistance[0]);
            this->log["Ry"].emplace_back(magnetoresistance[1]);
        }
        this->log["time"].emplace_back(t);
        this->logLength++;
    }

    void
    saveLogs()
    {
        if (this->fileSave == "")
        {
            // if there's an empty fn, don't save
            std::cout << "Ignoring file save to an empty filename" << std::endl;
            return;
        }
        std::ofstream logFile;
        logFile.open(this->fileSave);
        for (const auto &keyPair : this->log)
        {
            logFile << keyPair.first << ";";
        }
        logFile << "\n";
        for (unsigned int i = 0; i < logLength; i++)
        {
            for (const auto &keyPair : this->log)
            {
                logFile << keyPair.second[i] << ";";
            }
            logFile << "\n";
        }
        logFile.close();
    }

    void runSingleLayerRK4Iteration(T &t, T &timeStep)
    {
        /**
         * Single layer iteration. IEC interaction is turned off.
         * @param t: current time
         * @param timeStep: integration step
         * */
        CVector<T> null;
        this->layers[0].rk4_step(
            t, timeStep, null);
    }

    void runMultiLayerRK4Iteration(T &t, T &timeStep)
    {
        CVector<T> l1mag = this->layers[0].mag;
        CVector<T> l2mag = this->layers[1].mag;
        this->layers[0].rk4_step(
            t, timeStep, l2mag);
        this->layers[1].rk4_step(
            t, timeStep, l1mag);
    }

    std::vector<T> stripMagnetoResistance(std::vector<T> &Rx0,
                                          std::vector<T> &Ry0,
                                          std::vector<T> &AMR_X,
                                          std::vector<T> &SMR_X,
                                          std::vector<T> &AMR_Y,
                                          std::vector<T> &SMR_Y,
                                          std::vector<T> &AHE)
    {
        T Rx_acc = 0.0;
        T Ry_acc = 0.0;

        for (int i = 0; i < this->layers.size(); i++)
        {
            const T Rx = Rx0[i] + AMR_X[i] * pow(this->layers[i].mag.x, 2) + SMR_X[i] * pow(this->layers[i].mag.y, 2);
            const T Ry = Ry0[i] + 0.5 * AHE[i] * this->layers[i].mag.z +
                         (AMR_Y[i] - SMR_Y[i]) * this->layers[i].mag.x * this->layers[i].mag.y;
            Rx_acc += Rx;
            Ry_acc += Ry;
        }

        return {1 / Rx_acc, 1 / Ry_acc};
    }

    T calculateMagnetoresistance(T cosTheta)
    {
        return this->Rp + (((this->Rap - this->Rp) / 2.0) * (1.0 - cosTheta));
    }

    std::vector<T> getMagnetoresistance()
    {
        if (this->MR_mode == CLASSIC)
        {
            return {calculateMagnetoresistance(c_dot<T>(layers[0].mag, layers[1].mag))};
        }
        else
        {
            return stripMagnetoResistance(this->Rx0,
                                          this->Ry0,
                                          this->AMR_X,
                                          this->SMR_X,
                                          this->AMR_Y,
                                          this->SMR_Y,
                                          this->AHE);
        }
    }

    void runSimulation(T totalTime, T timeStep = 1e-13, T writeFrequency = 1e-11,
                       bool persist = true, bool log = false, bool calculateEnergies = false)
    {
        if (timeStep > writeFrequency)
        {
            std::runtime_error("The time step cannot be larger than write frequency!");
        }
        const unsigned int totalIterations = (int)(totalTime / timeStep);
        T t;
        const unsigned int writeEvery = (int)(writeFrequency / timeStep);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (unsigned int i = 0; i < totalIterations; i++)
        {
            t = i * timeStep;
            if (this->layerNo == 1)
            {
                runSingleLayerRK4Iteration(t, timeStep);
            }
            else
            {
                runMultiLayerRK4Iteration(t, timeStep);
            }

            if (!(i % writeEvery))
            {
                logLayerParams(t, calculateEnergies);
            }
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if (persist)
            saveLogs();
        if (log)
        {
            std::cout << "Steps in simulation: " << totalIterations << std::endl;
            std::cout << "Write every: " << writeEvery << std::endl;
            std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
        }
    }
};

#endif