#ifndef JUNCTION_H
#define JUNCTION_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdio.h>
#include <string>
#include <tuple>
#include <vector>

#include <fftw3.h>

#include "cvector.hpp"
#include "drivers.hpp"

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 221000.0
#define TtoAm 795774.715459
#define HBAR 6.62607015e-34 / (2. * M_PI)
#define ELECTRON_CHARGE 1.60217662e-19
#define BOLTZMANN_CONST 1.380649e-23
#define PERGYR MAGNETIC_PERMEABILITY *GYRO

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

CVector calculate_tensor_interaction(CVector m,
                                     std::vector<CVector> tensor,
                                     double Ms)
{
    CVector res(
        -Ms * tensor[0][0] * m[0] - Ms * tensor[0][1] * m[1] - Ms * tensor[0][2] * m[2],
        -Ms * tensor[1][0] * m[0] - Ms * tensor[1][1] * m[1] - Ms * tensor[1][2] * m[2],
        -Ms * tensor[2][0] * m[0] - Ms * tensor[2][1] * m[1] - Ms * tensor[2][2] * m[2]);
    return res;
}

CVector c_cross(CVector a, CVector b)
{
    CVector res(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);

    return res;
}

double c_dot(CVector a, CVector b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class EnergyDriver
{
public:
    static double calculateZeemanEnergy(CVector mag, CVector Hext, double cellVolume, double Ms)
    {
        return -MAGNETIC_PERMEABILITY * Ms * c_dot(mag, Hext) * cellVolume;
    }

    static double calculateAnisotropyEnergy(CVector mag, CVector anis, double K, double cellVolume)
    {
        const double sinSq = 1 - pow(c_dot(mag, anis) / (anis.length() * mag.length()), 2);
        return K * sinSq * cellVolume;
    }

    static double calculateIECEnergy(CVector mag, CVector other, double J, double cellSurface)
    {
        return -c_dot(mag, other) * J * cellSurface;
    }

    static double calculateDemagEnergy(CVector mag, CVector Hdemag, double Ms, double cellVolume)
    {
        return -0.5 * MAGNETIC_PERMEABILITY * Ms * c_dot(mag, Hdemag) * cellVolume;
    }
};
static std::default_random_engine generator;
class Layer
{
private:
    ScalarDriver currentDriver, IECDriver, anisotropyDriver;
    AxialDriver externalFieldDriver;

public:
    double cellVolume, cellSurface = 0.0;

    std::string id;

    CVector H_log, Hconst, mag, anis;
    CVector Hext, Hdipole, Hdemag, HIEC, HAnis, Hthermal;

    CVector IFlow = {0., 0., 1.0};
    double cellRadius = 35e-9;

    double Ms;
    double J_log = 0.0, K_log = 0.0, I_log = 0.0;
    double thickness;

    std::vector<CVector>
        demagTensor,
        dipoleTensor;

    // resting temperature in Kelvin
    double temperature;

    double Hstart = 0.0, Hstop = 0.0, Hstep = 0.0;

    bool includeSTT = false;

    // LLG params
    double damping;
    double SlonczewskiSpacerLayerParameter;
    double beta; // usually either set to 0 or to damping
    double spinPolarisation;

    std::normal_distribution<double> distribution;
    const double stochasticTorqueMean = 0.0;
    Layer() {}
    Layer(std::string id,
          CVector mag,
          CVector anis,
          double Ms,
          double thickness,
          double cellSurface,
          std::vector<CVector> demagTensor,
          std::vector<CVector> dipoleTensor,
          double temperature = 0.0,
          bool includeSTT = false,
          double damping = 0.011,
          double currentDensity = 1.,
          double SlonczewskiSpacerLayerParameter = 1.0,
          double beta = 0.0,
          double spinPolarisation = 0.8) : id(id),
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

        const double torqueStd = getStochasticStandardDeviation();
        std::cout << "Torque Std: " << torqueStd << std::endl;

        this->distribution = std::normal_distribution<double>(stochasticTorqueMean, torqueStd);
    }

    double getStochasticStandardDeviation()
    {
        const double scaledDamping = 2 * (this->damping) / (1 + pow(this->damping, 2));
        const double mainFactor = BOLTZMANN_CONST * this->temperature / (GYRO * this->Ms * this->cellVolume);
        return sqrt(scaledDamping * mainFactor);
    }

    void setCurrentDriver(ScalarDriver &driver)
    {
        this->currentDriver = driver;
    }

    void setExternalFieldDriver(AxialDriver &driver)
    {
        this->externalFieldDriver = driver;
    }
    void setAnisotropyDriver(ScalarDriver &driver)
    {
        this->anisotropyDriver = driver;
    }

    void setIECDriver(ScalarDriver &driver)
    {
        this->IECDriver = driver;
    }

    CVector calculateHeff(double time, double timeStep, CVector otherMag)
    {
        CVector Heff = {0., 0., 0.};

        this->Hext = calculateExternalField(time);
        this->Hdipole = calculate_tensor_interaction(otherMag, this->dipoleTensor, this->Ms);
        this->Hdemag = calculate_tensor_interaction(this->mag, this->demagTensor, this->Ms);
        this->HIEC = calculateIEC(time, otherMag);
        this->HAnis = calculateAnisotropy(time);
        this->Hthermal = calculateStochasticThermalField(timeStep);
        Heff = this->Hext +  // external
               this->HAnis + // anistotropy
               this->HIEC +  // IEC
               // demag
               // check the interaction here to be sure
               this->Hdemag +
               // dipole
               this->Hdipole +
               // stochastic field dependent on the temperature
               this->Hthermal;
        //    calculateStochasticThermalField(timeStep);

        return Heff;
    }

    CVector calculateStochasticThermalField(double timeStep)
    {
        // becomes zero if the temperature is 0
        CVector res(this->distribution, generator);
        res.normalize();
        // const double nom = sqrt((2 * this->damping * BOLTZMANN_CONST * this->temperature) /
        // (MAGNETIC_PERMEABILITY * GYRO * this->cellVolume * this->Ms * timeStep));
        // return res * nom;
        return res;
        // return CVector(0., 0., 0.);
    }

    CVector calculateExternalField(double time)
    {
        this->H_log =
            this->externalFieldDriver.getCurrentAxialDrivers(time);
        return this->H_log;
    }

    CVector calculateAnisotropy(double time)
    {
        this->K_log = this->anisotropyDriver.getCurrentScalarValue(time);
        const double nom = (2 * this->K_log) * c_dot(this->anis, this->mag) / (MAGNETIC_PERMEABILITY * this->Ms);
        return this->anis * nom;
    }

    CVector calculateIEC(double time, CVector coupledMag)
    {
        this->J_log = this->IECDriver.getCurrentScalarValue(time);
        const double nom = this->J_log / (MAGNETIC_PERMEABILITY * this->Ms * this->thickness);
        return (coupledMag - this->mag) * nom;
    }

    CVector llg(double time, CVector m, CVector coupledMag, CVector heff, double timeStep)
    {
        CVector prod, prod2, dmdt;
        // heff = calculateHeff(time, timeStep, coupledMag);
        prod = c_cross(m, heff);
        prod2 = c_cross(m, prod);
        dmdt = prod * -GYRO - prod2 * GYRO * this->damping;
        if (this->includeSTT)
        {
            // we will use coupledMag as the reference layer
            CVector prod3;
            this->I_log = this->currentDriver.getCurrentScalarValue(time);
            // damping-like torque factor
            const double aJ = HBAR * this->I_log /
                              (ELECTRON_CHARGE * MAGNETIC_PERMEABILITY * this->Ms * this->thickness);

            const double slonSq = pow(this->SlonczewskiSpacerLayerParameter, 2);
            // field like
            const double eta = (this->spinPolarisation * slonSq) / (slonSq + 1 + (slonSq - 1) * c_dot(m, coupledMag));
            const double sttTerm = GYRO * aJ * eta;

            prod3 = c_cross(m, coupledMag);
            // first term is "damping-like torque"
            // second term is "field-like torque"
            dmdt += c_cross(m, prod3) * -sttTerm + prod3 * sttTerm * this->beta;
        }
        return dmdt;
    }

    void rk4_step(double time, double timeStep, CVector coupledMag)
    {
        CVector k1, k2, k3, k4, m_t, heff;
        m_t = mag;
        heff = calculateHeff(time, timeStep, coupledMag);
        k1 = llg(time, m_t, coupledMag, heff, timeStep) * timeStep;
        k2 = llg(time + 0.5 * timeStep, m_t + k1 * 0.5, coupledMag, heff, timeStep) * timeStep;
        k3 = llg(time + 0.5 * timeStep, m_t + k2 * 0.5, coupledMag, heff, timeStep) * timeStep;
        k4 = llg(time + timeStep, m_t + k3, coupledMag, heff, timeStep) * timeStep;
        m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        mag = m_t;
    }

    double calculateLayerCriticalSwitchingCurrent(std::string plane)
    {
        const double fixedTerm = (2 * this->damping * ELECTRON_CHARGE * MAGNETIC_PERMEABILITY) / (HBAR * this->beta);
        const double layerParamTerm = (this->Ms * this->thickness);
        // static Hk
        const double Hk = calculateAnisotropy(0).length();
        if (plane == "IP")
        {
            return fixedTerm * layerParamTerm * Hk;
        }
        else if (plane == "PP")
        {
            const double Hdemag = calculate_tensor_interaction(this->mag, this->dipoleTensor, this->Ms).length();
            return fixedTerm * layerParamTerm * (Hdemag / 2 + Hk);
        }
        std::cout << "Warning -- UNKNOWN ENUM: " << plane << std::endl;
        return 0.0;
    }

    void calculatePinningCurrents()
    {
        const double qh = ELECTRON_CHARGE / HBAR;
        const double PMA_Ip = 6 * qh * this->damping * BOLTZMANN_CONST * this->temperature;
        // const double locDemag =
        // const double IMA_Ip = 2 * qh * sqrt(2 / M_PI) * sqrt(this->Hdemag * this->Ms * this->cellVolume * BOLTZMANN_CONST * this->temperature);

        std::cout << "PMA pinning current: " << PMA_Ip << std::endl;
        // std::cout << "IMA pinning current: " << IMA_Ip << std::endl;
    }
};

class Junction
{
    friend class Layer;
    std::vector<std::string> vectorNames = {"x", "y", "z"};

public:
    std::vector<Layer> layers;
    double Rp, Rap = 0.0;
    std::map<std::string, std::vector<double>> log;
    std::string fileSave;
    unsigned int logLength = 0;

    Junction(std::vector<Layer> layersToSet, std::string filename, double Rp = 100, double Rap = 200)
    {
        this->layers = std::move(layersToSet);
        this->Rp = Rp;
        this->Rap = Rap;
        this->fileSave = filename;
    }

    void clearLog()
    {
        this->log.clear();
        this->logLength = 0;
    }

    std::map<std::string, std::vector<double>> getLog()
    {
        return this->log;
    }

    Layer &findLayerByID(std::string lID)
    {
        for (Layer &l : layers)
        {
            if (l.id == lID)
            {
                return l;
            }
        }
        throw std::runtime_error("Failed to find the specified layer!");
    }

    double calculateMagnetoresistance(double cosTheta)
    {
        return this->Rp + (((this->Rap - this->Rp) / 2.0) * (1.0 - cosTheta));
    }

    typedef void (Layer::*scalarDriverSetter)(ScalarDriver &driver);
    typedef void (Layer::*axialDriverSetter)(AxialDriver &driver);
    void scalarlayerSetter(std::string layerID, scalarDriverSetter functor, ScalarDriver driver)
    {
        bool found = false;
        for (Layer &l : this->layers)
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
    void axiallayerSetter(std::string layerID, axialDriverSetter functor, AxialDriver driver)
    {
        bool found = false;
        for (Layer &l : this->layers)
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

    void setLayerExternalFieldDriver(std::string layerID, AxialDriver driver)
    {
        axiallayerSetter(layerID, &Layer::setExternalFieldDriver, driver);
    }
    void setLayerCurrentDriver(std::string layerID, ScalarDriver driver)
    {
        scalarlayerSetter(layerID, &Layer::setCurrentDriver, driver);
    }
    void setLayerAnisotropyDriver(std::string layerID, ScalarDriver driver)
    {
        scalarlayerSetter(layerID, &Layer::setAnisotropyDriver, driver);
    }
    void setLayerIECDriver(std::string layerID, ScalarDriver driver)
    {
        scalarlayerSetter(layerID, &Layer::setIECDriver, driver);
    }

    void logLayerParams(double t, double magnetoresistance, bool calculateEnergies = false)
    {
        for (Layer &layer : this->layers)
        {
            for (int i = 0; i < 3; i++)
            {
                this->log[layer.id + "_m" + vectorNames[i]].push_back(layer.mag[i]);
                this->log[layer.id + "_Hext" + vectorNames[i]].push_back(layer.H_log[i]);
            }
            this->log[layer.id + "_K"].push_back(layer.K_log);
            if (layer.includeSTT)
                this->log[layer.id + "_I"].push_back(layer.I_log);

            if (calculateEnergies)
            {
                this->log[layer.id + "_EZeeman"].push_back(EnergyDriver::calculateZeemanEnergy(layer.mag,
                                                                                               layer.Hext,
                                                                                               layer.cellVolume,
                                                                                               layer.Ms));
                this->log[layer.id + "_EAnis"].push_back(EnergyDriver::calculateAnisotropyEnergy(layer.mag,
                                                                                                 layer.anis,
                                                                                                 layer.K_log,
                                                                                                 layer.cellVolume));
                // this->log[layer.id + "_EIEC"] = EnergyDriver::calculateDemagEnergy(layer.mag,
                //                                                                    layer.other,
                //                                                                    layer.J_log,
                //                                                                    layer.cellSurface);
                this->log[layer.id + "_EDemag"].push_back(EnergyDriver::calculateDemagEnergy(layer.mag,
                                                                                             layer.Hdemag,
                                                                                             layer.Ms,
                                                                                             layer.cellVolume));
                this->log[layer.id + "_EDipole"].push_back(EnergyDriver::calculateDemagEnergy(layer.mag,
                                                                                              layer.Hdipole,
                                                                                              layer.Ms,
                                                                                              layer.cellVolume));
            }
        }

        this->log["R_free_bottom"].push_back(magnetoresistance);
        this->log["time"].push_back(t);

        this->logLength++;
    }

    void
    saveLogs()
    {
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

    std::map<std::string, double> calculateVoltageSpinDiode(double frequency, double power = 10e-6, const double minTime = 10e-9)
    {
        if (this->log.empty())
        {
            throw std::runtime_error("Empty log! Cannot proceed without running a simulation!");
        }
        const double omega = 2 * M_PI * frequency;
        const std::string res = "R_free_bottom";
        std::vector<double> &resistance = this->log[res];
        auto it = std::find_if(this->log["time"].begin(), this->log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });
        // turn into index
        const int thresIdx = (int)(this->log["time"].end() - it);
        const int cutSize = this->log["time"].size() - thresIdx;
        // Rpp
        const double RppMax = *std::max_element(resistance.begin() + thresIdx, resistance.end());
        const double RppMin = *std::min_element(resistance.begin() + thresIdx, resistance.end());
        const double avgR = std::accumulate(resistance.begin() + thresIdx, resistance.end(), 0.0) / cutSize;
        const double Iampl = sqrt(power / avgR);
        std::vector<double> voltage, current;
        std::transform(
            this->log["time"].begin() + thresIdx, this->log["time"].end(),
            std::back_inserter(current),
            [&Iampl, &omega](const double &time) { return Iampl * sin(omega * time); });

        for (unsigned int i = 0; i < cutSize; i++)
        {
            voltage.push_back(resistance[thresIdx + i] * current[i]);
        }
        const double Vmix = std::accumulate(voltage.begin(), voltage.end(), 0.0) / voltage.size();
        std::map<std::string, double> mRes = {{"Vmix", Vmix}, {"RMax", RppMax}, {"RMin", RppMin}, {"Rpp", (RppMax - RppMin)}};
        return mRes;
    }

    std::map<std::string, double>
    calculateFFT(double minTime = 10.0e-9, double timeStep = 1e-11)
    {

        if (this->log.empty())
        {
            throw std::runtime_error("Empty log! Cannot proceed without running a simulation!");
        }
        // std::cout << "FFT calculation" << std::endl;
        auto it = std::find_if(this->log["time"].begin(), this->log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });
        // turn into index
        const int thresIdx = (int)(this->log["time"].end() - it);
        const int cutSize = this->log["time"].size() - thresIdx;

        // std::cout << "Initiating FFT plan execution " << std::endl;

        // plan creation is not thread safe
        const double normalizer = timeStep * cutSize;
        const int maxIt = (cutSize % 2) ? cutSize / 2 : (cutSize - 1) / 2;
        std::vector<double> frequencySteps(maxIt);
        frequencySteps[0] = 0;
        for (int i = 1; i <= maxIt; i++)
        {
            frequencySteps[i - 1] = (i - 1) / normalizer;
        }

        std::map<std::string, std::vector<double>> result;
        std::map<std::string, double> maxAmpls;
        for (const auto &magTag : this->vectorNames)
        {
            // std::cout << "Doing FFT for: " << magTag << std::endl;
            std::vector<double> cutMag(this->log["free_m" + magTag].begin() + thresIdx, this->log["free_m" + magTag].end());
            // std::cout << "Sub size: " << cutMag.size() << std::endl;
            fftw_complex out[cutMag.size()];
            // define FFT plan
            fftw_plan plan = fftw_plan_dft_r2c_1d(cutMag.size(),
                                                  cutMag.data(),
                                                  out,
                                                  FFTW_ESTIMATE); // here it's weird, FFT_FORWARD produces an empty plan
            if (plan == NULL)
            {
                throw std::runtime_error("Plan creation for fftw failed, cannot proceed");
            }
            fftw_execute(plan);
            std::vector<double> amplitudes, phases;

            const int outBins = (cutMag.size() + 1) / 2;
            amplitudes.reserve(outBins);
            phases.reserve(outBins);

            double maxAmpl = 0.0;
            double maxFreq = 0.0;
            double maxPhase = 0.0;
            amplitudes[0] = out[0][0];
            phases[0] = 0;
            for (int i = 1; i < outBins; i++)
            {
                const auto tandem = out[i];
                const double real = tandem[0];
                const double img = tandem[1];
                const double ampl = sqrt(pow(real, 2) + pow(img, 2));
                phases[i] = tan(img / real);
                amplitudes[i] = ampl;
                if (ampl >= maxAmpl)
                {
                    maxAmpl = ampl;
                    maxFreq = frequencySteps[i];
                    maxPhase = phases[i];
                }
            }
            // std::cout << "Max frequency for the system is " << maxFreq << " with " << maxAmpl << std::endl;
            fftw_destroy_plan(plan);
            result[magTag + "_amplitude"] = amplitudes;
            result[magTag + "_phase"] = phases;
            maxAmpls[magTag + "_resonant"] = maxFreq;
            maxAmpls[magTag + "_amplitude"] = maxAmpl;
            maxAmpls[magTag + "_phase"] = maxPhase;
        }
        return maxAmpls;
    }

    void runSingleRK4Iteration(double t, double timeStep)
    {
        CVector l1mag = this->layers[0].mag;
        CVector l2mag = this->layers[1].mag;
        this->layers[0].rk4_step(
            t, timeStep, l2mag);
        this->layers[1].rk4_step(
            t, timeStep, l1mag);
    }

    double getMagnetoresistance()
    {
        return calculateMagnetoresistance(c_dot(layers[0].mag, layers[1].mag));
    }

    void runSimulation(double totalTime, double timeStep = 1e-13, double writeFrequency = 1e-11,
                       bool persist = true, bool log = false, bool calculateEnergies = false)
    {

        const unsigned int totalIterations = (int)(totalTime / timeStep);
        double t;
        const unsigned int writeEvery = (int)(writeFrequency / timeStep) - 1;

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (unsigned int i = 0; i < totalIterations; i++)
        {
            t = i * timeStep;

            runSingleRK4Iteration(t, timeStep);

            if (!(i % writeEvery))
            {
                const double magRes = getMagnetoresistance();
                logLayerParams(t, magRes, calculateEnergies);
            }
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if (persist)
            saveLogs();
        if (log)
            std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    }

    static std::vector<CVector> RK45(CVector mag, CVector mTop, CVector mBottom, CVector Hext, int layer, double dt, CVector HOe,
                                     std::vector<double> Ms,
                                     std::vector<double> Ku,
                                     std::vector<double> Ju,
                                     std::vector<CVector> Kdir,
                                     std::vector<double> th,
                                     std::vector<double> alpha,
                                     std::vector<CVector> demag);

    static CVector LLG(CVector mag, CVector mTop, CVector mBottom, CVector Hext, CVector HOe, int layer,
                       std::vector<double> Ms,
                       std::vector<double> Ku,
                       std::vector<double> Ju,
                       std::vector<CVector> Kdir,
                       std::vector<double> th,
                       std::vector<double> alpha,
                       std::vector<CVector> demag);
};

#endif