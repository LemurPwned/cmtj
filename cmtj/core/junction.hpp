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
#include <vector>

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

    static double calculateAnisotropyEnergy(CVector mag, CVector anisAxis, CVector K, double cellVolume)
    {
        // const double sinSq = 1 - pow(c_dot(mag, anis) / (anis.length() * mag.length()), 2);
        // return K * sinSq * cellVolume;
        const double Kval = 2 * c_dot(mag, anisAxis) * c_dot(mag, K) * cellVolume;
        return Kval;
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
    ScalarDriver currentDriver = NullDriver();
    ScalarDriver IECDriver = NullDriver();
    AxialDriver anisotropyDriver = NullAxialDriver();
    AxialDriver externalFieldDriver = NullAxialDriver();
    AxialDriver HoeDriver = NullAxialDriver();

public:
    double cellVolume, cellSurface = 0.0;

    std::string id;

    CVector H_log, Hoe_log, Hconst, mag, anisAxis, K_log;
    CVector Hext, Hdipole, Hdemag, HIEC, Hoe, HAnis, Hthermal, Hfl;

    double Ms;
    double J_log = 0.0, I_log = 0.0;
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
          double Ms,
          double thickness,
          double cellSurface,
          std::vector<CVector> demagTensor,
          std::vector<CVector> dipoleTensor,
          double temperature = 0.0,
          bool includeSTT = false,
          double damping = 0.011,
          double SlonczewskiSpacerLayerParameter = 1.0,
          double beta = 0.0,
          double spinPolarisation = 0.8,
          bool silent = true) : id(id),
                                mag(mag),
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
        const double torqueStd = this->getLangevinStochasticStandardDeviation();
        if (!silent)
        {
            std::cout << "Langevin torque std: " << torqueStd << std::endl;
            std::cout << "Cell Volume: " << cellVolume << std::endl;
        }
        this->distribution = std::normal_distribution<double>(stochasticTorqueMean, torqueStd);
    }

    double getLangevinStochasticStandardDeviation()
    {
        const double scaledDamping = 2 * (this->damping) / (1 + pow(this->damping, 2));
        const double mainFactor = BOLTZMANN_CONST * this->temperature / (GYRO * this->Ms * this->cellVolume);
        return sqrt(scaledDamping * mainFactor);
    }

    void setCurrentDriver(ScalarDriver &driver)
    {
        this->currentDriver = driver;
    }

    void setAnisotropyDriver(AxialDriver &driver)
    {
        this->anisotropyDriver = driver;
        this->anisAxis = driver.getUnitAxis();
    }

    void setIECDriver(ScalarDriver &driver)
    {
        this->IECDriver = driver;
    }

    void setExternalFieldDriver(AxialDriver &driver)
    {
        this->externalFieldDriver = driver;
    }
    void setLayerOerstedFieldDriver(AxialDriver &driver)
    {
        this->HoeDriver = driver;
    }

    CVector calculateHeff(double time, double timeStep, CVector otherMag)
    {
        CVector Heff = {0., 0., 0.};

        this->Hext = calculateExternalField(time);
        this->Hoe = calculateHOeField(time);
        this->Hdipole = calculate_tensor_interaction(otherMag, this->dipoleTensor, this->Ms);
        this->Hdemag = calculate_tensor_interaction(this->mag, this->demagTensor, this->Ms);
        this->HIEC = calculateIEC(time, otherMag);
        this->HAnis = calculateAnisotropy(time);
        this->Hfl = calculateLangevinStochasticField(timeStep);
        Heff = this->Hext +  // external
               this->HAnis + // anistotropy
               this->HIEC +  // IEC
               this->Hoe +   // Oersted field
               // demag
               // check the interaction here to be sure
               this->Hdemag +
               // dipole
               this->Hdipole +
               // stochastic field dependent on the temperature
               this->Hfl;

        return Heff;
    }

    CVector calculateHOeField(double time)
    {
        this->Hoe_log = this->HoeDriver.getCurrentAxialDrivers(time);
        return this->Hoe_log;
    }

    CVector calculateLangevinStochasticField(double timeStep)
    {
        // becomes zero if the temperature is 0
        CVector res(this->distribution, generator);
        // either of those expressions may be correct -- undecided for now
        const double nom = sqrt((2 * this->damping * BOLTZMANN_CONST * this->temperature) /
                                (MAGNETIC_PERMEABILITY * GYRO * this->cellVolume * this->Ms * timeStep));
        return res * nom;
        // return res;
    }

    CVector calculateExternalField(double time)
    {
        this->H_log =
            this->externalFieldDriver.getCurrentAxialDrivers(time);
        return this->H_log;
    }

    CVector calculateAnisotropy(double time)
    {
        this->K_log = this->anisotropyDriver.getCurrentAxialDrivers(time);
        return this->K_log * c_dot(this->anisAxis, this->mag) * 2 / (MAGNETIC_PERMEABILITY * this->Ms);
        // this->K_log = this->anisotropyDriver.getCurrentScalarValue(time);
        // const double nom = (2 * this->K_log) * c_dot(this->anis, this->mag) / (MAGNETIC_PERMEABILITY * this->Ms);
        // return this->anis * nom;
    }

    CVector calculateIEC(double time, CVector coupledMag)
    {
        this->J_log = this->IECDriver.getCurrentScalarValue(time);
        const double nom = this->J_log / (MAGNETIC_PERMEABILITY * this->Ms * this->thickness);
        // return (coupledMag - this->mag) * nom;
        // either of those may be correct -- undecided for now
        return coupledMag * nom;
    }

    CVector llg(double time, CVector m, CVector coupledMag, CVector heffEx, double timeStep)
    {
        CVector prod, prod2, dmdt, heff;
        heff = calculateHeff(time, timeStep, coupledMag);
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
        // heff = calculateHeff(time, timeStep, coupledMag);
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
    enum MRmode
    {
        NONE = 0,
        CLASSIC = 1,
        STRIP = 2
    };
    MRmode MR_mode;
    std::vector<Layer> layers;
    double Rp, Rap = 0.0;

    std::vector<double> Rx0, Ry0, AMR_X, AMR_Y, SMR_X, SMR_Y, AHE;
    std::map<std::string, std::vector<double>> log;
    std::string fileSave;
    unsigned int logLength = 0;
    int layerNo;
    Junction() {}
    Junction(std::vector<Layer> layersToSet, std::string filename = "")
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
    Junction(std::vector<Layer> layersToSet, std::string filename, double Rp, double Rap) : Junction(
                                                                                                layersToSet, filename)
    {
        if (this->layerNo != 2)
        {
            throw std::invalid_argument("This constructor supports only bilayers! Choose the other one with the strip resistance!");
        }
        // this->layers = std::move(layersToSet);
        this->Rp = Rp;
        this->Rap = Rap;
        // this->fileSave = filename;
        this->MR_mode = CLASSIC;
        // this->layerNo = this->layers.size();
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
    Junction(std::vector<Layer> layersToSet,
             std::string filename,
             std::vector<double> Rx0,
             std::vector<double> Ry0,
             std::vector<double> AMR_X,
             std::vector<double> AMR_Y,
             std::vector<double> SMR_X,
             std::vector<double> SMR_Y,
             std::vector<double> AHE) : Rx0(std::move(Rx0)),
                                        Ry0(std::move(Ry0)),
                                        AMR_X(std::move(AMR_X)),
                                        AMR_Y(std::move(AMR_Y)),
                                        SMR_X(std::move(SMR_Y)),
                                        SMR_Y(std::move(SMR_Y)),
                                        AHE(std::move(AHE))

    {
        if (this->layerNo != Rx0.size() != Ry0.size() != AMR_X.size() != AMR_Y.size() != AHE.size() != SMR_X.size() != SMR_Y.size())
        {
            throw std::invalid_argument("Layers and Rx0, Ry, AMR, AMR and SMR must be of the same size!");
        }
        // this->layers = std::move(layersToSet);
        // this->fileSave = filename;
        this->MR_mode = STRIP;
        // this->layerNo = this->layers.size();
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

    void setLayerAnisotropyDriver(std::string layerID, AxialDriver driver)
    {
        axiallayerSetter(layerID, &Layer::setAnisotropyDriver, driver);
    }
    void setLayerExternalFieldDriver(std::string layerID, AxialDriver driver)
    {
        axiallayerSetter(layerID, &Layer::setExternalFieldDriver, driver);
    }
    void setLayerOerstedFieldDriver(std::string layerID, AxialDriver driver)
    {
        axiallayerSetter(layerID, &Layer::setLayerOerstedFieldDriver, driver);
    }
    void setLayerCurrentDriver(std::string layerID, ScalarDriver driver)
    {
        scalarlayerSetter(layerID, &Layer::setCurrentDriver, driver);
    }
    void setLayerIECDriver(std::string layerID, ScalarDriver driver)
    {
        scalarlayerSetter(layerID, &Layer::setIECDriver, driver);
    }

    void logLayerParams(double t, bool calculateEnergies = false)
    {
        for (Layer &layer : this->layers)
        {
            for (int i = 0; i < 3; i++)
            {
                this->log[layer.id + "_m" + vectorNames[i]].push_back(layer.mag[i]);
                this->log[layer.id + "_Hext" + vectorNames[i]].push_back(layer.H_log[i]);
                this->log[layer.id + "_K" + vectorNames[i]].push_back(layer.K_log[i]);
            }

            if (layer.includeSTT)
                this->log[layer.id + "_I"].push_back(layer.I_log);

            if (calculateEnergies)
            {
                this->log[layer.id + "_EZeeman"].push_back(EnergyDriver::calculateZeemanEnergy(layer.mag,
                                                                                               layer.Hext,
                                                                                               layer.cellVolume,
                                                                                               layer.Ms));
                this->log[layer.id + "_EAnis"].push_back(EnergyDriver::calculateAnisotropyEnergy(layer.mag,
                                                                                                 layer.anisAxis,
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
        if (MR_mode == CLASSIC)
        {
            const auto magnetoresistance = getMagnetoresistance();
            this->log["R_free_bottom"].push_back(magnetoresistance[0]);
        }
        else if (MR_mode == STRIP)
        {
            const auto magnetoresistance = getMagnetoresistance();
            this->log["Rx"].push_back(magnetoresistance[0]);
            this->log["Ry"].push_back(magnetoresistance[1]);
        }
        this->log["time"].push_back(t);

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


    std::map<std::string, std::vector<double>>
    spectralFFT(double minTime = 10.0e-9, double timeStep = 1e-11)
    {

        if (this->log.empty())
        {
            throw std::runtime_error("Empty log! Cannot proceed without running a simulation!");
        }
        auto it = std::find_if(this->log["time"].begin(), this->log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });
        // turn into index
        const int thresIdx = (int)(this->log["time"].end() - it);
        const int cutSize = this->log["time"].size() - thresIdx;

        // plan creation is not thread safe
        const double normalizer = timeStep * cutSize;
        const int maxIt = (cutSize % 2) ? cutSize / 2 : (cutSize - 1) / 2;
        std::vector<double> frequencySteps(maxIt);
        frequencySteps[0] = 0;
        for (int i = 1; i <= maxIt; i++)
        {
            frequencySteps[i - 1] = (i - 1) / normalizer;
        }
        // plan creation is not thread safe
        std::map<std::string, std::vector<double>> spectralFFTResult;
        spectralFFTResult["frequencies"] = std::move(frequencySteps);

        for (const auto &l : this->layers)
        {
            for (const auto &magTag : this->vectorNames)
            {
                std::vector<double> cutMag(this->log[l.id + "_m" + magTag].begin() + thresIdx, this->log[l.id + "_m" + magTag].end());

                // define FFT plan
                fftw_complex out[cutMag.size()];
                fftw_plan plan = fftw_plan_dft_r2c_1d(cutMag.size(),
                                                      cutMag.data(),
                                                      out,
                                                      FFTW_ESTIMATE); // here it's weird, FFT_FORWARD produces an empty plan

                if (plan == NULL)
                {
                    throw std::runtime_error("Plan creation for fftw failed, cannot proceed");
                }
                fftw_execute(plan);
                const int outBins = (cutMag.size() + 1) / 2;
                std::vector<double> amplitudes; //, phases;

                amplitudes.push_back(out[0][0]);
                // phases[0] = 0;
                for (int i = 1; i < outBins; i++)
                {
                    const auto tandem = out[i];
                    const double real = tandem[0];
                    const double img = tandem[1];
                    // phases[i] = tan(img / real);
                    amplitudes.push_back(sqrt(pow(real, 2) + pow(img, 2)));
                }
                spectralFFTResult[l.id + "_m" + magTag + "_amplitude"] = std::move(amplitudes);
                // spectralFFTResult["phase"] = std::move(phases);
                fftw_destroy_plan(plan);
            }
        }
        return spectralFFTResult;
    }

    std::map<std::string, double>
    calculateFFT(double minTime = 10.0e-9, double timeStep = 1e-11)
    {

        if (this->log.empty())
        {
            throw std::runtime_error("Empty log! Cannot proceed without running a simulation!");
        }
        auto it = std::find_if(this->log["time"].begin(), this->log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });
        // turn into index
        const int thresIdx = (int)(this->log["time"].end() - it);
        const int cutSize = this->log["time"].size() - thresIdx;

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
            std::vector<double> cutMag(this->log["free_m" + magTag].begin() + thresIdx, this->log["free_m" + magTag].end());
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
            fftw_destroy_plan(plan);
            result[magTag + "_amplitude"] = amplitudes;
            result[magTag + "_phase"] = phases;
            maxAmpls[magTag + "_resonant"] = maxFreq;
            maxAmpls[magTag + "_amplitude"] = maxAmpl;
            maxAmpls[magTag + "_phase"] = maxPhase;
        }
        return maxAmpls;
    }

    void runSingleLayerRK4Iteration(double t, double timeStep)
    {
        /**
         * Single layer iteration. IEC interaction is turned off.
         * @param t: current time
         * @param timeStep: integration step
         * */
        this->layers[0].rk4_step(
            t, timeStep, CVector(0., 0., 0.));
    }

    void runMultiLayerRK4Iteration(double t, double timeStep)
    {
        CVector l1mag = this->layers[0].mag;
        CVector l2mag = this->layers[1].mag;
        this->layers[0].rk4_step(
            t, timeStep, l2mag);
        this->layers[1].rk4_step(
            t, timeStep, l1mag);
    }

    std::vector<double> stripMagnetoResistance(std::vector<double> Rx0,
                                               std::vector<double> Ry0,
                                               std::vector<double> AMR_X,
                                               std::vector<double> SMR_X,
                                               std::vector<double> AMR_Y,
                                               std::vector<double> SMR_Y,
                                               std::vector<double> AHE)
    {
        double Rx_acc = 0.0;
        double Ry_acc = 0.0;

        for (int i = 0; i > this->layers.size(); i++)
        {
            const double Rx = Rx0[i] + AMR_X[i] * pow(this->layers[i].mag.x, 2) + SMR_X[i] * pow(this->layers[i].mag.y, 2);
            const double Ry = Ry0[i] + 0.5 * AHE[i] * this->layers[i].mag.z + \
                        (AMR_Y[i] - SMR_Y[i]) * this->layers[i].mag.x * this->layers[i].mag.y;
            Rx_acc += Rx;
            Ry_acc += Ry;
        }

        return {1. / Rx_acc, 1. / Ry_acc};
    }

    double calculateMagnetoresistance(double cosTheta)
    {
        return this->Rp + (((this->Rap - this->Rp) / 2.0) * (1.0 - cosTheta));
    }

    std::vector<double> getMagnetoresistance()
    {
        if (this->MR_mode == CLASSIC)
        {
            return {calculateMagnetoresistance(c_dot(layers[0].mag, layers[1].mag))};
        }
        else
        {
            return stripMagnetoResistance(this->Rx0, this->Ry0,
                                          this->AMR_X,
                                          this->SMR_X,
                                          this->AMR_Y,
                                          this->SMR_Y,
                                          this->AHE);
        }
    }

    void runSimulation(double totalTime, double timeStep = 1e-13, double writeFrequency = 1e-11,
                       bool persist = true, bool log = false, bool calculateEnergies = false)
    {

        const unsigned int totalIterations = (int)(totalTime / timeStep);
        double t;
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