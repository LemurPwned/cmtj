#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>
#include <cstring>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>
#include <numeric>
#include <thread>
#include <future>
#include <tuple>
#include <random>
#include <complex>
#include <fftw3.h>
#include "cvector.hpp"

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 221000.0
#define DAMPING 0.011
#define TtoAm 795774.715459
#define HBAR 6.62607015e-34 / (2 * M_PI)
#define ELECTRON_CHARGE 1.60217662e-19
#define BOLTZMANN_CONST 1.380649e-23
#define SPIN_POLARISATION_EFF 0.8
#define CURRENT_DENSITY 1
#define B1 0.5
#define B2 0.3

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

enum Axis
{
    xaxis,
    yaxis,
    zaxis
};

enum ExcitationMode
{
    axial,
    step
};

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);
class Layer
{
private:
    double cellVolume, cellSurface = 0.0;

public:
    std::string id;

    CVector H_log, Hconst, anis, mag;

    double K = 0.0, J = 0.0, Ms = 0.0;
    double Kvar, Jvar, Hvar = 0.0;
    double K_frequency = 0.0, J_frequency = 0.0, H_frequency = 0.0;
    double J_log = 0.0, K_log = 0.0;

    // resting temperature in Kelvin
    double temperature = 0.0;
    double thickness = 0.0;
    Axis Hax = xaxis;

    double Hstart = 0.0, Hstop = 0.0, Hstep = 0.0;

    bool includeSTT = false;

    std::vector<CVector> demag_tensor, dipole_tensor;
    Layer(std::string id,
          CVector mag,
          CVector anis,
          double K,
          double Ms,
          double J,
          double thickness,
          double cellSurface,
          std::vector<CVector> demag_tensor,
          std::vector<CVector> dipole_tensor, double temperature = 0.0) : id(id),
                                                                          mag(mag), anis(anis), K(K), Ms(Ms), J(J),
                                                                          thickness(thickness),
                                                                          cellSurface(cellSurface),
                                                                          demag_tensor(demag_tensor),
                                                                          dipole_tensor(dipole_tensor),
                                                                          temperature(temperature)

    {
        this->cellVolume = this->cellSurface * this->thickness;
    }

    double sinusoidalUpdate(double amplitude, double frequency, double time, double phase = 0.0)
    {
        return amplitude * sin(2 * M_PI * time * frequency + phase);
    }

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

    CVector updateAxial(double amplitude, double frequency, double time, double phase, Axis axis)
    {
        CVector *result = new CVector();
        const double updateValue = sinusoidalUpdate(amplitude, frequency, time, phase);
        switch (axis)
        {
        case xaxis:
            result->x = updateValue;
        case yaxis:
            result->y = updateValue;
        case zaxis:
            result->z = updateValue;
        }
        return *result;
    }

    CVector updateAxialStep(double amplitude, double time, double timeStart, double timeStop, Axis axis)
    {
        CVector *result = new CVector();
        const double updateValue = stepUpdate(amplitude, time, timeStart, timeStop);
        switch (axis)
        {
        case xaxis:
            result->x = updateValue;
        case yaxis:
            result->y = updateValue;
        case zaxis:
            result->z = updateValue;
        }
        return *result;
    }

    CVector calculateHeff(double time, double timeStep, CVector otherMag)
    {
        CVector Heff = {0., 0., 0.};

        Heff = calculateExternalField(time) +
               calculateAnisotropy(time) +
               calculateIEC(time, otherMag) +
               // demag
               // check the interaction here to be sure
               calculate_tensor_interaction(otherMag, this->demag_tensor, this->Ms) +
               // dipole
               calculate_tensor_interaction(this->mag, this->dipole_tensor, this->Ms) +
               calculateStochasticThermalField(timeStep);

        return Heff;
    }

    CVector calculateStochasticThermalField(double timeStep)
    {
        CVector res(distribution, generator);
        const double nom = sqrt((2 * DAMPING * BOLTZMANN_CONST * this->temperature) /
                                (MAGNETIC_PERMEABILITY * GYRO * this->cellVolume * this->Ms * timeStep));
        return res * nom;
    }

    CVector calculateExternalField(double time)
    {

        this->H_log = this->Hconst + updateAxial(this->Hvar, this->H_frequency, time, 0, this->Hax);
        this->H_log += updateAxialStep(this->Hstep, time, this->Hstart, this->Hstop, this->Hax);
        return this->H_log;
    }

    CVector calculateAnisotropy(double time)
    {
        this->K_log = this->K + sinusoidalUpdate(this->Kvar, this->K_frequency, time, 0);
        const double nom = (2 * this->K_log) * c_dot(this->anis, this->mag) / (MAGNETIC_PERMEABILITY * this->Ms);
        return this->anis * nom;
    }

    CVector calculateIEC(double time, CVector coupledMag)
    {
        this->J_log = this->J + sinusoidalUpdate(this->Jvar, this->J_frequency, time, 0);
        const double nom = this->J_log / (MAGNETIC_PERMEABILITY * this->Ms * this->thickness);
        return (coupledMag - this->mag) * nom;
    }

    CVector llg(double time, CVector m, CVector coupledMag, CVector heff)
    {
        CVector prod, prod2, dmdt;
        // heff = calculateHeff(time, coupledMag, otherMs);
        prod = c_cross(m, heff);
        prod2 = c_cross(m, prod);
        dmdt = prod * -GYRO - prod2 * GYRO * DAMPING;
        if (this->includeSTT)
        {
            CVector prod3;
            // damping-like torque factor
            const double aJ = HBAR * SPIN_POLARISATION_EFF * CURRENT_DENSITY /
                              (2 * ELECTRON_CHARGE * MAGNETIC_PERMEABILITY * this->Ms * this->thickness);

            const double bJ = B1 * CURRENT_DENSITY + B2 * CURRENT_DENSITY * CURRENT_DENSITY;

            prod3 = c_cross(m, coupledMag);
            dmdt += c_cross(m, prod3) * GYRO * aJ;
            dmdt += prod3 * GYRO * bJ;
        }
        return dmdt;
    }

    void setGlobalExternalFieldValue(CVector &Hval)
    {
        this->Hconst = Hval;
    }

    void setCoupling(double amplitude)
    {
        this->J = amplitude;
    }
    void setAnisotropy(double amplitude)
    {
        this->K = amplitude;
    }

    void rk4_step(double time, double timeStep, CVector coupledMag)
    {
        CVector k1, k2, k3, k4, m_t, heff;
        m_t = mag;
        heff = calculateHeff(time, timeStep, coupledMag);
        k1 = llg(time, m_t, coupledMag, heff) * timeStep;
        k2 = llg(time + 0.5 * timeStep, m_t + k1 * 0.5, coupledMag, heff) * timeStep;
        k3 = llg(time + 0.5 * timeStep, m_t + k2 * 0.5, coupledMag, heff) * timeStep;
        k4 = llg(time + timeStep, m_t + k3, coupledMag, heff) * timeStep;
        m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        mag = m_t;
    }
};

class Junction
{
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
        return this->Rp + (((this->Rp + this->Rap) / 2.0) * (1.0 - cosTheta));
    }

    void setConstantExternalField(double Hval, Axis axis)
    {
        CVector fieldToSet(0, 0, 0);
        switch (axis)
        {
        case xaxis:
            fieldToSet.x = Hval;
            break;
        case yaxis:
            fieldToSet.y = Hval;
            break;
        case zaxis:
            fieldToSet.z = Hval;
            break;
        }
        for (Layer &l : this->layers)
        {
            l.setGlobalExternalFieldValue(fieldToSet);
        }
    }

    void setLayerAnisotropyUpdate(std::string layerID, double amplitude, double frequency, double phase)
    {
        Layer &l1 = findLayerByID(layerID);
        l1.Kvar = amplitude;
        l1.K_frequency = frequency;
    }
    void setLayerIECUpdate(std::string layerID, double amplitude, double frequency, double phase)
    {
        Layer &l1 = findLayerByID(layerID);
        l1.Jvar = amplitude;
        l1.J_frequency = frequency;
    }

    void setLayerStepUpdate(std::string layerID, double Hstep, double timeStart, double timeStop, Axis hax)
    {
        Layer &l1 = findLayerByID(layerID);
        l1.Hax = hax;
        l1.Hstep = Hstep;
        l1.Hstart = timeStart;
        l1.Hstop = timeStop;
    }

    void setLayerCoupling(std::string layerID, double J)
    {
        bool found = false;
        for (Layer l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                l.setAnisotropy(J);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }

    void setLayerAnisotropy(std::string layerID, double K)
    {
        bool found = false;
        for (Layer l : this->layers)
        {
            if (l.id == layerID || layerID == "all")
            {
                l.setAnisotropy(K);
                found = true;
            }
        }
        if (!found)
        {
            throw std::runtime_error("Failed to find a layer with a given id!");
        }
    }

    void logLayerParams(double t, double magnetoresistance)
    {
        for (int i = 0; i < 3; i++)
        {
            this->log["L1m" + vectorNames[i]].push_back(this->layers[0].mag[i]);
            this->log["L2m" + vectorNames[i]].push_back(this->layers[1].mag[i]);
            this->log["L1Hext" + vectorNames[i]].push_back(this->layers[0].H_log[i]);
            this->log["L2Hext" + vectorNames[i]].push_back(this->layers[1].H_log[i]);
        }
        this->log["L1K"].push_back(this->layers[0].K_log);
        this->log["L2K"].push_back(this->layers[1].K_log);
        this->log["R_free_bottom"].push_back(magnetoresistance);
        this->log["time"].push_back(t);

        this->logLength++;
    }

    void saveLogs()
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

    double calculateVoltageSpinDiode(double frequency, double power = 10e-6, const double minTime = 10e-9)
    {
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
        return Vmix;
    }

    std::map<std::string, std::vector<double>> calculateFFT(double minTime = 10.0e-9, double timeStep = 1e-11)
    {
        std::cout << "FFT calculation" << std::endl;
        auto it = std::find_if(this->log["time"].begin(), this->log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });
        // turn into index
        const int thresIdx = (int)(this->log["time"].end() - it);
        const int cutSize = this->log["time"].size() - thresIdx;

        std::cout << "Initiating FFT plan execution " << std::endl;

        // plan creation is not thread safe
        const double normalizer = timeStep * cutSize;
        const int maxIt = (cutSize % 2) ? cutSize / 2 : (cutSize - 1) / 2;
        std::cout << maxIt << "," << cutSize << std::endl;
        std::vector<double> frequencySteps(maxIt);
        frequencySteps[0] = 0;
        for (int i = 1; i <= maxIt; i++)
        {
            frequencySteps[i - 1] = (i - 1) / normalizer;
        }

        std::map<std::string, std::vector<double>> result;
        for (const auto &magTag : this->vectorNames)
        {
            std::cout << "Doing FFT for: " << magTag << std::endl;
            std::vector<double> cutMag(this->log["L1m" + magTag].begin() + thresIdx, this->log["L1m" + magTag].end());
            std::cout << "Sub size: " << cutMag.size() << std::endl;
            fftw_complex out[cutMag.size()];
            // define FFT plan
            fftw_plan plan = fftw_plan_dft_r2c_1d(cutMag.size(),
                                                  cutMag.data(),
                                                  out,
                                                  FFTW_ESTIMATE);
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
                }
            }
            std::cout << "Max frequency for the system is " << maxFreq << " with " << maxAmpl << std::endl;
            fftw_destroy_plan(plan);
            result[magTag + "_amplitude"] = amplitudes;
            result[magTag + "_phase"] = phases;
        }
        return result;
    }

    void runSimulation(double totalTime, double timeStep = 1e-13, bool persist = false, bool log = false)
    {

        const unsigned int totalIterations = (int)(totalTime / timeStep);
        double t;
        const unsigned int writeEvery = (int)(0.01 * 1e-9 / timeStep) - 1;
        std::vector<double> magnetoresistance;

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (unsigned int i = 0; i < totalIterations; i++)
        {
            t = i * timeStep;

            CVector l1mag = this->layers[0].mag;
            CVector l2mag = this->layers[1].mag;
            layers[0].rk4_step(
                t, timeStep, l2mag);
            layers[1].rk4_step(
                t, timeStep, l1mag);

            if (!(i % writeEvery))
            {
                const double magRes = calculateMagnetoresistance(c_dot(layers[0].mag, layers[1].mag));
                logLayerParams(t, magRes);
            }
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if (persist)
            saveLogs();
        if (log)
            std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    }
};

void customResultMap(std::map<std::string, std::vector<double>> resultMap, std::string filename)
{

    std::ofstream saveFile;

    saveFile.open(filename);
    int logLength = 0;
    for (const auto [key, _] : resultMap)
    {
        if (!logLength)
            logLength = resultMap[key].size();
        saveFile << key << ";";
    }
    saveFile << "\n";
    for (unsigned int i = 0; i < logLength; i++)
    {
        for (const auto [_, value] : resultMap)
        {
            saveFile << value[i] << ";";
        }
        saveFile << "\n";
    }
    saveFile.close();
}

void threadedSimulation(Junction cjx, double minField, double maxField, int numberOfPoints, std::ofstream &vsdFile)
{
    const int threadNum = std::thread::hardware_concurrency() - 2;
    std::vector<std::future<std::vector<std::tuple<double, double>>>> threadResults;
    threadResults.reserve(threadNum);

    const int pointsPerThread = numberOfPoints / threadNum + 1;
    const double spacing = (maxField - minField) / numberOfPoints;
    for (int i = 0; i < threadNum; i++)
    {
        const double threadMinField = pointsPerThread * i * spacing;
        const double threadMaxField = pointsPerThread * (i + 1) * spacing;
        threadResults.emplace_back(std::async([cjx, threadMinField, threadMaxField, spacing]() mutable {
            std::vector<std::tuple<double, double>>
                resAcc;
            const double freq = 7e9;
            for (double field = threadMinField; field <= threadMaxField;
                 field += spacing)
            {
                cjx.setConstantExternalField((field / 1000) * TtoAm, xaxis);
                cjx.setLayerAnisotropyUpdate("free", 12000, freq, 0);
                cjx.setLayerAnisotropyUpdate("bottom", 12000, freq, 0);
                cjx.runSimulation(20e-9);
                const auto vsd = cjx.calculateVoltageSpinDiode(freq);
                resAcc.push_back({field, vsd});
                cjx.log.clear();
            }

            return resAcc;
        }));
    }

    vsdFile << "H;Vmix\n";
    for (auto &result : threadResults)
    {
        for (const auto [field, vsdVal] : result.get())
        {
            vsdFile << field << ";" << vsdVal << "\n";
        }
    };
}

void threadedSimulation2(Junction mtj, double minField, double maxField, int numberOfPoints, std::ofstream &vsdFile,
                         std::tuple<double, double> runnableFunction(Junction mtj, double scanningParam))
{
    const int threadNum = std::thread::hardware_concurrency() - 2;
    std::vector<std::future<std::vector<std::tuple<double, double>>>> threadResults;
    threadResults.reserve(threadNum);

    const int pointsPerThread = numberOfPoints / threadNum + 1;
    const double spacing = (maxField - minField) / numberOfPoints;
    for (int i = 0; i < threadNum; i++)
    {
        const double threadMinField = pointsPerThread * i * spacing;
        const double threadMaxField = pointsPerThread * (i + 1) * spacing;
        threadResults.emplace_back(std::async([mtj, runnableFunction, threadMinField, threadMaxField, spacing]() mutable {
            std::vector<std::tuple<double, double>>
                resAcc;
            for (double field = threadMinField; field <= threadMaxField;
                 field += spacing)
            {
                auto subRes = runnableFunction(mtj, field);
                resAcc.push_back(subRes);
            }
            return resAcc;
        }));
    }

    // vsdFile.open("VSD-anisotropy.csv");
    vsdFile << "H;Vmix\n";
    for (auto &result : threadResults)
    {
        for (const auto [field, vsdVal] : result.get())
        {
            vsdFile << field << ";" << vsdVal << "\n";
        }
    };
}

int main(void)
{

    std::vector<CVector> dipoleTensor = {
        {6.8353909454237E-4, 0., 0.},
        {0., 0.00150694452305927, 0.},
        {0., 0., 0.99780951638608}};
    std::vector<CVector> demagTensor = {
        {5.57049776248663E-4, 0., 0.},
        {0., 0.00125355500286346, 0.},
        {0., 0.0, -0.00181060482770131}};

    Layer l1("free",              // id
             CVector(0., 0., 1.), // mag
             CVector(0, 0., 1.),  // anis
             900e3,               // K
             1200e3,              // Ms
             0.0,                 // J
             1.4e-9,              // thickness
             7e-10 * 7e-10,       // surface
             demagTensor,         // demag
             dipoleTensor);
    Layer l2("bottom",            // id
             CVector(0., 0., 1.), // mag
             CVector(0, 0., 1.),  // anis
             10000e3,             // K
             1000e3,              // Ms
             0.0,                 // J
             7e-10,               //thickness
             7e-10 * 7e-10,       // surface
             demagTensor,         // demag
             dipoleTensor);

    Junction mtj(
        {l1, l2}, "test2.csv");

    double minField = 000.0;
    double maxField = 400.0;
    int numPoints = 50;
    double spacing = (maxField - minField) / numPoints;
    std::cout << spacing << std::endl;
    std::ofstream vsdFile;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    vsdFile.open("VSD-IEC.csv");
    threadedSimulation(mtj, minField, maxField, numPoints, vsdFile);
    vsdFile.close();
    vsdFile.open("VSD-anisotropy.csv");
    threadedSimulation2(mtj, minField, maxField, numPoints, vsdFile,
                        [](Junction mtj, double field) mutable {
                            const double freq = 7e9;
                            mtj.setConstantExternalField((field / 1000) * TtoAm, xaxis);
                            mtj.setLayerAnisotropyUpdate("free", 12000, freq, 0);
                            mtj.setLayerAnisotropyUpdate("bottom", 12000, freq, 0);
                            mtj.runSimulation(20e-9);
                            const auto vsd = mtj.calculateVoltageSpinDiode(freq);
                            mtj.log.clear();
                            return std::tuple<double, double>{field, vsd};
                        });
    vsdFile.close();
    // mtj.setConstantExternalField(250_mT * TtoAm, xaxis);
    // mtj.setLayerStepUpdate("free", 10e-3 * TtoAm, 5_ns, 5.01_ns, xaxis);
    // mtj.setLayerStepUpdate("bottom", 10e-3 * TtoAm, 5_ns, 5.01_ns, xaxis);
    // mtj.runSimulation(20_ns, 1e-13, true);
    // mtj.calculateFFT(10_ns, 1e-11);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    return 0;
}