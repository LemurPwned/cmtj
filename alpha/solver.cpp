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

#include "cvector.hpp"

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 221000.0
#define DAMPING 0.011
#define TtoAm 795774.715459
#define HBAR 6.62607015e-34 / (2 * M_PI)
#define ELECTRON_CHARGE 1.60217662e-19

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

class Layer
{

public:
    std::string id;

    CVector H_log, Hconst, anis, mag = {0., 0., 0.};

    double K, J, Ms, thickness = 0.0;
    double Kvar, Jvar, Hvar = 0.0;
    double K_frequency, J_frequency, H_frequency = 0.0;
    double J_log, K_log = 0.0;

    Axis Hax = xaxis;

    std::vector<CVector> demag_tensor, dipole_tensor;
    Layer(std::string id,
          CVector mag,
          CVector anis,
          double K,
          double Ms,
          double J,
          double thickness,
          std::vector<CVector> demag_tensor,
          std::vector<CVector> dipole_tensor) : id(id), mag(mag), anis(anis), K(K), Ms(Ms), J(J),
                                                thickness(thickness), demag_tensor(demag_tensor),
                                                dipole_tensor(dipole_tensor) {}

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
        switch (axis)
        {
        case xaxis:
            result->x = sinusoidalUpdate(amplitude, frequency, time, phase);
        case yaxis:
            result->y = sinusoidalUpdate(amplitude, frequency, time, phase);
        case zaxis:
            result->z = sinusoidalUpdate(amplitude, frequency, time, phase);
        }
        return *result;
    }

    CVector calculateHeff(double time, CVector otherMag, double otherMs)
    {
        CVector Heff = {0., 0., 0.};

        Heff = calculateExternalField(time) +
               calculateAnisotropy(time) +
               calculateIEC(time, otherMag) +
               // demag
              // check the interaction here to be sure
               calculate_tensor_interaction(otherMag, this->demag_tensor, this->Ms) +
               // dipole
               calculate_tensor_interaction(this->mag, this->dipole_tensor, this->Ms);

        return Heff;
    }

    CVector calculateExternalField(double time)
    {
        this->H_log = this->Hconst + updateAxial(this->Hvar, this->H_frequency, time, 0, this->Hax);
        return this->H_log;
    }

    CVector calculateAnisotropy(double time)
    {
        this->K_log = this->K + sinusoidalUpdate(this->Kvar, this->K_frequency, time, 0);
        double nom = (2 * this->K_log) * c_dot(this->anis, this->mag) / (MAGNETIC_PERMEABILITY * this->Ms);
        return this->anis * nom;
    }

    CVector calculateIEC(double time, CVector coupledMag)
    {
        this->J_log = this->J + sinusoidalUpdate(this->Jvar, this->J_frequency, time, 0);
        double nom = this->J_log / (MAGNETIC_PERMEABILITY * this->Ms * this->thickness);
        return (coupledMag - this->mag) * nom;
    }

    CVector llg(double time, CVector m, CVector coupledMag, CVector heff, double otherMs)
    {
        CVector prod, prod2, dmdt;
        // heff = calculateHeff(time, coupledMag, otherMs);
        prod = c_cross(m, heff);
        prod2 = c_cross(m, prod);
        dmdt = prod * -GYRO - prod2 * GYRO * DAMPING;
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

    void rk4_step(double time, double time_step, CVector coupledMag, double otherMs)
    {
        CVector k1, k2, k3, k4, m_t, heff;
        m_t = mag;
        heff = calculateHeff(time, coupledMag, otherMs);
        k1 = llg(time, m_t, coupledMag, heff, otherMs) * time_step;
        k2 = llg(time + 0.5 * time_step, m_t + k1 * 0.5, coupledMag, heff, otherMs) * time_step;
        k3 = llg(time + 0.5 * time_step, m_t + k2 * 0.5, coupledMag, heff, otherMs) * time_step;
        k4 = llg(time + time_step, m_t + k3, coupledMag, heff, otherMs) * time_step;
        m_t = m_t + (k1 + (k2 * 2.0) + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        mag = m_t;
    }
};

class Junction
{
public:
    std::vector<Layer> layers;
    double Rp, Rap = 0.0;
    std::map<std::string, std::vector<double>> log;
    std::string fileSave;
    unsigned int logLength = 0;
    std::vector<std::string> vectorNames = {"x", "y", "z"};

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
        CVector fieldToSet(0,0,0);
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
            this->log["L1Hext" + vectorNames[i]].push_back(this->layers[0].Hconst[i]);
            this->log["L2Hext" + vectorNames[i]].push_back(this->layers[1].Hconst[i]);
        }
        this->log["L1K"].push_back(this->layers[0].K_log);
        this->log["L2K"].push_back(this->layers[1].K_log);
        this->log["R_free_bottom"].push_back(magnetoresistance);
        this->log["time"].push_back(t);

        logLength++;
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
    }

    double calculateVoltageSpinDiode(double frequency, double power = 10e-6, const double minTime = 10e-9)
    {

        double omega = 2 * M_PI * frequency;
        std::string res = "R_free_bottom";
        std::vector<double> &resistance = this->log[res];
        auto it = std::find_if(this->log["time"].begin(), this->log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });
        // turn into index
        int thresIdx = (int)(this->log["time"].end() - it);
        int cutSize = this->log["time"].size() - thresIdx;
        // Rpp
        double RppMax = *std::max_element(resistance.begin() + thresIdx, resistance.end());
        double RppMin = *std::min_element(resistance.begin() + thresIdx, resistance.end());
        double avgR = std::accumulate(resistance.begin() + thresIdx, resistance.end(), 0.0) / cutSize;
        double Iampl = sqrt(power / avgR);
        std::vector<double> voltage, current;
        std::transform(
            this->log["time"].begin() + thresIdx, this->log["time"].end(),
            std::back_inserter(current), 
            [&Iampl, &omega](const double &time) { return Iampl * sin(omega * time); });

        for (unsigned int i = 0; i < cutSize; i++)
        {
            voltage.push_back(resistance[thresIdx + i] * current[i]);
        }
        double Vmix = std::accumulate(voltage.begin(), voltage.end(), 0.0) / voltage.size();
        // std::cout << "Rpp: " << RppMax - RppMin << ", meanR: " << avgR << ", VSD: " << Vmix << std::endl;
        return Vmix;
    }

    void runSimulation(double totalTime, double timeStep = 1e-13)
    {

        unsigned int totalIterations = (int)(totalTime / timeStep);
        double t;
        std::vector<double> magnetoresistance;
        unsigned int writeEvery = (int)(0.01 * 1e-9 / timeStep) - 1;

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        for (unsigned int i = 0; i < totalIterations; i++)
        {
            t = i * timeStep;

            CVector l1mag = this->layers[0].mag;
            CVector l2mag = this->layers[1].mag;
            layers[0].rk4_step(
                t, timeStep, l2mag, layers[1].Ms);
            layers[1].rk4_step(
                t, timeStep, l1mag, layers[0].Ms);

            if (!(i % writeEvery))
            {
                double magRes = calculateMagnetoresistance(c_dot(layers[0].mag, layers[1].mag));

                logLayerParams(t, magRes);
            }
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        // std::cout << "Simulation time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    }
};

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

    Layer l1("free",
             CVector(0., 0., 1.),
             CVector(0, 0., 1.),
             900e3,
             1200e3,
             0.0, 1.4e-9, demagTensor, dipoleTensor);
    Layer l2("bottom",
             CVector(0., 0., 1.),
             CVector(0, 0., 1.),
             1000e3,
             1000e3,
             0.0, 7e-10, demagTensor, dipoleTensor);

    Junction mtj(
        {l1, l2}, "test.csv");



    Layer l3("free",
             CVector(0., 0., 1.),
             CVector(0, 0., 1.),
             900e3,
             1200e3,
             -3e-6, 1.4e-9, demagTensor, dipoleTensor);
    Layer l4("bottom",
             CVector(0., 0., 1.),
             CVector(0, 0., 1.),
             10000e3,
             1000e3,
             -3e-6, 7e-10, demagTensor, dipoleTensor);
    Junction mtj2(
        {l3, l4}, "test.csv");

    double minField = 000.0;
    double maxField = 400.0;
    int numPoints = 50;
    double spacing = (maxField - minField) / numPoints;
    std::cout << spacing << std::endl;
    std::ofstream vsdFile;
    vsdFile.open("VSD-anisotropy.csv");
    vsdFile << "H;Vmix\n";
    for (double field = minField; field <= maxField; field += spacing)
    {
        mtj.setConstantExternalField((field/1000) * TtoAm, xaxis);
        mtj.setLayerAnisotropyUpdate("free", 900, 7e9, 0);
        mtj.setLayerAnisotropyUpdate("bottom", 900, 7e9, 0);
        mtj.runSimulation(15e-9);
        double res = mtj.calculateVoltageSpinDiode(7e9);

        // clear logs
        mtj.log.clear();
        vsdFile << field << ";" << res << "\n";
    }
    vsdFile.close();

    std::ofstream vsdFile2;
    std::cout << "Finished ANIS\n"
              << std::endl;
    vsdFile2.open("VSD-IEC.csv");
    vsdFile2 << "H;Vmix\n";
    for (double field = minField; field <= maxField; field += spacing)
    {
        mtj2.setConstantExternalField((field/1000) * TtoAm, xaxis);
        mtj2.setLayerIECUpdate("free", 1e-6, 7e9, 0);
        mtj2.setLayerIECUpdate("bottom", 1e-6, 7e9, 0);
        mtj2.runSimulation(15e-9);
        double res = mtj2.calculateVoltageSpinDiode(7e9);

        // clear logs
        mtj2.log.clear();
        vsdFile2 << field << ";" << res << "\n";
    }
    vsdFile2.close();
    return 0;
}