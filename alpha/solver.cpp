#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <fstream>
#include <string>

#include "cvector.hpp"

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 2.21e5
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

    CVector Hconst, anis, mag;

    double K, J, Ms, thickness;
    double Kvar, Jvar, Hvar;
    double K_frequency, J_frequency, H_frequency = 0.0;
    double J_log, K_log, H_log;

    Axis Kax, Hax;

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
                                                dipole_tensor(dipole_tensor)
    {
    }

    double sinusoidalUpdate(double amplitude, double frequency, double time, double phase)
    {
        return amplitude * sin(2 * M_PI * time * frequency + phase);
    }
    double updateCoupling(double J, double frequency, double time, double phase)
    {
        return sinusoidalUpdate(J, frequency, time, phase);
    }
    CVector updateAxial(double H, double frequency, double time, double phase, Axis axis)
    {
        CVector *result = new CVector();
        switch (axis)
        {
        case xaxis:
            result->x = sinusoidalUpdate(H, frequency, time, phase);
        case yaxis:
            result->y = sinusoidalUpdate(H, frequency, time, phase);
        case zaxis:
            result->z = sinusoidalUpdate(H, frequency, time, phase);
        }
        return *result;
    }

    CVector Heff(double time, CVector otherMag, double otherMs)
    {
        double kVar, couplingVar = 0.0;
        CVector Heff = {0., 0., 0.};
        // kVar = updateAxial(this->K, this->K_frequency, time, 0, this->Kax);
        // HextVar = updateExternalField(time);
        K_log = K + kVar;

        // couplingVar = updateCoupling(time);

        // J_log = J + updateCoupling();

        Heff = Hconst + calculateAnisotropy() +
               calculateIEC(otherMag) +
               // demag
               calculate_tensor_interaction(this->mag, this->demag_tensor, otherMs) +
               calculate_tensor_interaction(this->mag, this->dipole_tensor, this->Ms);

        return Heff;
    }
    CVector calculateAnisotropy()
    {
        double nom = (2 * K_log) * c_dot(anis, mag) / (MAGNETIC_PERMEABILITY * Ms);
        return anis * nom;
    }

    CVector calculateIEC(CVector coupledMag)
    {
        double nom = this->J / (MAGNETIC_PERMEABILITY * this->Ms * this->thickness);
        return (coupledMag - this->mag) * nom;
    }

    CVector llg(double time, CVector m, CVector coupledMag, double otherMs)
    {
        CVector heff, prod, prod2, dmdt;
        heff = Heff(time, coupledMag, otherMs);
        prod = c_cross(m, heff);
        prod2 = c_cross(m, prod);
        dmdt = prod * -GYRO - c_cross(m, prod2) * GYRO * DAMPING;
        return dmdt;
    }

    void setGlobalExtenralFieldValue(CVector Hval)
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
        CVector k1, k2, k3, k4, m_t;
        m_t = mag;
        k1 = llg(time, m_t, coupledMag, otherMs) * time_step;
        k2 = llg(time + 0.5 * time_step, m_t + k1 * 0.5, coupledMag, otherMs) * time_step;
        k3 = llg(time + 0.5 * time_step, m_t + k2 * 0.5, coupledMag, otherMs) * time_step;
        k4 = llg(time + time_step, m_t + k3, coupledMag, otherMs) * time_step;
        m_t = m_t + (k1 + k2 * 2.0 + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        mag = m_t;
    }
};

class Junction
{
public:
    std::vector<Layer> layers;
    double Rp, Rap = 0.0;

    Junction(std::vector<Layer> layersToSet)
    {
        this->layers = std::move(layersToSet);
        Rp = 100;
        Rap = 200;
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
    // double coupleLayers(std::string l1ID, std::string l2ID, double couplingStrength)
    // {
    //     // Layer &l1 = findLayerByID(l1ID);
    //     // Layer &l2 = findLayerByID(l2ID);
    //     // l1.otherMs = l2.Ms;
    //     // l2.otherMs = l1.Ms;
    // }

    double calculateMagnetoresistance(double cosTheta)
    {
        return this->Rp + ((this->Rp + this->Rap) / 2.0 * (1.0 - cosTheta));
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

    void runSimulation(double totalTime, double timeStep)
    {

        int totalIterations = (int)(totalTime / timeStep);
        double t;
        std::vector<double> magnetoresistance;
        for (int i = 0; i < totalIterations; i++)
        {
            t = i * timeStep;

            CVector l1mag = this->layers[0].mag;
            CVector l2mag = this->layers[1].mag;
            layers[0].rk4_step(
                t, timeStep, l1mag, layers[1].Ms);
            layers[1].rk4_step(
                t, timeStep, l2mag, layers[0].Ms);

            for (Layer &l1 : layers)
            {
                for (Layer &l2 : layers)
                {
                    if (l1.id != l2.id)
                    {
                        double cosTheta = c_dot(l1.mag, l2.mag);
                        magnetoresistance.push_back(
                            calculateMagnetoresistance(cosTheta));
                    }
                }
            }
        }
    }
};

int main(void)
{
    return 0;
}