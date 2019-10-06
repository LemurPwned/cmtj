#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <fstream>
#include <string>

#include <math.h>

#define MAGNETIC_PERMEABILITY 12.57e-7
#define GYRO 2.21e5
#define DAMPING 0.011
#define TtoAm 795774.715459
#define HBAR 6.62607015e-34 / (2 * PI)
#define ELECTRON_CHARGE 1.60217662e-19

#ifndef _CVECTOR_ // Check if C Vector is defined
#define _CVECTOR_ // the next time CVector will be define
#define ZEROVECTOR CVector()

class CVector
{ // Private
public:
    double x, y, z;                  // Order list of 4 elements |x|y|z|w|
    static char sBuffer[38];         // holds the string of given vector
    char *toString();                // Return sBuffer with (x,y,z,w) values
    CVector(void);                   // zero Vector Constructor
    CVector(double, double, double); // Constructor
    CVector(CVector &);              // Copy Vector Constructor
    ~CVector();
    CVector operator+(CVector v)
    {
        return CVector(
            x + v.x,
            y + v.y,
            z + v.z);
    };  // Destructor
        // Addition
    CVector operator-(CVector v)
    {
        return CVector(
            x - v.x,
            y - v.y,
            z - v.z);
    };
    double dotProduct(CVector &);
    // CVector operator=(CVector &v){
    //     return CVector(x + v.x, y+ v.y, z + v.z);
    // }
    void operator=(CVector v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    }                               // Copy Vector
    int operator==(CVector &);      // Comparison Vector
    CVector operator*(CVector v){}; // Cross Product
    CVector operator*(double val)
    {
        return CVector(
            x * val,
            y * val,
            z * val);
    };
    CVector operator/(double val)
    {
        return CVector(
            x / val,
            y / val,
            z / val);
    };
    double operator[](int i)
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            return z;
    } // Scalar Multiplication
    double length()
    {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }; // Magnitude
    void normalize()
    {
        double mag = this->length();
        x = x / mag;
        y = y / mag;
        z = z / mag;
    };
};
#endif // _CVECTOR_
// class Layer:
//     def __init__(self, id_, start_mag, start_anisotropy, K, Ms, coupling, thickness,
//                     demagnetisation_tensor, dipole_tensor):

CVector calculate_tensor_interaction(CVector m,
                                     std::vector<CVector> tensor,
                                     double Ms)
{
    return CVector(
        -Ms * tensor[0][0] * m[0] - Ms * tensor[0][1] * m[1] - Ms * tensor[0][2] * m[2],
        -Ms * tensor[1][0] * m[0] - Ms * tensor[1][1] * m[1] - Ms * tensor[1][2] * m[2],
        -Ms * tensor[2][0] * m[0] - Ms * tensor[2][1] * m[1] - Ms * tensor[2][2] * m[2]);
}

CVector c_cross(CVector a, CVector b)
{
    return CVector(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);
}

double c_dot(CVector a, CVector b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Layer
{

public:
    int id;
    CVector mag, anis, Hext;
    double K, Ms, coupling, thickness;
    double coupling_log, K_log;
    std::vector<CVector> demag_tensor, dipole_tensor;
    Layer(int id,
          CVector mag,
          CVector anis,
          double K,
          double Ms,
          double coupling,
          double thickness,
          std::vector<CVector> demag_tensor,
          std::vector<CVector> dipole_tensor) : id(id), mag(mag), anis(anis), K(K), Ms(Ms), coupling(coupling),
                                                thickness(thickness), demag_tensor(demag_tensor),
                                                dipole_tensor(dipole_tensor)
    {
    }

    double updateAnisotropy(double time) {}
    double updateCoupling(double time) {}
    CVector updateExternalField(double time) {}

    CVector Heff(double time)
    {
        double kVar, couplingVar = 0.0;
        CVector HextVar, Heff = {0., 0., 0.};
        kVar = updateAnisotropy(time);
        couplingVar = updateCoupling(time);
        HextVar = updateExternalField(time);
        K_log = K + kVar;
        coupling_log = coupling + couplingVar;

        Heff = calculateAnisotropy() +
               calculate_tensor_interaction(mag, demag_tensor, Ms) +
               calculate_tensor_interaction(mag, dipole_tensor, Ms);
    }
    CVector calculateAnisotropy()
    {
        double nom = (2 * K_log) * c_dot(anis, mag) / (MAGNETIC_PERMEABILITY * Ms);
        return anis * nom;
    }

    CVector llg(double time, CVector m)
    {
        CVector heff, prod, dmdt;
        heff = Heff(time);
        prod = c_cross(m, heff);
        dmdt = prod * GYRO - c_cross(m, prod) * GYRO * DAMPING;
        return dmdt;
    }

    void rk4_step(double time, double time_step)
    {
        CVector k1, k2, k3, k4, m_t;
        m_t = mag;
        k1 = llg(time, m_t) * time_step;
        k2 = llg(time + 0.5 * time_step, m_t + k1 * 0.5) * time_step;
        k3 = llg(time + 0.5 * time_step, m_t + k2 * 0.5) * time_step;
        k4 = llg(time + time_step, m_t + k3) * time_step;
        m_t = m_t + (k1 + k2 * 2.0 + (k3 * 2.0) + k4) / 6.0;
        m_t.normalize();
        mag = m_t;
    }
};