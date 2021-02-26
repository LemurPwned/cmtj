#include "junction.hpp"
#include "stack.hpp"
#include <iostream>
#include <stdio.h>

int main(void)
{
    std::vector<CVector> demagTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0., 1.0}};

    std::vector<CVector> dipoleTensor = {
        {5.57049776248663E-4, 0., 0.},
        {0., 0.00125355500286346, 0.},
        {0., 0.0, -0.00181060482770131}};

    double damping = 0.03;
    double currentDensity = 1e10;
    double beta = 1;
    double spinPolarisation = 1.0;
    double STT_ = 1;
    double temperature = 1e4;

    double sttOn = true;

    double radius = 70e-9;
    double surface = 7e-10 * 7e-10; // M_PI * pow(radius, 2);
    double Ms = 1 * TtoAm;
    double topAnisotropy = 800e3;

    ScalarDriver anisotropyBottom(
        constant,
        1500e3);

    ScalarDriver anisotropyTop(
        constant,
        topAnisotropy // constant
    );

    Layer l1("free",              // id
             CVector(0., 0., -1), // mag
             CVector(0, 0., 1.),  // anis
             1200e3,                  // Ms
             1e-9,              // thickness
             surface,             // surface
             demagTensor,         // demag
             dipoleTensor,        // dipole
             temperature,         // temp
             sttOn,               // STT
             damping,             // damping
             STT_,                // Slonczewski spacer
             beta,                // beta
             spinPolarisation,    // spin polaristaion
             false                // silent
    );
    Layer l2("bottom",            // id
             CVector(0., 1., 1.), // mag
             CVector(0, 1., 1.),  // anis
             1000e3,                  // Ms
             3e-9,                // thickness
             surface,             // surface
             demagTensor,         // demag
             dipoleTensor,        // dipole
             temperature,         // temp
             false,               // STT
             damping,             // damping
             STT_,                // Slonczewski spacer
             beta,                // beta,
             spinPolarisation     // spin polaristaion
    );

    Junction mtj(
        {l1, l2}, "STT.csv");
    mtj.setLayerIECDriver("all", ScalarDriver::getConstantDriver(-2.5e-6));
    mtj.setLayerAnisotropyDriver("free", anisotropyTop);
    mtj.setLayerAnisotropyDriver("bottom", anisotropyBottom);
    mtj.setLayerCurrentDriver("free", ScalarDriver::getConstantDriver(currentDensity));
    mtj.runSimulation(150e-9, 1e-13, 1e-12, true, true, false);
}