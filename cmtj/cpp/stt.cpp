#include "junction.hpp"
#include "stack.hpp"
#include <stdio.h>
#include <iostream>

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
    double currentDensity = 0.6e11;
    double beta = 1;
    double spinPolarisation = -1.0;
    double STT_ = 1;
    double temperature = 300;

    double sttOn = true;

    double radius = 60e-9;
    double surface = M_PI * pow(radius, 2);
    double Ms = 1 * TtoAm;
    double topThickness = 1e-9;
    double topAnisotropy = 600e3;

    ScalarDriver anisotropyBottom(
        constant,
        900e3);

    ScalarDriver anisotropyTop(
        constant,
        topAnisotropy // constant
    );

    double Hkp = (2 * topAnisotropy) - 4 * M_1_PI * Ms;
    std::cout << "Hkp: " << Hkp << std::endl;
    std::cout << "4Ms: " << 4 * M_1_PI * Ms << std::endl;

    Layer l1("free",             // id
             CVector(0., 0., 1), // mag
             CVector(0, 0., 1.), // anis
             Ms,                 // Ms
             topThickness,       // thickness
             surface,            // surface
             demagTensor,        // demag
             dipoleTensor,       // dipole
             temperature,        // temp
             sttOn,              // STT
             damping,            // damping
             STT_,               // Slonczewski spacer
             beta,               // beta
             spinPolarisation,   // spin polaristaion
             false               // silent
    );

    auto currentDriver = ScalarDriver::getConstantDriver(currentDensity);
    l1.setCurrentDriver(currentDriver);
    l1.setAnisotropyDriver(anisotropyTop);

    Layer l2("bottom",            // id
             CVector(0., 0., 1.), // mag
             CVector(0, 0., 1.),  // anis
             Ms,                  // Ms
             1.4e-9,              // thickness
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

    l2.setAnisotropyDriver(anisotropyBottom);
    Junction mtj(
        {l1, l2}, "STT.csv");

    mtj.runSimulation(2e-9, 1e-13, 1e-13, true, true, false);
    std::cout << l1.calculateLayerCriticalSwitchingCurrent("IP") << std::endl;
    std::cout << l1.calculateLayerCriticalSwitchingCurrent("PP") << std::endl;
}