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
    double currentDensity = 1e10;
    double beta = 1;
    double spinPolarisation = -1.0;
    double STT_ = 1;
    double temperature = 300;

    double sttOn = false;

    ScalarDriver anisotropyBottom(
        constant,
        1500e3);

    ScalarDriver anisotropyTop(
        constant,
        800e3 // constant
    );

    Layer l1("free",             // id
             CVector(0., 0., 1), // mag
             CVector(0, 0., 1.), // anis
             1200e3,             // Ms
             1.4e-9,             // thickness
             7e-10 * 7e-10,      // surface
             demagTensor,        // demag
             dipoleTensor,       // dipole
             temperature,        // temp
             sttOn,              // STT
             damping,            // damping
             STT_,               // Slonczewski spacer
             spinPolarisation,   // spin polaristaion
             beta                // beta
    );


    auto currentDriver = ScalarDriver::getConstantDriver(1.0);
    l1.setCurrentDriver(currentDriver);
    l1.setAnisotropyDriver(anisotropyTop);

    Layer l2("bottom",            // id
             CVector(0., 1., 1.), // mag
             CVector(0, 1., 1.),  // anis
             1000e3,              // Ms
             1.4e-9,              // thickness
             7e-10 * 7e-10,       // surface
             demagTensor,         // demag
             dipoleTensor,        // dipole
             temperature,         // temp
             sttOn,               // STT
             damping,             // damping
             STT_,                // Slonczewski spacer
             spinPolarisation,    // spin polaristaion
             beta                 // beta
    );

    l2.setAnisotropyDriver(anisotropyBottom);
    l2.setCurrentDriver(currentDriver);
    Junction mtj(
        {l1, l2}, "STT.csv");

    mtj.runSimulation(50e-9, 1e-13, 1e-11, true, true, true);
    std::cout << l1.calculateLayerCriticalSwitchingCurrent("IP") << std::endl;
    std::cout << l1.calculateLayerCriticalSwitchingCurrent("PP") << std::endl;
}