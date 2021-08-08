#include "../../core/compute.hpp"
#include "../../core/junction.hpp"
#include <algorithm>
#include <iostream>
#include <stdio.h>

typedef Layer<double> DLayer;

std::vector<double> generateRange(double start, double stop, double step)
{
    std::vector<double> ranges;
    double current = start;
    while (current < stop)
    {
        current += step;
        ranges.push_back(current);
    }
    return ranges;
}

int main(void)
{

    std::vector<DVector> demagTensor = {
        {0.00022708623583019705, 0., 0.},
        {0., 0.0011629799534817735, 0.},
        {0., 0., 0.99861}};

    std::vector<DVector> dipoleTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0.0, 0.0}};

    double damping = 0.01;
    double surface = 70e-9*70e-9;
    double Ms = 1.07;
    double thickness = 1e-9;

    DLayer l1("free",             // id
              DVector(0., 0., 1), // mag
              DVector(0, -0.0871557, 0.996195),
              Ms,           // Ms
              thickness,    // thickness
              surface,      // surface
              demagTensor,  // demag
              dipoleTensor, // dipole
              damping       // damping
    );

    DLayer l2("bottom",            // id
              DVector(0., 0., 1.), // mag
              DVector(0.34071865, -0.08715574, 0.936116),
              Ms,           // Ms
              thickness,    // thickness
              surface,      // surface
              demagTensor,  // demag
              dipoleTensor, // dipole
              damping       // damping

    );

    Junction<double> mtj(
        {l1, l2},
        "",
        {100, 100},       // Rx0
        {0, 0},           // Rxy
        {10, 10},         // AMR_X
        {30 / 2, 30 / 2}, // AMR_Y
        {0, 0},           // SMR_X
        {0, 0},           // SMR_y
        {0, 0}            // AHE
    );
    mtj.setLayerAnisotropyDriver("free", ScalarDriver<double>::getConstantDriver(305e3));
    mtj.setLayerAnisotropyDriver("bottom", ScalarDriver<double>::getConstantDriver(728e3));
    mtj.setIECDriver("free", "bottom", ScalarDriver<double>::getConstantDriver(4e-5));
    // mtj.setLayerTemperatureDriver("all", ScalarDriver<double>::getConstantDriver(300));
    const double hmin = -800e3;
    const double hmax = 800e3;
    const int hsteps = 100;

    const double tStart = 1e-9;
    const double time = 4e-9;
    const double tStep = 4e-12;
    std::ofstream saveFile;
    saveFile.open("VSD_res.csv");
    saveFile << "H;f;Vmix;\n";

    const double HoePulseAmplitude = 5e2;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto frequencies = generateRange(0e9, 48e9, 1e9);
    const auto Hdist = generateRange(hmin, hmax, (hmax - hmin) / hsteps);
    std::cout << "Generated frequency range" << std::endl;
    for (auto &f : frequencies)
    {
        std::cout << "Computing " << f << std::endl;
        for (auto &H : Hdist)
        {
            mtj.clearLog();
            const AxialDriver<double> HDriver(
                ScalarDriver<double>::getConstantDriver(H * sqrt(2) / 2),
                ScalarDriver<double>::getConstantDriver(H * sqrt(2) / 2),
                NullDriver<double>());

            const AxialDriver<double> HoeDriver(
                NullDriver<double>(),
                ScalarDriver<double>::getSineDriver(0, HoePulseAmplitude, f, 0),
                NullDriver<double>());

            mtj.setLayerExternalFieldDriver(
                "all",
                HDriver);
            mtj.setLayerOerstedFieldDriver(
                "all",
                HoeDriver);
            mtj.runSimulation(
                time,
                tStep, tStep, false, false, false);

            auto vsdMix = ComputeFunctions<double>::calculateVoltageSpinDiode(
                mtj.getLog(), "Rx", f, 1, tStart);
            saveFile << H << ";" << f << ";" << vsdMix["Vmix"] << ";" << std::endl;
        }
    }

    saveFile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total result retrieval time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
}