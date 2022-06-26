#include "../../core/compute.hpp"
#include "../../core/junction.hpp"
#include <algorithm>
#include <iostream>
#include <stdio.h>

typedef Layer<float> FLayer;
typedef CVector<float> FVector;

std::vector<float> generateRange(float start, float stop, float step)
{
    std::vector<float> ranges;
    float current = start;
    while (current < stop)
    {
        current += step;
        ranges.push_back(current);
    }
    return ranges;
}

int main(void)
{
    std::vector<FVector> demagTensor = {
        {0.00022708623583019705, 0., 0.},
        {0., 0.0011629799534817735, 0.},
        {0., 0., 0.99861}};

    std::vector<FVector> dipoleTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0.0, 0.0}};

    float damping = 0.01;

    const float temperature = 0.0;

    float surface = 1e-7;
    float Ms = 1.07;
    float thickness = 1e-9;

    FLayer l1("free",             // id
              FVector(0., 0., 1), // mag
              FVector(0, -0.0871557, 0.996195),
              Ms,           // Ms
              thickness,    // thickness
              surface,      // surface
              demagTensor,  // demag
              damping       // damping
    );

    FLayer l2("bottom",            // id
              FVector(0., 0., 1.), // mag
              FVector(0.34071865, -0.08715574, 0.936116),
              Ms,           // Ms
              thickness,    // thickness
              surface,      // surface
              demagTensor,  // demag
              damping       // damping

    );

    Junction<float> mtj(
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
    mtj.setLayerAnisotropyDriver("free", ScalarDriver<float>::getConstantDriver(305e3));
    mtj.setLayerAnisotropyDriver("bottom", ScalarDriver<float>::getConstantDriver(728e3));
    mtj.setIECDriver("free", "bottom", ScalarDriver<float>::getConstantDriver(4e-5));
    const float hmin = -800e3;
    const float hmax = 800e3;
    const int hsteps = 100;

    const float tStart = 2e-9;
    const float time = 4e-9;
    const float tStep = 4e-12;
    std::ofstream saveFile;
    saveFile.open("VSD_res.csv");
    saveFile << "H;f;Vmix;\n";

    const float HoePulseAmplitude = 5e2;

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
            const AxialDriver<float> HDriver(
                ScalarDriver<float>::getConstantDriver(H * sqrt(2) / 2),
                ScalarDriver<float>::getConstantDriver(H * sqrt(2) / 2),
                NullDriver<float>());

            const AxialDriver<float> HoeDriver(
                NullDriver<float>(),
                ScalarDriver<float>::getSineDriver(0, HoePulseAmplitude, f, 0),
                NullDriver<float>());

            mtj.setLayerExternalFieldDriver(
                "all",
                HDriver);
            mtj.setLayerOerstedFieldDriver(
                "all",
                HoeDriver);
            mtj.runSimulation(
                time,
                tStep, tStep, false, false, false);

            auto vsdMix = ComputeFunctions<float>::calculateVoltageSpinDiode(
                mtj.getLog(), "Rx", f, 1, tStart);

            saveFile << H << ";" << f << ";" << vsdMix["Vmix"] << ";" << std::endl;
        }
    }

    saveFile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total result retrieval time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
}
