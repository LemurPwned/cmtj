#include "../../core/compute.hpp"
#include "../../core/stack.hpp"
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
    double surface = 70e-9 * 70e-9;
    double Ms = 1.07;
    double thickness = 1e-9;
    const double fL = 215e-10;
    const double dL = 300e-10;
    DLayer l1 = DLayer::LayerSOT("free",                    // id
                                 DVector(0., 0., 1.),       // mag
                                 DVector(0., -0.087, 0.99), // anis
                                 Ms,                        // Ms
                                 thickness,                 // thickness
                                 surface,                   // surface
                                 demagTensor,               // demag
                                 dipoleTensor,              // dipole
                                 damping,                   // damping
                                 fL,
                                 dL);

    DLayer l1m = DLayer::LayerSOT("free",              // id
                                  DVector(0., 0., 1.), // mag
                                  DVector(0., 0, 1.),  // anis
                                  Ms * 1.05,           // Ms
                                  thickness * 0.98,    // thickness
                                  surface,             // surface
                                  demagTensor,         // demag
                                  dipoleTensor,        // dipole
                                  damping,             // damping
                                  fL,
                                  dL);

    DVector p(1, 0, 0);
    l1.setReferenceLayer(p);
    l1m.setReferenceLayer(p);

    const double f = 7e9;
    Junction<double> mtj1(
        {l1}, "", 100, 150);
    mtj1.setLayerAnisotropyDriver("free", ScalarDriver<double>::getConstantDriver(305e3));
    mtj1.setLayerCurrentDriver("free", ScalarDriver<double>::getSineDriver(1e8, 1e3, f, 0));

    Junction<double> mtj2(
        {l1m}, "", 100, 150);

    mtj2.setLayerAnisotropyDriver("free", ScalarDriver<double>::getConstantDriver(345e3));
    mtj2.setLayerCurrentDriver("free", ScalarDriver<double>::getSineDriver(1e8, 1e3, f, 0));

    // set external field
    const double H = 95492.965855;
    const double theta = 70 * M_PI / 180;
    const double phi = 0 * M_PI / 180;
    const AxialDriver<double> HDriver(
        ScalarDriver<double>::getConstantDriver(H * sin(theta) * cos(phi)),
        ScalarDriver<double>::getConstantDriver(H * sin(theta) * sin(phi)),
        ScalarDriver<double>::getConstantDriver(H * cos(theta)));

    mtj1.setLayerExternalFieldDriver(
        "all",
        HDriver);
    mtj2.setLayerExternalFieldDriver(
        "all",
        HDriver);

    std::vector<Junction<double>> s = {mtj1, mtj2};
    SeriesStack<double> stack(s);
    stack.setCouplingStrength(0.1);
    const double jmin = 0;
    const double jmax = 3e10;
    const int jsteps = 100;

    const double tStart = 5e-9;
    const double time = 20e-9;
    const double tStep = 4e-12;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto Jdist = generateRange(jmin, jmax, (jmax - jmin) / jsteps);
    std::cout << "Generated range" << std::endl;
    std::ofstream saveFile;

    saveFile.open("Jspectrum.csv");
    saveFile << "J;frequency;amplitude\n";
    for (int i = 0; i < Jdist.size(); i++)
    {
        // we need to reset the magnetisation
        stack.setMagnetisation(0, "free", DVector(0, 0, 1));
        stack.setMagnetisation(1, "free", DVector(0, 0, 1));
        stack.setCouplingStrength(0.1);
        stack.clearLogs();
        // std::cout << stack.junctionList[0].logLength << std::endl;
        stack.setCoupledCurrentDriver(ScalarDriver<double>::getConstantDriver(Jdist[i]));
        stack.runSimulation(time, tStep, tStep);
        if (i == 0)
        {
            stack.saveLogs("test.csv");
        }
        auto log = stack.getLog();
        // for (int j = 0; j < log["time"].size(); j++)
        // {
        //     saveFile << Jdist[i] << ";" << log["time"][j] << ";" << log["Resistance"][j] << std::endl;
        // }
        auto Jspectrum = ComputeFunctions<double>::spectralFFT(
            stack.getLog(), {"Resistance"}, tStart, tStep);
        for (int j = 1; j < floor(Jspectrum["frequencies"].size() / 4); j++)
        {
            saveFile << Jdist[i] / 1e8 << ";" << Jspectrum["frequencies"][j] / 1e9 << ";" << Jspectrum["Resistance_amplitude"][j] << std::endl;
        }
    }

    saveFile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total result retrieval time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
}