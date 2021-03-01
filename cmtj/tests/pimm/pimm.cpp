#include "../../core/compute.hpp"
#include "../../core/junction.hpp"
#include <algorithm>
#include <iostream>
#include <stdio.h>

enum H_MODE
{
    MAG = 1,
    THETA,
    PHI
};
typedef std::tuple<std::vector<CVector>, std::vector<double>, double> Hspace;
Hspace calculateHdistribution(double Hmag, double theta, double phi,
                              double minVal, double maxVal, int steps, H_MODE mode)
{
    std::vector<CVector> H;
    double step = (maxVal - minVal) / steps;
    std::vector<double> valueSpace; // # (steps);
    for (int i = 0; i < steps; i++)
    {
        valueSpace.push_back(minVal + step * i);
    }
    for (const auto &v : valueSpace)
    {
        if (mode == MAG)
        {
            Hmag = v;
        }
        else if (mode == THETA)
        {
            theta = v;
        }
        else if (mode == PHI)
        {
            phi = v;
        }
        H.push_back(CVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)) * Hmag);
    }
    return std::make_tuple(H, valueSpace, step);
}

int main(void)
{
    std::vector<CVector> demagTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0., 1.0}};

    std::vector<CVector> dipoleTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0.0, 0.0}};

    double damping = 0.035;

    double sttOn = false;
    const double temperature = 0.0;

    double surface = 1e-7;
    double Ms = 1.07 * TtoAm;
    double thickness = 1e-9;

    Layer l1("free",             // id
             CVector(0., 0., 1), // mag
             Ms,                 // Ms
             thickness,          // thickness
             surface,            // surface
             demagTensor,        // demag
             dipoleTensor,       // dipole
             temperature,        // temp
             false,              // STT
             damping             // damping
    );

    Layer l2("bottom",            // id
             CVector(0., 0., 1.), // mag
             Ms,                  // Ms
             thickness,           // thickness
             surface,             // surface
             demagTensor,         // demag
             dipoleTensor,        // dipole
             temperature,         // temp
             false,               // STT
             damping              // damping

    );

    auto iecD = ScalarDriver::getConstantDriver(1e-9);
    l1.setIECDriver(iecD);
    iecD = ScalarDriver::getConstantDriver(4e-5);
    l2.setIECDriver(iecD);
    Junction mtj(
        {l1, l2}, "", 100, 105);
    mtj.setLayerAnisotropyDriver("free", AxialDriver(CVector(0, 0, 305e3)));
    mtj.setLayerAnisotropyDriver("bottom", AxialDriver(CVector(0, 0, 728e3)));

    const double hmin = -800e3;
    const double hmax = 800e3;
    const int hsteps = 40;

    const double tStart = 5e-9;
    const double time = 10e-9;
    const double tStep = 1e-12;

    const double theta = 90;
    const double phi = 45;

    std::ofstream saveFile;
    saveFile.open("PIMM_res.csv");
    saveFile << "H;f_x;f_y;f_z;\n";
    auto HspaceVals = calculateHdistribution(0, theta, phi, hmin, hmax, hsteps, MAG);

    auto Hdist = std::get<0>(HspaceVals);
    auto itValues = std::get<1>(HspaceVals);
    auto step = std::get<2>(HspaceVals);
    int indx = 0;
    CVector HoeDir(0, 1, 0);
    const double HoePulseAmplitude = 397.88;
    const double pulseStart = 0.0e-9;
    const double pulseStop = 0.1e-9;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const std::vector<std::string> tagIds = {"free_mx", "free_my", "free_mz"};
    for (auto &H : Hdist)
    {
        // std::cout << "H: (" << H.x << "," << H.y << "," << H.z << ") " << H.length() << std::endl;
        mtj.clearLog();
        const AxialDriver HDriver(
            ScalarDriver::getConstantDriver(H.x),
            ScalarDriver::getConstantDriver(H.y),
            ScalarDriver::getConstantDriver(H.z));

        AxialDriver HoeDriver(
            ScalarDriver::getStepDriver(0, HoePulseAmplitude, pulseStart, pulseStop),
            ScalarDriver::getStepDriver(0, HoePulseAmplitude, pulseStart, pulseStop),
            ScalarDriver::getStepDriver(0, HoePulseAmplitude, pulseStart, pulseStop));

        HoeDriver.applyMask(HoeDir);

        mtj.setLayerExternalFieldDriver(
            "all",
            HDriver);
        mtj.setLayerOerstedFieldDriver(
            "all",
            HoeDriver);
        mtj.runSimulation(
            time,
            tStep, tStep, false, false, false);
        auto fftResult = ComputeFunctions::spectralFFT(
            mtj.getLog(), tagIds, tStart, tStep);

        saveFile << itValues[indx] << ";";
        for (const std::string tag : {"x", "y", "z"})
        {
            const double maxAmplitude = *std::max_element(
                fftResult["free_m" + tag + "_amplitude"].begin() + 1,
                fftResult["free_m" + tag + "_amplitude"].end());
            saveFile << maxAmplitude << ";";
        }
        saveFile << std::endl;
        indx += 1;
    }

    saveFile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total result retrieval time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
}