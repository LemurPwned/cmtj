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
typedef std::tuple<std::vector<CVector<double>>, std::vector<double>, double> Hspace;
Hspace calculateHdistribution(double Hmag, double theta, double phi,
                              double minVal, double maxVal, int steps, H_MODE mode)
{
    std::vector<CVector<double>> H;
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
        H.push_back(CVector<double>(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)) * Hmag);
    }
    return std::make_tuple(H, valueSpace, step);
}

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
    std::vector<CVector<double>> demagTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0., 1.0}};

    std::vector<CVector<double>> dipoleTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0.0, 0.0}};

    double damping = 0.01;

    double sttOn = false;
    const double temperature = 0.0;

    double surface = 0.0;
    double Ms = 1.07;
    double thickness = 1e-9;

    Layer<double> l1("free",                                   // id
                     CVector<double>(0., 0., 1),               // mag
                     CVector<double>(0, -0.0871557, 0.996195), // anis
                     Ms,                                       // Ms
                     thickness,                                // thickness
                     surface,                                  // surface
                     demagTensor,                              // demag
                     dipoleTensor,                             // dipole
                     temperature,                              // temp
                     false,                                    // STT
                     damping                                   // damping
    );

    Layer<double> l2("bottom",                                           // id
                     CVector<double>(0., 0., 1.),                        // mag
                     CVector<double>(0.34071865, -0.08715574, 0.936116), // anis
                     Ms,                                                 // Ms
                     thickness,                                          // thickness
                     surface,                                            // surface
                     demagTensor,                                        // demag
                     dipoleTensor,                                       // dipole
                     temperature,                                        // temp
                     false,                                              // STT
                     damping                                             // damping

    );

    Junction<double> mtj(
        {l1, l2}, "", 100, 105);
    mtj.setLayerAnisotropyDriver("free", ScalarDriver<double>::getConstantDriver(305e3));
    mtj.setLayerAnisotropyDriver("bottom", ScalarDriver<double>::getConstantDriver(728e3));
    mtj.setIECDriver("free", "bottom", ScalarDriver<double>::getConstantDriver(4e-5));

    const double hmin = -800e3;
    const double hmax = 800e3;
    const int hsteps = 100;

    const double time = 8e-9;
    const double tStep = 1e-13;

    const double theta = 90;
    const double phi = 45;

    std::ofstream saveFile;
    saveFile.open("PIMM_res.csv");
    auto Hdist = generateRange(hmin, hmax, (hmax - hmin) / hsteps);
    int indx = 0;
    CVector<double> HoeDir(0, 0, 1);
    const double HoePulseAmplitude = 10000;
    const double pulseStart = 0.0e-9;
    const double pulseStop = 1e-11;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const std::vector<std::string> tagIds = {"bottom_mz", "free_mz"};

    CVector<double> m_init_free(1, 1, 0);
    CVector<double> m_init_bottom(1, 1, 0);
    std::string freeLayer = "free";
    std::string bottomLayer = "bottom";
    std::string all = "all";
    const AxialDriver<double> HoeDriver(
        NullDriver<double>(),
        NullDriver<double>(),
        ScalarDriver<double>::getStepDriver(0, HoePulseAmplitude, pulseStart, pulseStop));
    mtj.setLayerOerstedFieldDriver(
        all,
        HoeDriver);
    for (auto &H : Hdist)
    {
        mtj.clearLog();
        mtj.setLayerMagnetisation(freeLayer, m_init_free);
        mtj.setLayerMagnetisation(bottomLayer, m_init_bottom);
        const AxialDriver<double> HDriver(
            ScalarDriver<double>::getConstantDriver(H * sqrt(2) / 2),
            ScalarDriver<double>::getConstantDriver(H * sqrt(2) / 2),
            NullDriver<double>());
        mtj.setLayerExternalFieldDriver(
            all,
            HDriver);

        mtj.runSimulation(
            time,
            1e-12, 1e-12, false, false, false);
        m_init_free = mtj.layers[0].mag;
        m_init_bottom = mtj.layers[1].mag;
        // write a sum of mzs to a file
        for (int i = 0; i < mtj.log["time"].size(); i++)
        {
            saveFile << ";" << (mtj.log["free_mz"][i] + mtj.log["bottom_mz"][i]) / 2;
        }
        if (indx == Hdist.size())
            break;
        saveFile << "\n";
        indx += 1;
    }

    saveFile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total result retrieval time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
}