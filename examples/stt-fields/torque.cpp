#include "../../core/compute.hpp"
#include "../../core/junction.hpp"
#include <algorithm>
#include <iostream>
#include <stdio.h>

typedef Layer<double> DLayer;

std::vector<double> generateRange(double start, double stop, double step, bool back)
{
    std::vector<double> ranges;
    double current = start;
    while (current < stop)
    {
        current += step;
        ranges.push_back(current);
    }
    if (back)
    {
        current = stop;
        while (current > start)
        {
            current -= step;
            ranges.push_back(current);
        }
    }
    return ranges;
}

int main(void)
{

    std::vector<DVector> demagTensor = {
        {0.00022708623583019705, 0., 0.},
        {0., 0.0011629799534817735, 0.},
        {0., 0., 0.999598310229469}};

    std::vector<DVector> dipoleTensor = {
        {0.0, 0., 0.},
        {0., 0.0, 0.},
        {0., 0.0, 0.0}};

    double damping = 0.004;
    double surface = 1;
    double Ms = 0.54;
    double thickness = 1.45e-9;

    const double Irf = 0.05 * 945.9653094793309;
    const double Hdl = 1200;
    const double Hfl = 430;
    std::cout << "Hdl: " << Hdl << " Hfl: " << Hfl << std::endl;

    DLayer l1("free",             // id
              DVector(0., 0., 1), // mag
              DVector(0.0, 0., 0.94),
              Ms,          // Ms
              thickness,   // thickness
              surface,     // surface
              demagTensor, // demag
              damping      // damping
    );

    DVector p(0, 1, 0);
    l1.setReferenceLayer(p);
    const double l = 3e-5;
    const double w = 2e-5;
    const double ratio = w / l;

    Junction<double> mtj(
        {l1},
        "",
        {186},            // Rx0
        {100},            // Rxy
        {-0.02},          // AMR_X
        {-0.02 * -ratio}, // AMR_Y
        {-0.25},          // SMR_X
        {-0.25 * ratio},  // SMR_y
        {-2.7}            // AHE
    );
    mtj.setLayerAnisotropyDriver("free", ScalarDriver<double>::getConstantDriver(1.8e5));

    const double hmin = -800e3;
    const double hmax = 800e3;
    const int hsteps = 50;

    const double theta = 89 * M_PI / 180;
    const double phi = 1 * M_PI / 180;

    const double tStart = 000e-9;
    const double time = 1200e-9;
    const double tStep = 8e-12;
    std::ofstream saveFile;
    saveFile.open("Torque_res.csv");
    saveFile << "H;Vmix;phase\n";
    // saveFile << "H;Vmix;indx\n";

    const double HoePulseAmplitude = 1.5e3 * Irf;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto frequencies = {0.8e9};
    auto Hdist = generateRange(hmin, hmax, (hmax - hmin) / hsteps, true);

    const std::string resTag = "Ry";

    std::cout << "Generated frequency range" << std::endl;
    // bottom, top mag
    // std::reverse(Hdist.begin(), Hdist.end());
    for (auto &f : frequencies)
    {
        std::cout << "Computing " << f << std::endl;
        for (auto &H : Hdist)
        {
            mtj.clearLog();
            const AxialDriver<double> HDriver(
                ScalarDriver<double>::getConstantDriver(H * sin(theta) * cos(phi)),
                ScalarDriver<double>::getConstantDriver(H * sin(theta) * sin(phi)),
                ScalarDriver<double>::getConstantDriver(H * cos(theta)));

            const AxialDriver<double> HoeDriver(
                NullDriver<double>(),
                ScalarDriver<double>::getSineDriver(0, -HoePulseAmplitude, f, 0),
                NullDriver<double>());

            // mtj.setLayerOerstedFieldDriver("all", HoeDriver);

            mtj.setLayerDampingLikeTorqueDriver(
                "free", ScalarDriver<double>::getSineDriver(
                            0, -Hdl, f, 0));
            mtj.setLayerFieldLikeTorqueDriver(
                "free", ScalarDriver<double>::getSineDriver(
                            0, Hfl, f, 0));

            mtj.setLayerExternalFieldDriver(
                "all",
                HDriver);

            mtj.runSimulation(
                time,
                tStep, tStep, false, false, false);

            auto log = mtj.getLog();
            // compute the mixing voltage
            std::vector<double> mixingVoltage;
            for (std::size_t i = 0; i < log[resTag].size(); ++i)
            {
                mixingVoltage.push_back(log[resTag][i] * Irf * sin(2 * M_PI * f * log["time"][i]));
                // saveFile << -H << ';' << log[resTag][i] * Irf * sin(2 * M_PI * f * log["time"][i]) << ";" << i << std::endl;
            }

            log["mixing_voltage"] = mixingVoltage;
            // calculate the FFT
            auto spectrum = ComputeFunctions<double>::spectralFFT(
                log, {"mixing_voltage"}, tStart, tStep);

            // find 1f and 2f spectra
            auto it1f = std::lower_bound(spectrum["frequencies"].begin(), spectrum["frequencies"].end(), f);
            auto it2f = std::lower_bound(spectrum["frequencies"].begin(), spectrum["frequencies"].end(), 2 * f);
            // auto it1f = std::find(spectrum["frequencies"].begin(), spectrum["frequencies"].end(), f);
            // auto it2f = std::find(spectrum["frequencies"].begin(), spectrum["frequencies"].end(), 2 * f);
            if (it1f == spectrum["frequencies"].end() || it2f == spectrum["frequencies"].end())
            {
                throw std::runtime_error("Increase T to fit in 2f and 1f frequencies!");
            }
            const int indx1f = it1f - spectrum["frequencies"].begin();
            const int indx2f = it2f - spectrum["frequencies"].begin();

            saveFile << -H << ";" << spectrum["mixing_voltage_amplitude"][indx1f] << ";" <<  -1*spectrum["mixing_voltage_amplitude"][indx2f] * cos(spectrum["mixing_voltage_phase"][indx2f]) << std::endl;
        }
    }

    saveFile.close();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total result retrieval time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
}