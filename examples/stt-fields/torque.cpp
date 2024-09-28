#include "../../core/compute.hpp"
#include "../../core/junction.hpp"
#include <algorithm>
#include <iostream>
#include <stdio.h>

typedef Layer<double> DLayer;

std::vector<double> generateRange(double start, double stop, double step,
                                  bool back) {
  std::vector<double> ranges;
  double current = start;
  while (current < stop) {
    ranges.push_back(current);
    current += step;
  }
  if (back) {
    current = stop;
    while (current > start) {
      ranges.push_back(current);
      current -= step;
    }
  }
  return ranges;
}

int main(void) {

  std::vector<DVector> demagTensor = {
      {0.00024164288391924, 2.71396011566517e-10, 5.95503928124313e-14},
      {2.71396011566517e-10, 0.000160046006320031, 1.32504057070646e-14},
      {5.95503928124313e-14, 1.32504057070646e-14, 0.999598310229469}};
  // std::vector<DVector> demagTensor = {
  //     {0.0, 0., 0.},
  //     {0., 0.0, 0.},
  //     {0., 0., 0.98}};

  double damping = 0.004;
  double surface = 1;
  double Ms = 1.5; // 0.54 0.52
  double thickness = 1.45e-9;

  const double Irf = 5e-3; // 0.0065 / 2.1
  const double Hdl = -600; // 1200
  const double Hfl = -447; // 430
  std::cout << "Hdl: " << Hdl << " Hfl: " << Hfl << std::endl;

  DLayer l1("free",               // id
            DVector(.0, 0., 1.),  // mag
            DVector(0.0, .0, 1.), // 0.94 // 0.85
            Ms,                   // Ms
            thickness,            // thickness
            surface,              // surface
            demagTensor,          // demag
            damping               // damping
  );

  DVector p(0, 1, 0);
  l1.setReferenceLayer(p);
  const double l = 2e-5;
  const double w = 3e-5;
  const double ratio = w / l;

  // Junction<double> mtj(
  //     {l1},
  //     "",
  //     {186},            // Rx0
  //     {100},            // Rxy
  //     {-0.02},          // AMR_X
  //     {-0.02 * -ratio}, // AMR_Y
  //     {-0.25},          // SMR_X
  //     {-0.25 * ratio},  // SMR_y
  //     {-2.3}            // AHE
  // );

  Junction<double> mtj({l1}, "", {304.7}, // Rx0
                       {3},               // Rxy
                       {-0.466},          // AMR_X
                       {-0.466 * -ratio}, // AMR_Y
                       {-0.053},          // SMR_X
                       {-0.053 * ratio},  // SMR_y
                       {-5.7}             // AHE
  );

  double Ku = 1.e6; // 1.8e5 0.85
  mtj.setLayerAnisotropyDriver("free",
                               ScalarDriver<double>::getConstantDriver(Ku));

  const double hmin = -700e3;
  const double hmax = -hmin;
  const int hsteps = 80;

  const double theta = 89 * M_PI / 180;
  const double phi = 89 * M_PI / 180;

  const double tStart = 000e-9;
  const double time = 1200e-9;
  const double tStep = 1e-11;
  std::ofstream saveFile;
  saveFile.open("Torque_res.csv");
  saveFile << "H;Vmix;phase\n";
  // saveFile << "H;Vmix;indx\n";

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();
  const auto frequencies = {0.8e9};
  auto Hdist = generateRange(hmin, hmax, (hmax - hmin) / hsteps, false);

  const std::string resTag = "Ry";

  std::cout << "Generated frequency range" << std::endl;
  // bottom, top mag
  // std::reverse(Hdist.begin(), Hdist.end());
  for (auto &f : frequencies) {
    std::cout << "Computing " << f << std::endl;
    for (auto &H : Hdist) {
      mtj.clearLog();
      const AxialDriver<double> HDriver(
          ScalarDriver<double>::getConstantDriver(H * sin(theta) * cos(phi)),
          ScalarDriver<double>::getConstantDriver(H * sin(theta) * sin(phi)),
          ScalarDriver<double>::getConstantDriver(H * cos(theta)));
      // const AxialDriver<double> HoeDriver(
      // ScalarDriver<double>::getSineDriver(0, 1000, f, 0),
      // NullDriver<double>(), NullDriver<double>());
      // mtj.setLayerOerstedFieldDriver("all", HoeDriver);
      // mtj.setLayerCurrentDriver("all",
      //                           ScalarDriver<double>::getSineDriver(
      //                               0, jrf, f, 0));

      mtj.setLayerDampingLikeTorqueDriver(
          "free", ScalarDriver<double>::getSineDriver(0, Hdl, f, 0));
      mtj.setLayerFieldLikeTorqueDriver(
          "free", ScalarDriver<double>::getSineDriver(0, Hfl, f, 0));

      mtj.setLayerExternalFieldDriver("all", HDriver);

      mtj.runSimulation(time, tStep, tStep, false, false, false);

      auto log = mtj.getLog();
      // compute the mixing voltage
      std::vector<double> mixingVoltage;
      for (std::size_t i = 0; i < log[resTag].size(); ++i) {
        const double multiplierHann =
            0.5 * (1 - cos(2 * M_PI * i / (log[resTag].size() - 1)));
        const double v =
            log[resTag][i] * Irf * sin(2 * M_PI * f * log["time"][i]);
        mixingVoltage.push_back(v * multiplierHann);
        // mixingVoltage.push_back(v);
        // saveFile << -H << ';' << v << ";" << i << std::endl;
      }

      // const double maxV = *std::max_element(mixingVoltage.begin(),
      // mixingVoltage.end()); const double maxR =
      // *std::max_element(log[resTag].begin(), log[resTag].end()); std::cout <<
      // "Max V: " << maxV << "; Max R: " << maxR << std::endl;
      log["mixing_voltage"] = mixingVoltage;
      // calculate the FFT
      auto spectrum = ComputeFunctions<double>::spectralFFT(
          log, {"mixing_voltage"}, tStart, tStep);

      // find 1f and 2f spectra
      auto it1f = std::lower_bound(spectrum["frequencies"].begin(),
                                   spectrum["frequencies"].end(), f);
      auto it2f = std::lower_bound(spectrum["frequencies"].begin(),
                                   spectrum["frequencies"].end(), 2 * f);
      if (it1f == spectrum["frequencies"].end() ||
          it2f == spectrum["frequencies"].end()) {
        throw std::runtime_error("Increase T to fit in 2f and 1f frequencies!");
      }
      const int indx1f = it1f - spectrum["frequencies"].begin();
      const int indx2f = it2f - spectrum["frequencies"].begin();

      saveFile << H << ";" << spectrum["mixing_voltage_amplitude"][indx1f]
               << ";"
               << spectrum["mixing_voltage_amplitude"][indx2f] *
                      cos(spectrum["mixing_voltage_phase"][indx2f])
               << std::endl;
    }
  }

  saveFile.close();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout
      << "Total result retrieval time = "
      << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
      << std::endl;
}
