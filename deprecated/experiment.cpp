#include "junction.hpp"

void threadedSimulation(Junction cjx, double minField, double maxField,
                        int numberOfPoints, std::ofstream &vsdFile) {
  const int threadNum = std::thread::hardware_concurrency() - 2;
  std::vector<std::future<std::vector<std::tuple<double, double>>>>
      threadResults;
  threadResults.reserve(threadNum);

  const int pointsPerThread = numberOfPoints / threadNum + 1;
  const double spacing = (maxField - minField) / numberOfPoints;
  for (int i = 0; i < threadNum; i++) {
    const double threadMinField = pointsPerThread * i * spacing;
    const double threadMaxField = pointsPerThread * (i + 1) * spacing;
    threadResults.emplace_back(std::async([cjx, threadMinField, threadMaxField,
                                           spacing]() mutable {
      std::vector<std::tuple<double, double>> resAcc;
      const double freq = 7e9;
      for (double field = threadMinField; field < threadMaxField;
           field += spacing) {
        cjx.setConstantExternalField((field / 1000) * TtoAm, xaxis);
        // cjx.setLayerAnisotropyUpdate("free", 12000, freq, 0);
        // cjx.setLayerAnisotropyUpdate("bottom", 12000, freq, 0);
        cjx.setLayerCoupling("free", -3e-6);
        cjx.setLayerCoupling("bottom", -3e-6);
        cjx.setLayerIECUpdate("free", 1e-6, freq, 0);
        cjx.setLayerIECUpdate("bottom", 1e-6, freq, 0);
        cjx.runSimulation(20e-9);
        std::map<std::string, double> vsd = cjx.calculateVoltageSpinDiode(freq);
        resAcc.push_back({field, vsd["Vmix"]});
        cjx.log.clear();
      }

      return resAcc;
    }));
  }

  vsdFile << "H;Vmix\n";
  for (auto &result : threadResults) {
    for (const auto [field, vsdVal] : result.get()) {
      vsdFile << field << ";" << vsdVal << "\n";
    }
  };
}

void parameterScanVSD(Junction junction, std::string filename) {
  double minField = 000.0;
  double maxField = 500.0;
  int numPoints = 50;
  double spacing = (maxField - minField) / numPoints;
  std::ofstream vsdFile;
  vsdFile.open(filename);
  vsdFile << "F;J;eJ;H;Vmix;Rpp\n";
  std::vector<double> eJs = {1e-6, 9e-7, 8e-7, 7e-7, 6e-7,
                             5e-7, 4e-7, 3e-7, 2e-7, 1e-7};
  std::vector<double> Js = {-2e-6, -3e-6, -4e-6, -5e-6, -6e-6};
  std::vector<double> Fs = {5e9, 5.5e9, 6e9, 6.5e9, 7e9, 7.5e9, 8e9};
  // std::vector<double> Fs = {7e9};

  for (const double &F : Fs) {
    for (const double &J : Js) {
      for (const double &eJ : eJs) {
        for (double field = minField; field < maxField; field += spacing) {
          junction.setConstantExternalField((field / 1000) * TtoAm, xaxis);
          junction.setLayerCoupling("free", J);
          junction.setLayerCoupling("bottom", J);
          junction.setLayerIECUpdate("free", eJ, F, 0);
          junction.setLayerIECUpdate("bottom", eJ, F, 0);
          junction.runSimulation(20e-9, 1e-13, false, false);
          std::map<std::string, double> vsd =
              junction.calculateVoltageSpinDiode(F);
          vsdFile << F << ";" << J << ";" << eJ << ";" << field << ";"
                  << vsd["Vmix"] << ";" << vsd["Rpp"] << std::endl;
          junction.log.clear();
        }
      }
    }
  }

  vsdFile.close();
}

void parameterScanFFT(Junction &junction, std::string filename) {
  // std::vector<double> Js = {-4e-6, -3e-4, -2e-4, -1e-4, 0.0, 1e-4, 2e-4,
  // 3e-4, 4e-4}; std::vector<double> Js = {-9e-6, -8e-6, -7e-6, -6e-6, -5e-6,
  // -4e-6, -3e-6, -2e-6, -1e-6, 0.0,
  //                           1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6,
  //                           9e-6};
  std::vector<double> Js = {-9e-5, -8e-5, -7e-5, -6e-5, -5e-5, -4e-5, -3e-5,
                            -2e-5, -1e-5, 0.0,   1e-5,  2e-5,  3e-5,  4e-5,
                            5e-5,  6e-5,  7e-5,  8e-5,  9e-5};
  std::vector<double> Ks = {900e3, 950e3, 1000e3, 1050e3, 1100e3, 1150e3};
  std::vector<double> fixedFields = {150_mT, 200_mT, 250_mT, 300_mT,
                                     350_mT, 400_mT, 450_mT};

  std::ofstream fftFile;
  fftFile.open(filename);
  fftFile << "K;J;H;Mx;Ax;Px;My;Ay;Py;Mz;Az;Pz;\n";
  for (const double &J : Js) {
    for (const double &K : Ks) {
      for (const double &H : fixedFields) {
        junction.setConstantExternalField(H * TtoAm, xaxis);
        junction.setLayerCoupling("free", J);
        junction.setLayerCoupling("bottom", J);
        junction.setLayerAnisotropy("free", K);
        junction.setLayerStepUpdate("free", 1e-4 * TtoAm, 5.0_ns, 5.001_ns,
                                    xaxis);
        junction.setLayerStepUpdate("bottom", 1e-4 * TtoAm, 5.0_ns, 5.001_ns,
                                    xaxis);
        junction.runSimulation(20_ns, 1e-13, false, false);
        auto resMap = junction.calculateFFT(5e-9, 1e-11);
        junction.log.clear();
        fftFile << K << ";" << J << ";" << H << ";";
        for (const std::string majorKey : {"x", "y", "z"}) {
          for (const std::string minorKey :
               {"_resonant", "_amplitude", "_phase"}) {
            fftFile << resMap[majorKey + minorKey] << ";";
          }
        }
        fftFile << "\n";
      }
    }
  }
  fftFile.close();
}

int main(void) {

  std::vector<CVector> dipoleTensor = {{6.8353909454237E-4, 0., 0.},
                                       {0., 0.00150694452305927, 0.},
                                       {0., 0., 0.99780951638608}};
  std::vector<CVector> demagTensor = {{5.57049776248663E-4, 0., 0.},
                                      {0., 0.00125355500286346, 0.},
                                      {0., 0.0, -0.00181060482770131}};

  Layer l1("free",              // id
           CVector(0., 0., 1.), // mag
           CVector(0, 0., 1.),  // anis
           900e3,               // K
           1200e3,              // Ms
           -2.5e-6,             // J
           1.4e-9,              // thickness
           7e-10 * 7e-10,       // surface
           demagTensor,         // demag
           dipoleTensor);
  Layer l2("bottom",            // id
           CVector(0., 0., 1.), // mag
           CVector(0, 0., 1.),  // anis
           1500e3,              // K
           1000e3,              // Ms
           -2.5e-6,             // J
           7e-9,                // thickness
           7e-10 * 7e-10,       // surface
           demagTensor,         // demag
           dipoleTensor);

  Junction mtj({l1, l2}, "test2.csv");

  // double minField = 000.0;
  // double maxField = 600.0;
  // int numPoints = 50;
  // double spacing = (maxField - minField) / numPoints;
  // std::cout << spacing << std::endl;
  // std::ofstream vsdFile;
  // std::chrono::steady_clock::time_point begin =
  // std::chrono::steady_clock::now();

  // auto r = ComputeUtil::parallelFieldScan(mtj, minField, maxField, numPoints,
  //                                          [](Junction &mtj, const double
  //                                          field) mutable {
  //                                              mtj.log.clear();
  //                                              const double freq = 7e9;
  //                                              mtj.setConstantExternalField((field
  //                                              / 1000) * TtoAm, xaxis);
  //                                              mtj.setLayerCoupling("free",
  //                                              -3e-6);
  //                                              mtj.setLayerCoupling("bottom",
  //                                              -3e-6);
  //                                              mtj.setLayerIECUpdate("free",
  //                                              1e-6, freq, 0);
  //                                              mtj.setLayerIECUpdate("bottom",
  //                                              1e-6, freq, 0);
  //                                              mtj.runSimulation(20e-9);
  //                                              std::map<std::string, double>
  //                                              vsd =
  //                                              mtj.calculateVoltageSpinDiode(freq);
  //                                              mtj.log.clear();
  //                                              return std::make_tuple(field,
  //                                              vsd);
  //                                          });

  // ComputeUtil::customResultMap(r, "VSD-anisotropy.csv");
  // vsdFile.close();
  // mtj.setLayerAnisotropy("free", 200);
  // std::cout << mtj.layers[0].K << std::endl;
  // parameterScanFFT(mtj, "Jhigh_wide_H2_scan.csv");
  // parameterScanVSD(mtj, "VSD_scan_benchmark.csv");

  // vsdFile.open("VSD-IEC2.csv");
  // threadedSimulation(mtj, minField, maxField, numPoints, vsdFile);
  // vsdFile.close();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout
      << "Simulation time = "
      << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
      << "[s]" << std::endl;
  // parameterScanVSD(mtj, "VSD-anisotropy-lowIEC2.csv");

  // end = std::chrono::steady_clock::now();
  // std::cout << "Simulation time = " <<
  // std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() <<
  // "[s]" << std::endl;

  return 0;
}
