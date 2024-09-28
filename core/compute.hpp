
#ifndef COMPUTE_FUNCTIONS_H
#define COMPUTE_FUNCTIONS_H

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdio.h>
#include <string>
#include <unordered_map>

/**
 * Provides a static interface for computing useful magnetic properties
 * such as Voltage Spin Diode Effect or FFT on the magnetoresistance.
 */
template <typename T> class ComputeFunctions {
public:
  /**
   * Computes the Voltage Spin diode.
   * @param log: This is the log from the simulation.
   * @param resTag: Tag to fetch the resistance from the simulation.
   * @param frequency: excitation frequency for the current.
   * @param power: power assumed in the system (this is somewhat an arbitrary
   * value).
   * @param minTime: time after which to take the log, preferably when the
   * magnetisation is stable.
   */
  static std::unordered_map<std::string, T> calculateVoltageSpinDiode(
      std::unordered_map<std::string, std::vector<T>> &log,
      const std::string &resTag, T frequency, T power = 10e-6,
      const T minTime = 10e-9) {
    if (log.empty()) {
      throw std::invalid_argument(
          "Empty log! Cannot proceed without running a simulation!");
    }

    if (log.find(resTag) == log.end()) {
      // not found
      throw std::invalid_argument("Tag was not found in the junction log: " +
                                  resTag);
    }
    const T omega = 2 * M_PI * frequency;
    std::vector<T> &resistance = log[resTag];
    auto it = std::find_if(
        log["time"].begin(), log["time"].end(),
        [&minTime](const auto &value) { return value >= minTime; });
    // turn into index
    const int thresIdx = (int)(it - log["time"].begin());
    const int cutSize = log["time"].size() - thresIdx;

    // Rpp
    const T RppMax =
        *std::max_element(resistance.begin() + thresIdx, resistance.end());
    const T RppMin =
        *std::min_element(resistance.begin() + thresIdx, resistance.end());
    const T avgR =
        std::accumulate(resistance.begin() + thresIdx, resistance.end(), 0.0) /
        cutSize;
    const T Iampl = sqrt(power / avgR);
    std::vector<T> voltage, current;
    std::transform(
        log["time"].begin() + thresIdx, log["time"].end(),
        std::back_inserter(current),
        [&Iampl, &omega](const T &time) { return Iampl * sin(omega * time); });

    for (int i = 0; i < cutSize; i++) {
      voltage.push_back(resistance[thresIdx + i] * current[i]);
    }
    const T Vmix =
        std::accumulate(voltage.begin(), voltage.end(), 0.0) / voltage.size();
    std::unordered_map<std::string, T> mRes = {{"Vmix", Vmix},
                                               {"RMax", RppMax},
                                               {"RMin", RppMin},
                                               {"Rpp", (RppMax - RppMin)}};
    return mRes;
  }

  /**
   * Computes the FFT on a given tag.
   * @param log: This is the log from the simulation.
   * @param fftIds: a vector of ids (log keys) for which FFT is to be computed.
   * @param minTime: minimum waiting time (10e-9) by default. Set it so that
   * non-harmonic oscillations are not included into FFT computation.
   * @param timeStep: integration step (1e-11) by default .
   */
  static std::unordered_map<std::string, std::vector<T>>
  spectralFFT(std::unordered_map<std::string, std::vector<T>> &log,
              const std::vector<std::string> &fftIds, T minTime = 10.0e-9,
              T timeStep = 1e-11) {

    if (minTime >= log["time"][log["time"].size() - 1]) {
      throw std::invalid_argument(
          "The minTime parameter is larger than the simulation time!");
    }
    if (log.empty()) {
      throw std::invalid_argument(
          "Empty log! Cannot proceed without running a simulation!");
    }

    auto it = std::find_if(
        log["time"].begin(), log["time"].end(),
        [&minTime](const auto &value) { return value >= minTime; });
    const int thresIdx = (int)(it - log["time"].begin());
    const int cutSize = log["time"].size() - thresIdx;
    // plan creation is not thread safe
    const T normalizer = timeStep * cutSize;
    const int maxIt = (cutSize % 2) ? cutSize / 2 : (cutSize - 1) / 2;
    std::vector<T> frequencySteps(maxIt);
    for (int i = 1; i <= maxIt; i++) {
      frequencySteps[i - 1] = (i - 1) / normalizer;
    }
    // plan creation is not thread safe
    std::unordered_map<std::string, std::vector<T>> spectralFFTResult;
    spectralFFTResult["frequencies"] = std::move(frequencySteps);

    for (const auto &tag : fftIds) {
      if (log.find(tag) == log.end()) {
        // not found
        throw std::invalid_argument(
            "FFT id tag was not found in the junction log: " + tag);
      }
      std::vector<T> cutMag(log[tag].begin() + thresIdx, log[tag].end());
      // define FFT plan
      std::complex<T> *out = new std::complex<T>[cutMag.size()];
      fftw_plan plan = fftw_plan_dft_r2c_1d(
          cutMag.size(), cutMag.data(), reinterpret_cast<fftw_complex *>(out),
          FFTW_ESTIMATE); // here it's weird, FFT_FORWARD produces an empty plan

      if (plan == NULL) {
        throw std::runtime_error(
            "Plan creation for fftw failed, cannot proceed");
      }
      fftw_execute(plan);
      const int outBins = (cutMag.size() + 1) / 2;
      std::vector<T> amplitudes, phases;
      const double norm = (double)cutSize / 2;
      amplitudes.push_back(out[0].real());
      phases.push_back(0.);
      for (int i = 1; i < outBins - 1; i++) {
        const auto tandem = out[i];
        T real = tandem.real() / norm; //  [0];
        T img = tandem.imag() / norm;  // [1];
        amplitudes.push_back(sqrt(pow(real, 2) + pow(img, 2)));
        phases.push_back(atan2(img, real));
      }
      spectralFFTResult[tag + "_amplitude"] = std::move(amplitudes);
      spectralFFTResult[tag + "_phase"] = std::move(phases);
      fftw_destroy_plan(plan);
    }
    return spectralFFTResult;
  }

  static std::unordered_map<std::string, std::vector<T>>
  spectralFFTMixed(std::unordered_map<std::string, std::vector<T>> &log,
                   const std::vector<std::string> &tagsToMix,
                   T timeStep = 1e-11) {
    const int cutSize = log["time"].size();
    if (log.empty()) {
      throw std::invalid_argument(
          "Empty log! Cannot proceed without running a simulation!");
    }
    // plan creation is not thread safe
    const T normalizer = timeStep * cutSize;
    const int maxIt = (cutSize % 2) ? cutSize / 2 : (cutSize - 1) / 2;
    std::vector<T> frequencySteps(maxIt);
    frequencySteps[0] = 0;
    for (int i = 1; i <= maxIt; i++) {
      frequencySteps[i - 1] = (i - 1) / normalizer;
    }
    // plan creation is not thread safe
    std::unordered_map<std::string, std::vector<T>> spectralFFTResult;
    spectralFFTResult["frequencies"] = std::move(frequencySteps);

    std::vector<T> mixedSignal(log["time"].size(), 0);
    for (const auto &tag : tagsToMix) {
      if (log.find(tag) == log.end())
        // not found
        throw std::invalid_argument(
            "FFT id tag was not found in the junction log: " + tag);
      for (unsigned int i = 0; i < log["time"].size(); i++) {
        mixedSignal[i] += log[tag][i];
      }
    }

    // define FFT plan
    std::complex<T> *out = new std::complex<T>[mixedSignal.size()];
    fftw_plan plan = fftw_plan_dft_r2c_1d(
        mixedSignal.size(), mixedSignal.data(),
        reinterpret_cast<fftw_complex *>(out),
        FFTW_ESTIMATE); // here it's weird, FFT_FORWARD produces an empty plan

    if (plan == NULL) {
      throw std::runtime_error("Plan creation for fftw failed, cannot proceed");
    }
    fftw_execute(plan);
    const int outBins = (mixedSignal.size() + 1) / 2;
    std::vector<T> amplitudes;
    amplitudes.push_back(out[0].real());
    const double norm = (double)cutSize / 2;
    for (int i = 1; i < outBins; i++) {
      const auto tandem = out[i];
      T real = tandem.real() / norm; //  [0];
      T img = tandem.imag() / norm;  // [1];
      amplitudes.push_back(sqrt(pow(real, 2) + pow(img, 2)));
    }
    spectralFFTResult["mixed_amplitude"] = std::move(amplitudes);
    fftw_destroy_plan(plan);

    return spectralFFTResult;
  }
};

#endif
