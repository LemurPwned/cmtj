
#ifndef COMPUTE_FUNCTIONS_H
#define COMPUTE_FUNCTIONS_H
#include <algorithm>
#include <fftw3.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdio.h>
#include <string>

/**
 * Provides a static interface for computing useful magnetic properties
 * such as Voltage Spin Diode Effect or FFT on the magnetoresistance.
 */
class ComputeFunctions
{
public:
    /**
     * Computes the Voltage Spin diode.
     * @param log: This is the log from the simulation.
     * @param resTag: Tag to fetch the resistance from the simulation.
     * @param frequency: excitation frequency for the current.
     * @param power: power assumed in the system (this is somewhat an arbitrary value).
     * @param minTime: time after which to take the log, preferably when the magnetisation is stable.
     */
    static std::map<std::string, double> calculateVoltageSpinDiode(
        std::map<std::string, std::vector<double>> &log,
        const std::string resTag,
        double frequency, double power = 10e-6, const double minTime = 10e-9)
    {
        if (log.empty())
        {
            throw std::invalid_argument("Empty log! Cannot proceed without running a simulation!");
        }
        const double omega = 2 * M_PI * frequency;
        std::vector<double> &resistance = log[resTag];
        auto it = std::find_if(log["time"].begin(), log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });
        // turn into index
        const int thresIdx = (int)(log["time"].end() - it);
        const int cutSize = log["time"].size() - thresIdx;
        // Rpp
        const double RppMax = *std::max_element(resistance.begin() + thresIdx, resistance.end());
        const double RppMin = *std::min_element(resistance.begin() + thresIdx, resistance.end());
        const double avgR = std::accumulate(resistance.begin() + thresIdx, resistance.end(), 0.0) / cutSize;
        const double Iampl = sqrt(power / avgR);
        std::vector<double> voltage, current;
        std::transform(
            log["time"].begin() + thresIdx, log["time"].end(),
            std::back_inserter(current),
            [&Iampl, &omega](const double &time) { return Iampl * sin(omega * time); });

        for (unsigned int i = 0; i < cutSize; i++)
        {
            voltage.push_back(resistance[thresIdx + i] * current[i]);
        }
        const double Vmix = std::accumulate(voltage.begin(), voltage.end(), 0.0) / voltage.size();
        std::map<std::string, double> mRes = {{"Vmix", Vmix}, {"RMax", RppMax}, {"RMin", RppMin}, {"Rpp", (RppMax - RppMin)}};
        return mRes;
    }

    /**
     * Computes the FFT on a given tag.
     * @param log: This is the log from the simulation.
     * @param fftIds: a vector of ids (log keys) for which FFT is to be computed.
     * @param minTime: minimum waiting time (10e-9) by default. Set it so that non-harmonic
     * oscillations are not included into FFT computation.
     * @param timeStep: integration step (1e-11) by default .
     */
    static std::map<std::string, std::vector<double>>
    spectralFFT(std::map<std::string, std::vector<double>> &log,
                const std::vector<std::string> &fftIds,
                double minTime = 10.0e-9, double timeStep = 1e-11)
    {

        if (minTime >= log["time"][log["time"].size() - 1])
        {
            throw std::invalid_argument("The minTime parameter is larger than the simulation time!");
        }
        if (log.empty())
        {
            throw std::invalid_argument("Empty log! Cannot proceed without running a simulation!");
        }
        const std::vector<std::string> vectorNames = {"x", "y", "z"};

        // std::cout << log["free_mx"][0] << " " << log["time"].size() << std::endl;

        auto it = std::find_if(log["time"].begin(), log["time"].end(),
                               [&minTime](const auto &value) { return value >= minTime; });

        const int thresIdx = (int)(log["time"].end() - it);
        const int cutSize = log["time"].size() - thresIdx;

        // plan creation is not thread safe
        const double normalizer = timeStep * cutSize;
        const int maxIt = (cutSize % 2) ? cutSize / 2 : (cutSize - 1) / 2;
        std::vector<double> frequencySteps(maxIt);
        frequencySteps[0] = 0;
        for (int i = 1; i <= maxIt; i++)
        {
            frequencySteps[i - 1] = (i - 1) / normalizer;
        }
        // plan creation is not thread safe
        std::map<std::string, std::vector<double>> spectralFFTResult;
        spectralFFTResult["frequencies"] = std::move(frequencySteps);

        for (const auto &tag : fftIds)
        {
            if (log.find(tag) == log.end())
            {
                // not found
                throw std::invalid_argument("FFT id tag was not found in the junction log: " + tag);
            }
            std::vector<double> cutMag(log[tag].begin() + thresIdx, log[tag].end());
            // define FFT plan
            fftw_complex out[cutMag.size()];
            fftw_plan plan = fftw_plan_dft_r2c_1d(cutMag.size(),
                                                  cutMag.data(),
                                                  out,
                                                  FFTW_ESTIMATE); // here it's weird, FFT_FORWARD produces an empty plan

            if (plan == NULL)
            {
                throw std::runtime_error("Plan creation for fftw failed, cannot proceed");
            }
            fftw_execute(plan);
            const int outBins = (cutMag.size() + 1) / 2;
            std::vector<double> amplitudes;
            amplitudes.push_back(out[0][0]);
            for (int i = 1; i < outBins; i++)
            {
                const auto tandem = out[i];
                const double real = tandem[0];
                const double img = tandem[1];
                amplitudes.push_back(sqrt(pow(real, 2) + pow(img, 2)));
            }
            spectralFFTResult[tag + "_amplitude"] = std::move(amplitudes);
            fftw_destroy_plan(plan);
        }
        return spectralFFTResult;
    }
};

#endif