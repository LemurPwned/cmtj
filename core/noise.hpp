/**
 * @file noise.hpp
 * @author Jakub
 * @brief One F generator, based on the Pink Noise generator from the Music DSP
 * https://www.musicdsp.org/en/latest/Synthesis/220-trammell-pink-noise-c-class.html
 * Second version is custom and gives better results, but builds on the initial one.
 * @version 1.0
 * @date 2022-03-22
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _PinkNoise_H
#define _PinkNoise_H

#include <algorithm>               // for generate, sort, unique
#include <iterator>                // for distance
#include <cstdlib>                 // for rand, srand, NULL, RAND_MAX, size_t
#include <ctime>                   // for time
#include <numeric>                 // for accumulate
#include <random>                  // for uniform_real_distribution, geometr...
#include <vector>                  // for vector

#define PINK_NOISE_NUM_STAGES 3
template <typename T>
class PinkNoise {
public:
    PinkNoise() {
        srand(time(NULL)); // initialize random generator
        clear();
    }

    void clear() {
        for (size_t i = 0; i < PINK_NOISE_NUM_STAGES; i++)
            state[i] = 0.0;
    }

    T tick() {
        static const T RMI2 = 2.0 / T(RAND_MAX); // + 1.0; // change for range [0,1)
        static const T offset = A[0] + A[1] + A[2];

        // unrolled loop
        T temp = T(rand());
        state[0] = P[0] * (state[0] - temp) + temp;
        temp = T(rand());
        state[1] = P[1] * (state[1] - temp) + temp;
        temp = T(rand());
        state[2] = P[2] * (state[2] - temp) + temp;
        return (A[0] * state[0] + A[1] * state[1] + A[2] * state[2]) * RMI2 - offset;
    }

protected:
    T state[PINK_NOISE_NUM_STAGES];
    static constexpr T A[PINK_NOISE_NUM_STAGES] = { 0.02109238, 0.07113478, 0.68873558 };
    static constexpr T P[PINK_NOISE_NUM_STAGES] = { 0.3190,  0.7756,  0.9613 };
};



template<typename T = double>
class OneFNoise {
private:
    int sources;
    std::vector<T> state;
    std::geometric_distribution<int> geom_distr;
    // Mersenne twister is higher quality than the default one
    std::mt19937 generator;
    std::uniform_real_distribution<T> float_dist;
    std::vector<int> trials;

    T scale = 1;
    T sumTrack = 0;
public:
    OneFNoise(int sources, T bias, T scale) : sources(sources), geom_distr(bias), scale(scale) {
        this->state.resize(sources); // fill it with 0s
        this->trials.resize(sources);
        this->float_dist = std::uniform_real_distribution<T>(0, 1);
        // start off with random values in the state
        std::generate(this->state.begin(), this->state.end(), [&] { return this->float_dist(generator);});
        // try out the binding stuff
    }
    /**
     * @brief This function works faster if p is a large number (p > 0.5)
     *
     * @return T sum of the state
     */
    T tick() {
        std::generate(this->trials.begin(), this->trials.end(), [&] { return this->geom_distr(generator);});
        std::sort(this->trials.begin(), this->trials.end());
        const auto uniq = std::unique(this->trials.begin(), this->trials.end());
        // compute the distance of the last unique element
        // this basically takes only the unique elements of the trials
        // because if we repeatedly change the same index, we don't get any advantage
        const auto lastIndx = std::distance(this->trials.begin(), uniq);
        for (int i = 0; i < lastIndx; ++i) {
            const auto t = this->trials[i];
            if (t < this->sources) {
                this->state[t] = this->float_dist(generator);
            }
        }
        return this->scale * std::accumulate(state.begin(), state.end(), 0.);
    }

    /**
     * @brief This function works faster if the p is a small number (p < 0.5)
     *
     * @return T sum of the state
     */
    T tick2() {
        std::generate(this->trials.begin(), this->trials.end(), [&] { return this->geom_distr(generator);});
        for (const auto& t : this->trials) {
            if (t < this->sources) {
                this->state[t] = this->float_dist(generator);
            }
        }
        return this->scale * std::accumulate(state.begin(), state.end(), 0.);
    }
};


#endif
