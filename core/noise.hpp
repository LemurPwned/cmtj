/**
 * @file noise.hpp
 * @author Jakub
 * @brief One F generator, based on the Pink Noise generator from the Music DSP
 * https://www.musicdsp.org/en/latest/Synthesis/220-trammell-pink-noise-c-class.html
 * Second version is custom and gives better results, but builds on the initial
 * one.
 * @version 1.0
 * @date 2022-03-22
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _PinkNoise_H
#define _PinkNoise_H

#include "../third_party/kissfft/kissfft.hh"
#include "cvector.hpp"
#include <algorithm> // for generate, sort, unique
#include <complex>
#include <cstdlib>  // for rand, srand, NULL, RAND_MAX, size_t
#include <ctime>    // for time
#include <iterator> // for distance
#include <memory>
#include <numeric> // for accumulate
#include <random>  // for uniform_real_distribution, geometr...
#include <vector>  // for vector

#define PINK_NOISE_NUM_STAGES 3
template <typename T> class PinkNoise {
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
    return (A[0] * state[0] + A[1] * state[1] + A[2] * state[2]) * RMI2 -
           offset;
  }

protected:
  T state[PINK_NOISE_NUM_STAGES];
  static constexpr T A[PINK_NOISE_NUM_STAGES] = {0.02109238, 0.07113478,
                                                 0.68873558};
  static constexpr T P[PINK_NOISE_NUM_STAGES] = {0.3190, 0.7756, 0.9613};
};

template <typename T = double> class NullTicker {
public:
  explicit NullTicker() {}
  ~NullTicker() {}
  virtual T tick() { return 0; }
};

template <typename T = double> class OneFNoise {
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
  OneFNoise(int sources, T bias, T scale)
      : sources(sources), geom_distr(bias), scale(scale) {
    this->state.resize(sources); // fill it with 0s
    this->trials.resize(sources);
    this->float_dist = std::uniform_real_distribution<T>(0, 1);
    // start off with random values in the state
    std::generate(this->state.begin(), this->state.end(),
                  [&] { return this->float_dist(generator); });
    // try out the binding stuff
  }
  /**
   * @brief This function works faster if p is a large number (p > 0.5)
   *
   * @return T sum of the state
   */
  T tick() {
    std::generate(this->trials.begin(), this->trials.end(),
                  [&] { return this->geom_distr(generator); });
    std::sort(this->trials.begin(), this->trials.end());
    const auto uniq = std::unique(this->trials.begin(), this->trials.end());
    // compute the distance of the last unique element
    // this basically takes only the unique elements of the trials
    // because if we repeatedly change the same index, we don't get any
    // advantage
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
    std::generate(this->trials.begin(), this->trials.end(),
                  [&] { return this->geom_distr(generator); });
    for (const auto &t : this->trials) {
      if (t < this->sources) {
        this->state[t] = this->float_dist(generator);
      }
    }
    return this->scale * std::accumulate(state.begin(), state.end(), 0.);
  }
};

// std::mt19937 generator(std::random_device{}());
template <typename T = double> class BufferedAlphaNoise : public NullTicker<T> {
protected:
  std::vector<std::complex<float>> bufferWhite, bufferColoured;
  std::vector<std::complex<float>> bufferWhiteComplex, bufferColouredComplex;
  std::vector<float> result;
  unsigned int bufferSize;
  std::function<float()> gaussPDF;
  T alpha = 1.;
  T scale = 1.;
  std::shared_ptr<kissfft<float>> fwd,
      inv; // configs for forward and inverse real fft
  unsigned int internalCounter = 0;
  unsigned int refills = 0;
  std::mt19937 generator;

public:
  /**
   * @brief Construct a new Buffered Alpha Noise object
   *
   * @param bufferSize the size of the buffer
   * @param alpha the alpha parameter 1/f^alpha
   * @param std the standard deviation of the gaussian
   * @param scale the scaling parameter
   */
  BufferedAlphaNoise(unsigned int bufferSize, T alpha, T std, T scale)
      : bufferSize(bufferSize), alpha(alpha), scale(scale) {
    this->generator = std::mt19937(std::random_device{}());
    this->bufferColoured.resize(2 * bufferSize);
    this->bufferWhite.resize(2 * bufferSize);
    this->result.resize(bufferSize);
    this->bufferColouredComplex.resize(2 * bufferSize);
    this->bufferWhiteComplex.resize(2 * bufferSize);

    // these are the filter weights -- we only have to fill it once per alpha
    // and N
    this->bufferColoured[0] = 1.0;
    for (unsigned int i = 1; i < this->bufferSize; ++i) {
      const float weight = (float)(0.5 * alpha + ((float)(i - 1))) / ((float)i);
      this->bufferColoured[i] = this->bufferColoured[i - 1] * weight;
    }

    this->gaussPDF =
        std::bind(std::normal_distribution<float>(0, std), std::ref(generator));

    this->fwd = std::shared_ptr<kissfft<float>>(
        new kissfft<float>(2 * this->bufferSize, false));
    this->inv = std::shared_ptr<kissfft<float>>(
        new kissfft<float>(2 * this->bufferSize, true));
  }

  ~BufferedAlphaNoise() {}

  void fillBuffer() {
    // this is actual generation function
    // generate random white as a baseline, only N values, rest is 0 padded
    std::generate(this->bufferWhite.begin(),
                  this->bufferWhite.begin() + this->bufferSize, this->gaussPDF);

    for (unsigned int i = this->bufferSize; i < 2 * this->bufferSize; ++i) {
      this->bufferColoured[i] = 0;
      this->bufferWhite[i] = 0;
    }
    // perform the fft
    this->fwd->transform(&this->bufferWhite[0], &this->bufferWhiteComplex[0]);
    this->fwd->transform(&this->bufferColoured[0],
                         &this->bufferColouredComplex[0]);

    // multiply the two
    for (unsigned int i = 0; i < this->bufferSize; ++i) {
      this->bufferColouredComplex[i] =
          this->bufferColouredComplex[i] * this->bufferWhiteComplex[i];
    }
    // invert
    this->bufferColouredComplex[0] =
        this->bufferColouredComplex[0] / std::complex<float>(2.0, 0);
    this->bufferColouredComplex[this->bufferSize - 1] =
        this->bufferColouredComplex[this->bufferSize - 1] /
        std::complex<float>(2.0, 0);
    for (unsigned int i = this->bufferSize; i < 2 * this->bufferSize; ++i) {
      this->bufferColouredComplex[i] = 0.;
    }
    this->inv->transform(&this->bufferColouredComplex[0],
                         &this->bufferWhiteComplex[0]);

    std::transform(this->bufferWhiteComplex.begin(),
                   this->bufferWhiteComplex.begin() + this->bufferSize,
                   this->result.begin(), [&](std::complex<float> x) {
                     return x.real() / (this->bufferSize);
                   });
  }

  const std::vector<float> &getFullBuffer() { return this->result; }

  // overload from null ticker
  T tick() override {
    // we measure only up to a buffer size, not 2x buffer size
    if (this->internalCounter == 0) {
      this->fillBuffer();
    }
    const auto ret = this->result[this->internalCounter];
    this->internalCounter = (this->internalCounter + 1) % this->bufferSize;
    return this->scale * ret;
  }
};

template <typename T = double> class VectorAlphaNoise {
private:
  T scale = 1.;
  // 3 components of type BufferedAlphaNoise, or NullTicker
  std::unique_ptr<NullTicker<T>> components_x, components_y, components_z;
  CVector<T> prevSample, currentSample;
  bool normalized = true;

public:
  VectorAlphaNoise(unsigned int bufferSize, T alpha, T std, T scale,
                   Axis axis = Axis::all)
      : scale(scale) {
    // initialize the as null tickers
    this->components_x = std::unique_ptr<NullTicker<T>>(new NullTicker<T>());
    this->components_y = std::unique_ptr<NullTicker<T>>(new NullTicker<T>());
    this->components_z = std::unique_ptr<NullTicker<T>>(new NullTicker<T>());

    switch (axis) {
    case Axis::all:
      this->components_x = std::unique_ptr<BufferedAlphaNoise<T>>(
          new BufferedAlphaNoise<T>(bufferSize, alpha, std, 1.));
      this->components_y = std::unique_ptr<BufferedAlphaNoise<T>>(
          new BufferedAlphaNoise<T>(bufferSize, alpha, std, 1.));
      this->components_z = std::unique_ptr<BufferedAlphaNoise<T>>(
          new BufferedAlphaNoise<T>(bufferSize, alpha, std, 1.));
      this->normalized = true;
      break;
    case Axis::xaxis:
      this->components_x = std::unique_ptr<BufferedAlphaNoise<T>>(
          new BufferedAlphaNoise<T>(bufferSize, alpha, std, 1.));
      this->normalized = false;
      break;
    case Axis::yaxis:
      this->components_y = std::unique_ptr<BufferedAlphaNoise<T>>(
          new BufferedAlphaNoise<T>(bufferSize, alpha, std, 1.));
      this->normalized = false;
      break;
    case Axis::zaxis:
      this->components_z = std::unique_ptr<BufferedAlphaNoise<T>>(
          new BufferedAlphaNoise<T>(bufferSize, alpha, std, 1.));
      this->normalized = false;
      break;
    default:
      throw std::runtime_error("Invalid axis specified: " +
                               std::to_string(static_cast<int>(axis)));
    }
  }

  CVector<T> tickVector() {
    // TODO  -- if normalized, generate only 2 values and compute the third from
    // the normalization
    this->prevSample = this->currentSample;
    this->currentSample =
        CVector<T>(this->components_x->tick(), this->components_y->tick(),
                   this->components_z->tick());
    if (this->normalized)
      this->currentSample.normalize();
    return this->currentSample * this->scale;
  }

  T tick() { return this->components_x->tick() * this->scale; }

  CVector<T> getPrevSample() { return this->prevSample; }

  T getScale() { return this->scale; }
};

#endif
