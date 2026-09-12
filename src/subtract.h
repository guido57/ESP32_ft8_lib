#ifndef SUBTRACT_H
#define SUBTRACT_H

#include <stdlib.h>
#include <vector>
#include <cstddef>
#include <cstdint>

int
blocksize(int rate);

void synthesize(
    float *dst,
    size_t nsamples,
    // const std::vector<int>& re79,
    const uint8_t * tones,
    const std::vector<double>& amps,
    const std::vector<double>& phases,
    double hz0,
    double off_sec,
    double sign,
    int sample_rate);


void subtract(const uint8_t * tones,
              double hz0,
              double hz1,
              double off_sec,
              float * samples_,
              size_t num_samples,
              int sample_rate);

void diagnose_subtraction(
    float *samples,
    size_t num_samples,
    const uint8_t *tones,
    double freq,
    double delay,
    int sample_rate);

float refine_ft8_delay(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print);

    struct complex_amp_t
{
    double amplitude;
    double phase;
    double ci;
    double cq;
};

    complex_amp_t estimate_amplitude_phase(
    float *samples,
    int num_samples,
    uint8_t *tones,
    float delay,
    float freq);


    float refine_ft8_frequency(
    const float *samples,
    int num_samples,
    const uint8_t *tones,
    float delay,
    float freq_coarse);

    float refine_ft8_frequency_costas(
    const float *samples,
    int num_samples,
    const uint8_t *tones,
    float delay,
    float freq_coarse);



#endif
