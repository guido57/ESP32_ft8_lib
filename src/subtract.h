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

#endif