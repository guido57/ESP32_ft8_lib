#ifndef FT8_DECODER_API_H
#define FT8_DECODER_API_H

#include <stdbool.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

// Metadata passed to one FT8/FT4 decode call.
typedef struct {
    bool is_ft8;
    float base_freq_mhz;
    struct tm utc;
    double utc_frac_sec;
} ft8_decode_context_t;

// Decode a single audio slot from floating-point samples.
// sample_rate is in Hz, base_freq_mhz is in MHz.
int ft8_decode_slot(const float* signal, int sample_rate, int num_samples, const ft8_decode_context_t* ctx);

#ifdef __cplusplus
}
#endif

#endif
