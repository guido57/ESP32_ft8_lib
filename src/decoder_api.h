#ifndef FT8_DECODER_API_H
#define FT8_DECODER_API_H

#include <stdbool.h>
#include <stdint.h>
#include <time.h>

#ifdef __cplusplus
// extern "C" {
#endif

// Metadata passed to one FT8/FT4 decode call.
typedef struct {
    bool is_ft8;
    float base_freq_mhz;
    struct tm utc;
    double utc_frac_sec;
    int cand_to_subtract;       // Candidate to subtract
    float freq_hz_subtract;       // Frequency in Hz to subtract
    float time_delay_subtract; // Time delay in seconds
    // 0 uses the build-time FT8_DECODE_PASSES setting. A positive value may
    // cap the pass count so a real-time consumer can catch up under load.
    int max_decode_passes;
    // Internal pipeline controls used by the embedded early pass.
    bool is_early_pass;
    int max_candidates;
} ft8_decode_context_t;

// Decode a single audio slot from floating-point samples.
// sample_rate is in Hz, base_freq_mhz is in MHz.
int ft8_decode_slot(float* signal, int sample_rate, int num_samples, const ft8_decode_context_t* ctx);

// Opaque, slot-sized capture buffer.  Samples are retained so the final
// decoder can perform residual/wide subtraction over the whole slot.
typedef struct ft8_stream_decoder ft8_stream_decoder_t;

ft8_stream_decoder_t* ft8_stream_open(int sample_rate, const ft8_decode_context_t* ctx);
int ft8_stream_append_i16(ft8_stream_decoder_t* stream, const int16_t* signal, int num_samples);
void ft8_stream_set_max_decode_passes(ft8_stream_decoder_t* stream, int max_decode_passes);
ft8_stream_decoder_t* ft8_stream_snapshot_early(const ft8_stream_decoder_t* stream);
int ft8_stream_finalize(ft8_stream_decoder_t* stream);
void ft8_stream_close(ft8_stream_decoder_t* stream);

#ifdef __cplusplus
// }
#endif

#endif
