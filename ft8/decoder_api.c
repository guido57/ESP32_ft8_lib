#include "ft8/decoder_api.h"

#include "common/wave.h"

int ft8_decode_slot(const float* signal, int sample_rate, int num_samples, const ft8_decode_context_t* ctx)
{
    if (signal == NULL || ctx == NULL) {
        return -1;
    }

    return process_buffer(signal,
                          sample_rate,
                          num_samples,
                          ctx->is_ft8,
                          ctx->base_freq_mhz,
                          &ctx->utc,
                          ctx->utc_frac_sec);
}
