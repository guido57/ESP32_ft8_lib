#include <stdio.h>
#include "decoder_api.h"
#include "ft8/ft8_config.h"
#include "common/wave.h"

int ft8_decode_slot(float* signal, int sample_rate, int num_samples, const ft8_decode_context_t* ctx)
{
    if (signal == NULL || ctx == NULL) {
        return -1;
    }

    return process_buffer(signal,
                          sample_rate,
                          num_samples,
                          ctx->is_ft8,
                          ctx->cand_to_subtract,
                          ctx->freq_hz_subtract,
                          ctx->time_delay_subtract,
                          0.0,
                          &ctx->utc,
                          ctx->utc_frac_sec);
                         
                        
}
