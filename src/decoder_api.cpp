#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "decoder_api.h"
#include "ft8/ft8_config.h"
#include "common/wave.h"

#if defined(ARDUINO_ARCH_ESP32)
#include <esp_heap_caps.h>
#endif

struct ft8_stream_decoder {
    float* samples;
    int sample_rate;
    int capacity;
    int count;
    ft8_decode_context_t context;
};

static float* allocate_slot_samples(size_t count)
{
#if defined(ARDUINO_ARCH_ESP32)
    // The residual decoder needs a mutable full-slot waveform.  Keep it out
    // of scarce internal RAM on PSRAM-equipped ESP32-S3 boards.
    float* samples = static_cast<float*>(heap_caps_malloc(
        count * sizeof(*samples), MALLOC_CAP_SPIRAM | MALLOC_CAP_8BIT));
    if (samples != NULL)
        return samples;
#endif
    return static_cast<float*>(malloc(count * sizeof(float)));
}

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
                          ctx->base_freq_mhz,
                          &ctx->utc,
                          ctx->utc_frac_sec,
                          ctx->max_decode_passes,
                          ctx->is_early_pass,
                          ctx->max_candidates);
}

ft8_stream_decoder_t* ft8_stream_open(int sample_rate, const ft8_decode_context_t* ctx)
{
    if (sample_rate <= 0 || ctx == NULL)
        return NULL;

    const double slot_seconds = ctx->is_ft8 ? 15.0 : 7.5;
    const int capacity = (int)(slot_seconds * sample_rate + sample_rate + 0.5);
    ft8_stream_decoder_t* stream = static_cast<ft8_stream_decoder_t*>(calloc(1, sizeof(*stream)));
    if (stream == NULL)
        return NULL;

    stream->samples = allocate_slot_samples((size_t)capacity);
    if (stream->samples == NULL) {
        free(stream);
        return NULL;
    }

    stream->capacity = capacity;
    stream->sample_rate = sample_rate;
    stream->context = *ctx;
    return stream;
}

int ft8_stream_append_i16(ft8_stream_decoder_t* stream, const int16_t* signal, int num_samples)
{
    if (stream == NULL || signal == NULL || num_samples < 0)
        return -1;

    const int available = stream->capacity - stream->count;
    const int accepted = (num_samples < available) ? num_samples : available;
    for (int i = 0; i < accepted; ++i)
        stream->samples[stream->count + i] = (float)signal[i] * (1.0f / 32768.0f);
    stream->count += accepted;
    return accepted;
}

void ft8_stream_set_max_decode_passes(ft8_stream_decoder_t* stream, int max_decode_passes)
{
    if (stream != NULL)
        stream->context.max_decode_passes = max_decode_passes;
}

ft8_stream_decoder_t* ft8_stream_snapshot_early(
    const ft8_stream_decoder_t* stream)
{
    if (stream == NULL || stream->count == 0)
        return NULL;

    ft8_stream_decoder_t* early =
        static_cast<ft8_stream_decoder_t*>(calloc(1, sizeof(*early)));
    if (early == NULL)
        return NULL;

    early->samples = allocate_slot_samples((size_t)stream->count);
    if (early->samples == NULL) {
        free(early);
        return NULL;
    }

    memcpy(early->samples, stream->samples,
           (size_t)stream->count * sizeof(*early->samples));
    early->capacity = stream->count;
    early->count = stream->count;
    early->sample_rate = stream->sample_rate;
    early->context = stream->context;
    early->context.max_decode_passes = 1;
    early->context.is_early_pass = true;
    early->context.max_candidates = 6;
    return early;
}

int ft8_stream_finalize(ft8_stream_decoder_t* stream)
{
    if (stream == NULL || stream->count == 0)
        return -1;
    return ft8_decode_slot(stream->samples, stream->sample_rate,
                           stream->count, &stream->context);
}

void ft8_stream_close(ft8_stream_decoder_t* stream)
{
    if (stream == NULL)
        return;
    free(stream->samples);
    free(stream);
}
