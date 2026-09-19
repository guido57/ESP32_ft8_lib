// Actual decoding code factored out by KA9Q June 2025
// File/directory handling code is now in main.c

#define _GNU_SOURCE 1
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <libgen.h>
#include <assert.h>
#include <time.h>
#include "ft8/decode.h"
#include "ft8/constants.h"
#include "ft8/ft8_config.h"
#include "subtract.h"

#include "common/wave.h"
#include "common/debug.h"

extern "C" __attribute__((weak)) void ft8_on_message_decoded(
    const char* phase,
    const struct tm* utc,
    double tbase_sec,
    float base_freq_mhz,
    const message_t* msg)
{
  (void)phase;
  (void)utc;
  (void)tbase_sec;
  (void)base_freq_mhz;
  (void)msg;
}

extern "C" __attribute__((weak)) void ft8_on_decode_cycle_complete()
{
}

// Early pass 0 and full residual pass 1 are separate decoder calls. Keep
// callback delivery unique across both calls for the same UTC slot.
static struct tm s_callback_slot = {};
static bool s_callback_slot_valid = false;
static char s_callback_messages[FT8_MAX_DECODED_MSGS][25] = {};
static int s_callback_message_count = 0;

static bool should_emit_decoded_callback(const struct tm* slot,
                                         const char* text)
{
  const bool same_slot = s_callback_slot_valid &&
      s_callback_slot.tm_year == slot->tm_year &&
      s_callback_slot.tm_yday == slot->tm_yday &&
      s_callback_slot.tm_hour == slot->tm_hour &&
      s_callback_slot.tm_min == slot->tm_min &&
      s_callback_slot.tm_sec == slot->tm_sec;
  if (!same_slot) {
    s_callback_slot = *slot;
    s_callback_slot_valid = true;
    s_callback_message_count = 0;
  }

  for (int i = 0; i < s_callback_message_count; ++i) {
    if (strcmp(s_callback_messages[i], text) == 0)
      return false;
  }
  if (s_callback_message_count < FT8_MAX_DECODED_MSGS) {
    strncpy(s_callback_messages[s_callback_message_count], text,
            sizeof(s_callback_messages[s_callback_message_count]) - 1);
    s_callback_messages[s_callback_message_count]
                       [sizeof(s_callback_messages[0]) - 1] = '\0';
    ++s_callback_message_count;
  }
  return true;
}
#include "fft/kiss_fftr.h"
#include "fft/kiss_fft.h"
#if defined(ARDUINO_ARCH_ESP32)
#include <esp_heap_caps.h>
#else
#define heap_caps_malloc(size, caps) malloc(size)   
#define heap_caps_free(ptr) free(ptr)
#endif


#ifndef LOG_LEVEL
#define LOG_LEVEL LOG_DEBUG
#endif

extern "C" void pack_bits(const uint8_t bit_array[], int num_bits, uint8_t packed[]);
extern "C" void ft8_encode(const uint8_t* payload, uint8_t* tones);
extern "C" void kiss_fftr(const kiss_fftr_cfg cfg, const kiss_fft_scalar *timedata, kiss_fft_cpx *freqdata);    

const int kMin_score = FT8_MIN_SCORE; // Minimum sync score threshold for candidates
const int kMax_candidates = FT8_MAX_CANDIDATES; // for 12 kHz sample rate; scaled for other sample rates
const int kLDPC_iterations = FT8_LDPC_ITERATIONS;

#ifndef FT8_DECODE_PASSES
#define FT8_DECODE_PASSES 3
#endif

// This used to be 50. We're now looking at some wider bandwidths *and* FT8 is pretty popular
// Making this bigger seems to only cost memory, which I now allocate from the heap, so what the hell
const int kMax_decoded_messages = FT8_MAX_DECODED_MSGS;

const int kFreq_osr = FT8_FREQ_OSR; // Frequency oversampling rate (bin subdivision)
const int kTime_osr = FT8_TIME_OSR; // Time oversampling rate (symbol subdivision)

static float hann_i(int i, int N)
{
    float x = sinf((float)M_PI * i / N);
    return x * x;
}

static float hamming_i(int i, int N)
{
    const float a0 = (float)25 / 46;
    const float a1 = 1 - a0;

    float x1 = cosf(2 * (float)M_PI * i / N);
    return a0 - a1 * x1;
}

static float blackman_i(int i, int N)
{
    const float alpha = 0.16f; // or 2860/18608
    const float a0 = (1 - alpha) / 2;
    const float a1 = 1.0f / 2;
    const float a2 = alpha / 2;

    float x1 = cosf(2 * (float)M_PI * i / N);
    float x2 = 2 * x1 * x1 - 1; // Use double angle formula

    return a0 - a1 * x1 + a2 * x2;
}

double elapsed_ms(const struct timespec *t0, const struct timespec *t1)
{
  long sec = t1->tv_sec - t0->tv_sec;
  long nsec = t1->tv_nsec - t0->tv_nsec;
  return (double)sec * 1000.0 + (double)nsec / 1000000.0;
}

void waterfall_init(waterfall_t* me, int max_blocks, int num_bins, int time_osr, int freq_osr)
{
    size_t mag_size = max_blocks * time_osr * freq_osr * num_bins * sizeof(me->mag[0]);
    me->max_blocks = max_blocks;
    me->num_blocks = 0;
    me->num_bins = num_bins;
    me->time_osr = time_osr;
    me->freq_osr = freq_osr;
    me->block_stride = (time_osr * freq_osr * num_bins);
    me->mag = (uint8_t  *)malloc(mag_size);
    LOG(LOG_DEBUG, "[ft8] waterfall storage=%zu bytes\n", mag_size);
}

void waterfall_free(waterfall_t* me)
{
    free(me->mag);
}

/// Configuration options for FT4/FT8 monitor
typedef struct
{
    float f_min;             ///< Lower frequency bound for analysis
    float f_max;             ///< Upper frequency bound for analysis
    int sample_rate;         ///< Sample rate in Hertz
    int time_osr;            ///< Number of time subdivisions
    int freq_osr;            ///< Number of frequency subdivisions
    ftx_protocol_t protocol; ///< Protocol: FT4 or FT8
} monitor_config_t;

/// FT4/FT8 monitor object that manages DSP processing of incoming audio data
/// and prepares a waterfall object
typedef struct
{
    float symbol_period; ///< FT4/FT8 symbol period in seconds
    int block_size;      ///< Number of samples per symbol (block)
    int subblock_size;   ///< Analysis shift size (number of samples)
    int nfft;            ///< FFT size
    float fft_norm;      ///< FFT normalization factor
    float* window;       ///< Window function for STFT analysis (nfft samples)
    float* last_frame;   ///< Current STFT analysis frame (nfft samples)
    kiss_fft_scalar* timedata; ///< FFT input scratch buffer (nfft samples)
    kiss_fft_cpx* freqdata;    ///< FFT output scratch buffer (nfft/2+1 samples)
    waterfall_t wf;      ///< Waterfall object
    float max_mag;       ///< Maximum detected magnitude (debug stats)

    // KISS FFT housekeeping variables
    void* fft_work;        ///< Work area required by Kiss FFT
    kiss_fftr_cfg fft_cfg; ///< Kiss FFT housekeeping object
} monitor_t;

// Iinitialize a monitor_t structure based on the provided configuration (monitor_config_t). 
// It calculates DSP parameters, allocates memory for FFT processing and windowing, 
// sets up a waterfall display, and logs key initialization details, 
//preparing the monitor for signal processing tasks.

#if defined ARDUINO_ARCH_ESP32
#include "esp_heap_caps.h"
#endif

void monitor_init(monitor_t* me, const monitor_config_t* cfg, int num_samples)
{
    float slot_time =
        (cfg->protocol == PROTO_FT4) ? FT4_SLOT_TIME : FT8_SLOT_TIME;

    float symbol_period =
        (cfg->protocol == PROTO_FT4) ? FT4_SYMBOL_PERIOD : FT8_SYMBOL_PERIOD;

    me->symbol_period = symbol_period;

    // Samples in one FT8/FT4 symbol
    me->block_size = (int)lroundf(cfg->sample_rate * symbol_period);

    // Time oversampling:
    // 1 -> 2048 samples
    // 2 -> 1024 samples
    // 4 ->  512 samples
    // 8 ->  256 samples
    me->subblock_size = me->block_size / cfg->time_osr;

    // Frequency oversampling:
    // 1 -> 2048-point FFT
    // 2 -> 4096-point FFT
    me->nfft = me->block_size * cfg->freq_osr;

    me->fft_norm = 2.0f / me->nfft;

    LOG(LOG_INFO,
        "[ft8] monitor sample_rate=%d samples=%d time_osr=%d freq_osr=%d\n",
        cfg->sample_rate,
        num_samples,
        cfg->time_osr,
        cfg->freq_osr);

    LOG(LOG_INFO,
        "[ft8] monitor block=%d (%.3fms) subblock=%d (%.3fms) fft=%d (%.3fms)\n",
        me->block_size,
        1000.0f * symbol_period,
        me->subblock_size,
        1000.0f * symbol_period / cfg->time_osr,
        me->nfft,
        1000.0f * me->nfft / cfg->sample_rate);

    /*
     * FFT working buffers in internal DRAM.
     */
    me->window =
        (float *)heap_caps_malloc(
            me->nfft * sizeof(float),
            MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT);

    me->last_frame =
        (float *)heap_caps_malloc(
            me->nfft * sizeof(float),
            MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT);

    me->timedata =
        (kiss_fft_scalar *)heap_caps_malloc(
            me->nfft * sizeof(kiss_fft_scalar),
            MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT);

    me->freqdata =
        (kiss_fft_cpx *)heap_caps_malloc(
            (me->nfft / 2 + 1) * sizeof(kiss_fft_cpx),
            MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT);

        
    if (!me->window || !me->last_frame ||
        !me->timedata || !me->freqdata)
    {
#if defined(ARDUINO_ARCH_ESP32)
        LOG(LOG_ERROR,
            "[ft8] monitor FFT buffers unavailable internal_free=%u largest=%u\n",
            (unsigned)heap_caps_get_free_size(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT),
            (unsigned)heap_caps_get_largest_free_block(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT));
#else
        LOG(LOG_ERROR, "[ft8] monitor FFT buffer allocation failed\n");
#endif
        return;
    }

    /*
     * Start with an empty analysis frame.
     *
     * This is important for freq_osr=2 because the first FFT
     * requires 4096 samples while the input block contains only
     * 2048 samples.
     */
    memset(me->last_frame, 0,
           me->nfft * sizeof(float));

    /*
     * The window is invariant for the lifetime of a monitor.  Include the
     * FFT normalization here, rather than multiplying by it for every input
     * sample of every overlapping FFT frame.
     */
    for (int i = 0; i < me->nfft; ++i)
        me->window[i] = me->fft_norm * hann_i(i, me->nfft);

    /*
     * KISS FFT work area.
     */
    size_t fft_work_size = 0;

    kiss_fftr_alloc(
        me->nfft,
        0,
        NULL,
        &fft_work_size);

    LOG(LOG_DEBUG,
        "[ft8] monitor FFT workspace=%zu bytes\n",
        fft_work_size);

    me->fft_work =
        heap_caps_malloc(
            fft_work_size,
            MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT);

    if (!me->fft_work)
    {
#if defined(ARDUINO_ARCH_ESP32)
        LOG(LOG_ERROR,
            "[ft8] monitor FFT workspace unavailable need=%zu internal_free=%u largest=%u\n",
            fft_work_size,
            (unsigned)heap_caps_get_free_size(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT),
            (unsigned)heap_caps_get_largest_free_block(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT));
#else
        LOG(LOG_ERROR,
            "[ft8] monitor FFT workspace allocation failed\n");
#endif
        return;
    }

    me->fft_cfg =
        kiss_fftr_alloc(
            me->nfft,
            0,
            me->fft_work,
            &fft_work_size);

    /*
     * One waterfall block corresponds to one FT8 symbol period.
     * Each block contains time_osr FFTs.
     */
    const int max_blocks =
        (int)(slot_time / symbol_period);

    const int num_bins =
        (int)(cfg->sample_rate * symbol_period / 2);

    waterfall_init(
        &me->wf,
        max_blocks,
        num_bins,
        cfg->time_osr,
        cfg->freq_osr);

    me->wf.protocol = cfg->protocol;

    me->max_mag = -120.0f;

}

void monitor_free(monitor_t* me)
{
    waterfall_free(&me->wf);

    heap_caps_free(me->fft_work);

    heap_caps_free(me->freqdata);

    heap_caps_free(me->timedata);

    heap_caps_free(me->last_frame);

    heap_caps_free(me->window);
}

#if defined(ARDUINO_ARCH_ESP32)
// The embedded receiver has one finalizer task. Keep its relatively large
// FFT/waterfall workspace for the lifetime of the application instead of
// allocating and freeing it every 15 seconds. Repeated large internal-DRAM
// allocations eventually fail when the rest of the application fragments the
// heap, even though the total free heap still looks adequate.
static monitor_t s_embedded_monitor = {};
static monitor_config_t s_embedded_monitor_cfg = {};
static bool s_embedded_monitor_ready = false;

static bool same_monitor_config(const monitor_config_t* lhs,
                                const monitor_config_t* rhs)
{
    return lhs->f_min == rhs->f_min && lhs->f_max == rhs->f_max &&
           lhs->sample_rate == rhs->sample_rate &&
           lhs->time_osr == rhs->time_osr &&
           lhs->freq_osr == rhs->freq_osr &&
           lhs->protocol == rhs->protocol;
}

static monitor_t* acquire_embedded_monitor(const monitor_config_t* cfg,
                                           int num_samples)
{
    if (s_embedded_monitor_ready &&
        same_monitor_config(&s_embedded_monitor_cfg, cfg))
        return &s_embedded_monitor;

    if (s_embedded_monitor_ready) {
        monitor_free(&s_embedded_monitor);
        memset(&s_embedded_monitor, 0, sizeof(s_embedded_monitor));
        s_embedded_monitor_ready = false;
    }

    monitor_init(&s_embedded_monitor, cfg, num_samples);
    if (s_embedded_monitor.fft_cfg == NULL ||
        s_embedded_monitor.wf.mag == NULL) {
        monitor_free(&s_embedded_monitor);
        memset(&s_embedded_monitor, 0, sizeof(s_embedded_monitor));
        return NULL;
    }

    s_embedded_monitor_cfg = *cfg;
    s_embedded_monitor_ready = true;
    LOG(LOG_INFO,
        "[ft8] monitor workspace ready internal_free=%u largest=%u\n",
        (unsigned)heap_caps_get_free_size(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT),
        (unsigned)heap_caps_get_largest_free_block(MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT));
    return &s_embedded_monitor;
}
#endif

// Funzione helper veloce per approssimare (2 * 10 * log10f(x) + 240) usando i cicli di clock hardware.
// Evita del tutto la pesantissima funzione log10f() software.
static inline int fast_db_scale(float mag2) 
{
    // Convertiamo il float in un intero interpretando i bit (Type Punning)
    // Questo ci permette di estrarre l'esponente binario direttamente in 1-2 cicli di clock
    union { float f; uint32_t i; } u;
    u.f = mag2;
    
    // Se il valore è infinitesimo o zero, restituiamo il minimo del range clampato (0 dB o -120dB scaled)
    if (u.i < 0x31800000) return 0; // Sotto ~1E-9

    // Estraiamo l'esponente del float (bit 23-30)
    int exponent = (int)((u.i >> 23) & 0xFF) - 127;
    
    // Estraiamo la mantissa come frazione lineare
    float mantissa = (float)(u.i & 0x7FFFFF) / 8388608.0f;

    // log2(x) = esponente + mantissa (approssimazione lineare ottima per i dB)
    float log2_approx = (float)exponent + mantissa;

    // Convertiamo da log2 a log10 per la formula dei dB:
    // db = 10 * log10(x) = 10 * log2(x) * log10(2) = 10 * log2(x) * 0.30103f = 3.0103f * log2(x)
    // La formula originale fa: scaled = (int)(2 * db + 240) = (int)(2 * 3.0103f * log2_approx + 240)
    int scaled = (int)(6.0206f * log2_approx + 240.0f);
    
    return scaled;
}
#include <string.h> // Necessario per memmove

#include <string.h>
#include <math.h>

static inline float fast_log2f(float x)
{
    union {
        float f;
        uint32_t i;
    } v = { x };

    int e = (int)(v.i >> 23) - 127;

    v.i = (v.i & 0x7FFFFF) | 0x3F800000;

    float m = v.f - 1.0f;

    return (float)e +
           0.00903028f +
           m * (1.3211363f - 0.33688046f * m);
}

void monitor_process(monitor_t* me, const float* frame)
{
    struct timespec t0, t_sub_start, t_input, t_fft, t_mag, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    if (me->wf.num_blocks >= me->wf.max_blocks)
        return;

    if (me->timedata == NULL ||
        me->freqdata == NULL ||
        me->last_frame == NULL)
        return;

    const int nfft         = me->nfft;
    const int block_size   = me->block_size;
    const int subblock     = me->subblock_size;
    const int time_osr     = me->wf.time_osr;
    const int freq_osr     = me->wf.freq_osr;
    const int num_bins     = me->wf.num_bins;

    /*
     * One monitor_process() call receives exactly one symbol
     * worth of new samples.
     */
    int frame_pos = 0;

    int offset =
        me->wf.num_blocks * me->wf.block_stride;

    float local_max_mag = me->max_mag;
    double input_ms = 0.0;
    double fft_ms = 0.0;
    double mag_ms = 0.0;

    float * __restrict last_frame = me->last_frame;
    float * __restrict timedata   = me->timedata;
    const float * __restrict window = me->window;
    uint8_t * __restrict wf_mag = me->wf.mag;

    typedef struct {
        float r;
        float i;
    } cpx_t;

    const cpx_t *freqdata =
        (const cpx_t *)me->freqdata;

    for (int time_sub = 0;
         time_sub < time_osr;
         ++time_sub)
    {
        clock_gettime(CLOCK_MONOTONIC, &t_sub_start);
        /*
         * Shift the previous FFT frame left by subblock samples.
         *
         * Example:
         *
         * nfft=4096, subblock=1024
         *
         * old:
         * [---------------- 4096 ----------------]
         *
         * new:
         * [---------- 3072 ------------][--1024--]
         *
         * The last 1024 samples are replaced by the
         * new input samples.
         */
        const int keep = nfft - subblock;

        if (keep > 0)
        {
            memmove(
                last_frame,
                last_frame + subblock,
                keep * sizeof(float));
        }

        memcpy(
            last_frame + keep,
            frame + frame_pos,
            subblock * sizeof(float));

        frame_pos += subblock;

        /*
         * Window and normalize.
         *
         * Native high-resolution mode deliberately keeps the analysis window
         * to one FT8 symbol and zero-pads it to the requested FFT length.
         * Frequency oversampling then interpolates the spectrum without
         * mixing several different FSK symbols into one FFT.
         */
#if defined(NATIVE_BUILD) && defined(FT8_NATIVE_SYMBOL_WINDOW)
        memset(timedata, 0, nfft * sizeof(*timedata));
        const float *const symbol = last_frame + nfft - block_size;
        for (int i = 0; i < block_size; ++i)
            timedata[i] = me->fft_norm * hann_i(i, block_size) * symbol[i];
#else
        for (int i = 0; i < nfft; ++i)
        {
            timedata[i] = window[i] * last_frame[i];
        }
#endif

        clock_gettime(CLOCK_MONOTONIC, &t_input);

        /*
         * FFT.
         */
#if defined ARDUINO_ARCH_ESP32
        // Benchmark normalization: subtract-test uses KISS because its
        // ESP-DSP transform corrupts the time_osr=16 waterfall.  Use the
        // identical, correctness-preserving backend here for a fair P0/P1
        // comparison.
        kiss_fftr(
            me->fft_cfg,
            timedata,
            me->freqdata);


#else

        kiss_fftr(
            me->fft_cfg,
            timedata,
            me->freqdata);

#endif

        clock_gettime(CLOCK_MONOTONIC, &t_fft);

        /*
         * Convert FFT bins into waterfall magnitudes.
         *
         * freq_osr=1:
         *
         *   FFT bins:
         *   0,1,2,...2047
         *
         * freq_osr=2:
         *
         *   FFT bins:
         *   0,1,2,...4095
         *
         * We keep the two interleaved frequency grids.
         */
        for (int freq_sub = 0; freq_sub < freq_osr; ++freq_sub)
        {
            for (int bin = 0;
                 bin < num_bins;
                 ++bin)
            {
                const int src_bin =
                    bin * freq_osr + freq_sub;

                const float real_part =
                    freqdata[src_bin].r;

                const float imag_part =
                    freqdata[src_bin].i;

                const float mag2 =
                    real_part * real_part +
                    imag_part * imag_part;

                const float db =
                    3.0102999566f *
                    fast_log2f(1E-12f + mag2);

                int scaled =
                    (int)(2.0f * db + 240.0f);

                if (scaled < 0)
                    scaled = 0;
                else if (scaled > 255)
                    scaled = 255;

                wf_mag[offset++] =
                    (uint8_t)scaled;

                if (db > local_max_mag)
                    local_max_mag = db;
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &t_mag);
        input_ms += elapsed_ms(&t_sub_start, &t_input);
        fft_ms += elapsed_ms(&t_input, &t_fft);
        mag_ms += elapsed_ms(&t_fft, &t_mag);
    }

    me->max_mag = local_max_mag;

    /*
     * One complete 160-ms block has now been added.
     */
    ++me->wf.num_blocks;

    clock_gettime(CLOCK_MONOTONIC, &t1);

    // LOG(LOG_DEBUG,
    //     "[ft8] monitor_process: %.3f ms (input %.3f, FFT %.3f, mag %.3f)\n",
    //     elapsed_ms(&t0, &t1), input_ms, fft_ms, mag_ms);
}


void monitor_reset(monitor_t* me)
{
    me->wf.num_blocks = 0;
    me->max_mag = -120.0f;
    if (me->last_frame != NULL)
        memset(me->last_frame, 0, me->nfft * sizeof(*me->last_frame));
}

// Used to sort messages by ascending frequency, and to push empty entries to end
int mcompare(void const *a, void const *b){
  message_t const *ma = *(message_t const **)a;
  message_t const *mb = *(message_t const **)b;
  // Null entries go to end of list
  if(ma == NULL && mb == NULL)
    return 0;
  else if(ma == NULL)
    return +1;
  else if(mb == NULL)
    return -1;
  if(ma->freq_hz > mb->freq_hz)
    return +1;
  else if(ma->freq_hz < mb->freq_hz)
    return -1;
  return 0;
}

static int decode_clock_seconds(const struct tm* slot_utc, float time_sec)
{
  // A candidate's time is measured from the start of the FT8/FT4 slot.
  // Work directly in UTC clock seconds so the host's local time zone cannot
  // alter the displayed value.
  int clock_seconds = slot_utc->tm_hour * 3600 + slot_utc->tm_min * 60 +
                      slot_utc->tm_sec + (int)lroundf(time_sec);
  clock_seconds %= 24 * 3600;
  if (clock_seconds < 0)
    clock_seconds += 24 * 3600;
  return clock_seconds;
}

static float mag_u8_to_power(uint8_t mag)
{
  // In monitor_process(): mag ~= 2*db + 240, where db = 10*log10(power)
  float db = 0.5f * ((float)mag - 240.0f);
  return powf(10.0f, db / 10.0f);
}

static float estimate_global_noise_power(const waterfall_t *wf)
{
  int hist[256] = {0};
  int total = wf->num_blocks * wf->block_stride;
  if (total <= 0)
    return 1e-12f;

  const uint8_t *p = wf->mag;
  for (int i = 0; i < total; ++i)
    ++hist[p[i]];

  // Use a lower percentile as robust noise-floor proxy, avoiding signal peaks.
  int target = (int)(0.35f * total);
  int acc = 0;
  int mag_floor = 0;
  for (int m = 0; m < 256; ++m)
  {
    acc += hist[m];
    if (acc >= target)
    {
      mag_floor = m;
      break;
    }
  }

  return mag_u8_to_power((uint8_t)mag_floor);
}


/**
 * Estimate the SNR of a candidate signal in dB, referenced to a 2500 Hz bandwidth.
 *
 * wf: the waterfall structure containing the signal magnitudes
 * cand: the candidate signal for which to estimate the SNR
 * plain: the decoded plain bits corresponding to the candidate
 * noise_power: the estimated global noise power
 * Returns the estimated SNR in dB, referenced to a 2500 Hz bandwidth.
 */
static float estimate_candidate_snr_db_2500(const waterfall_t *wf, const candidate_t *cand, const uint8_t *plain, float noise_power)
{
  float sum_sig = 0.0f;
  int n_sig = 0;
  int base = (cand->time_offset * wf->time_osr) + cand->time_sub;
  base = (base * wf->freq_osr) + cand->freq_sub;
  base = (base * wf->num_bins) + cand->freq_offset;

  if (wf->protocol == PROTO_FT8)
  {
    for (int k = 0; k < FT8_ND; ++k)
    {
      int sym_idx = k + ((k < 29) ? 7 : 14);
      int block_abs = cand->time_offset + sym_idx;
      if (block_abs < 0 || block_abs >= wf->num_blocks)
        continue;

      int bit_idx = 3 * k;
      int g = ((plain[bit_idx + 0] & 1) << 2) |
              ((plain[bit_idx + 1] & 1) << 1) |
              ((plain[bit_idx + 2] & 1) << 0);
      int tone = kFT8_Gray_map[g & 7];

      const uint8_t *p = wf->mag + base + (sym_idx * wf->block_stride);
      sum_sig += mag_u8_to_power(p[tone]);
      ++n_sig;
    }
  }
  else
  {
    for (int k = 0; k < FT4_ND; ++k)
    {
      int sym_idx = k + ((k < 29) ? 5 : ((k < 58) ? 9 : 13));
      int block_abs = cand->time_offset + sym_idx;
      if (block_abs < 0 || block_abs >= wf->num_blocks)
        continue;

      int bit_idx = 2 * k;
      int g = ((plain[bit_idx + 0] & 1) << 1) |
              ((plain[bit_idx + 1] & 1) << 0);
      int tone = kFT4_Gray_map[g & 3];

      const uint8_t *p = wf->mag + base + (sym_idx * wf->block_stride);
      sum_sig += mag_u8_to_power(p[tone]);
      ++n_sig;
    }
  }

  if (n_sig == 0)
    return -99.0f;

  float sig = sum_sig / n_sig;
  float noise = (noise_power > 0.0f) ? noise_power : 1e-12f;
  if (noise <= 0.0f)
    return -99.0f;

  // Convert local tone-group SNR to WSJT-style 2500 Hz reference.
  // Use per-tone bin bandwidth (~6.25 Hz for FT modes), then reference to 2500 Hz.
  float snr_bin_db = 10.0f * log10f(sig / noise + 1e-12f);
  float snr_2500_db = snr_bin_db - 10.0f * log10f(2500.0f / 6.25f);
  return snr_2500_db;
}

static void
get_ft8_tones_from_plain174(const uint8_t *plain174,
                            uint8_t *tones)
{
    uint8_t a91[FTX_LDPC_K_BYTES];

    pack_bits(plain174, FTX_LDPC_K, a91);

    uint8_t payload[10];
    memcpy(payload, a91, 10);

    ft8_encode(payload, tones);
}

float rms(const float *samples, size_t num_samples)
{
    float sum2 = 0.0f;

    for (size_t i = 0; i < num_samples; ++i)
    {
        const float x = samples[i];
        sum2 += x * x;
    }

    return sqrtf(sum2 / (float)num_samples);
}

// ------------------------------------------------------------------------------------
// Process "num_samples" samples stored in "signal"
// ------------------------------------------------------------------------------------
int process_buffer(float *samples,int sample_rate, int num_samples, 
    bool is_ft8, int cand_to_subtract, float freq_hz_subtract, float time_delay_subtract, float base_freq, struct tm const *tmp, double sec, int max_decode_passes, bool is_early_pass, int max_candidates){
  assert(samples != NULL && tmp != NULL);

#if defined(NATIVE_BUILD) && !defined(FT8_NATIVE_ALTERNATE_CANCEL)
#define FT8_NATIVE_ALTERNATE_CANCEL 0
#endif
#if defined(NATIVE_BUILD) && FT8_NATIVE_ALTERNATE_CANCEL
  float *const native_original = (float *)malloc((size_t)num_samples * sizeof(*samples));
  if (native_original == NULL)
    return -1;
  memcpy(native_original, samples, (size_t)num_samples * sizeof(*samples));
  const double primary_subtract_ramp = subtract_ramp;
#endif


  // Compute Waterfall accumulation (FFT)
  monitor_config_t const mon_cfg = {
        .f_min = 100,
        .f_max = (sample_rate)/2.0f - 500.0f, // allow room for the receiver filter rolloff
        .sample_rate = sample_rate,
        .time_osr = FT8_TIME_OSR,
        .freq_osr = FT8_FREQ_OSR,
        .protocol = is_ft8 ? PROTO_FT8 : PROTO_FT4
  };

  // Hash table for decoded messages (to check for duplicates)
  int num_decoded = 0;
  // Pointer to kMax_decoded_messages-element array of message_t structures
  message_t *decoded; 
  // Pointer to kMax_decoded_messsages-element array of pointers to message_t structures
  message_t **decoded_hashtable;
  // Pointer to kMax_decoded_messages-element array of message_t structures
  decoded = (message_t *) calloc(sizeof(message_t), kMax_decoded_messages);
  // Pointer to kMax_decoded_messsages-element array of pointers to message_t structures
  decoded_hashtable = (message_t **) calloc(sizeof(message_t *), kMax_decoded_messages);

#if defined(ARDUINO_ARCH_ESP32)
  monitor_t* mon_ptr = acquire_embedded_monitor(&mon_cfg, num_samples);
  if (mon_ptr == NULL) {
    LOG(LOG_ERROR, "[ft8] monitor initialization failed\n");
    free(decoded);
    free(decoded_hashtable);
    return -1;
  }
  monitor_t& mon = *mon_ptr;
#else
  monitor_t mon = {0};
#endif

  // Initial waterfall and DDC-residual decodes share one result set.  This
  // prevents a residual alias from being reported as an additional message.
  auto store_decoded_message = [&](const message_t& candidate) {
    int idx_hash = candidate.hash % kMax_decoded_messages;
    for (int probes = 0; probes < kMax_decoded_messages; ++probes) {
      if (decoded_hashtable[idx_hash] == NULL) {
        decoded[idx_hash] = candidate;
        decoded_hashtable[idx_hash] = &decoded[idx_hash];
        ++num_decoded;
        return true;
      }
      if (decoded_hashtable[idx_hash]->hash == candidate.hash &&
          0 == strcmp(decoded_hashtable[idx_hash]->text, candidate.text)) {
        return false;
      }
      idx_hash = (idx_hash + 1) % kMax_decoded_messages;
    }
    LOG(LOG_WARN, "[ft8] decoded-message table full; dropping text=%s\n",
        candidate.text);
    return false;
  };

  // The DDC windows overlap heavily.  Locate a message before doing any
  // optional, expensive presentation-only work (notably the local SNR
  // estimator), while keeping insertion and duplicate handling in the one
  // store function above.
  auto find_decoded_message = [&](const message_t& candidate) -> message_t* {
    int idx_hash = candidate.hash % kMax_decoded_messages;
    for (int probes = 0; probes < kMax_decoded_messages; ++probes) {
      if (decoded_hashtable[idx_hash] == NULL)
        return NULL;
      if (decoded_hashtable[idx_hash]->hash == candidate.hash &&
          0 == strcmp(decoded_hashtable[idx_hash]->text, candidate.text))
        return decoded_hashtable[idx_hash];
      idx_hash = (idx_hash + 1) % kMax_decoded_messages;
    }
    return NULL;
  };

  // ---------------------------------------------------
  //  ONE PASS: standard waterfall decode followed by local DDC-SIC.
  //  Do not re-run a waterfall on a conventionally subtracted waveform.
  // ---------------------------------------------------
  const int primary_passes =
      (max_decode_passes > 0 && max_decode_passes < FT8_DECODE_PASSES)
          ? max_decode_passes
          : FT8_DECODE_PASSES;
  const int total_passes = primary_passes
#if defined(NATIVE_BUILD) && FT8_NATIVE_ALTERNATE_CANCEL
      * 2
#endif
      ;

  // This allocation needs a contiguous internal-DRAM block.  Reusing one
  // monitor across passes avoids a second allocation failing after pass-0
  // subtraction and other temporary allocations have fragmented the heap.
#if !defined(ARDUINO_ARCH_ESP32)
  monitor_init(&mon, &mon_cfg, num_samples);
  if (mon.fft_cfg == NULL || mon.wf.mag == NULL) {
    LOG(LOG_ERROR, "[ft8] monitor initialization failed\n");
    monitor_free(&mon);
    free(decoded);
    free(decoded_hashtable);
#if defined(NATIVE_BUILD) && FT8_NATIVE_ALTERNATE_CANCEL
    subtract_ramp = primary_subtract_ramp;
    free(native_original);
#endif
    return -1;
  }
#endif

  for(int pass=0; pass<total_passes; ++pass){
    const int decoded_before_pass = num_decoded;
#if defined(NATIVE_BUILD) && FT8_NATIVE_ALTERNATE_CANCEL
    const bool alternate_cancel = pass >= primary_passes;
    const int branch_pass = pass % primary_passes;
    if (pass == primary_passes) {
      memcpy(samples, native_original, (size_t)num_samples * sizeof(*samples));
      subtract_ramp = 0.22;
    }
#else
    const bool alternate_cancel = false;
    const int branch_pass = pass;
#endif
#if defined(NATIVE_DIAGNOSTIC)
    printf("========================================\n");
    printf("Starting pass %d\n", pass);
    printf("========================================\n");
#endif
    
    struct timespec t_wf0 = {0};
    struct timespec t_wf1 = {0};
    struct timespec t_dec0 = {0};
    struct timespec t_dec1 = {0};
    struct timespec t_ndec = {0};
    struct timespec t_sub0 = {0};
    struct timespec t_sub1 = {0};
    struct timespec t_cand0 = {0};
    struct timespec t_cand1 = {0};
    
    clock_gettime(CLOCK_MONOTONIC, &t_wf0);

    // Rebuild the waterfall from the current (possibly subtracted) waveform
    // without reallocating the monitor's FFT storage.
    monitor_reset(&mon);
    LOG(LOG_DEBUG, "[ft8] stage=%s pass=%d waterfall start blocks=%d block_samples=%d\n",
        is_early_pass ? "early" : "full", pass,
        mon.wf.max_blocks, mon.block_size);

    for (int frame_pos = 0; frame_pos + mon.block_size <= num_samples; frame_pos += mon.block_size){
        // Process the waveform data frame by frame - you could have a live loop here with data from an audio device
        // (cool, now that we can get sample timings - KA9Q)
        monitor_process(&mon, samples + frame_pos);
    }

    clock_gettime(CLOCK_MONOTONIC, &t_wf1);
    LOG(LOG_INFO, "[ft8] stage=%s pass=%d waterfall blocks=%d elapsed=%.1fms max=%.1fdB\n",
        is_early_pass ? "early" : "full", pass, mon.wf.num_blocks,
        elapsed_ms(&t_wf0, &t_wf1), mon.max_mag);
    
    float const noise_power = estimate_global_noise_power(&mon.wf);
    
    // Find top candidates by Costas sync score and localize them in time and frequency
    int candidate_size = (mon_cfg.f_max * kMax_candidates) / 3000; // Scale by bandwidth relative to the original 3 kHz
    if (max_candidates > 0 && candidate_size > max_candidates)
      candidate_size = max_candidates;
    candidate_t candidate_list[candidate_size];
    int num_candidates = ft8_find_sync(&mon.wf, candidate_size, candidate_list, kMin_score);
    clock_gettime(CLOCK_MONOTONIC, &t_dec0);
    LOG(LOG_INFO, "[ft8] stage=%s pass=%d candidates=%d threshold=%d search=%.1fms\n",
        is_early_pass ? "early" : "full", pass, num_candidates, kMin_score,
        elapsed_ms(&t_wf1, &t_dec0));

    // ==================================================
    // Go over candidates and attempt to decode messages
    // ==================================================
    for (int idx = 0; idx < num_candidates; ++idx){
        clock_gettime(CLOCK_MONOTONIC, &t_cand0);
        candidate_t localized_candidate = candidate_list[idx];
        const candidate_t* cand = &localized_candidate;
        if (cand->score < kMin_score)
            continue;

        float const freq_hz = (cand->freq_offset + (float)cand->freq_sub / mon.wf.freq_osr) / mon.symbol_period;
        float const time_sec = (cand->time_offset + (float)cand->time_sub / mon.wf.time_osr) * mon.symbol_period;
        message_t message = {0};      // Written by ft8_decode()
        decode_status_t status = {0}; // ditto
        uint8_t plain174[FTX_LDPC_N];
        const bool decoded_ok = ft8_decode(
            &mon.wf, cand, &message, kLDPC_iterations, &status, plain174);
        if (!decoded_ok){
            if (status.ldpc_errors > 0){
            LOG(LOG_DEBUG, "[ft8] candidate=%d LDPC_errors=%d freq=%.1fHz time=%.3fs\n",
                idx, status.ldpc_errors, freq_hz, time_sec);
            }else if (status.crc_calculated != status.crc_extracted){
            LOG(LOG_DEBUG, "[ft8] candidate=%d CRC mismatch\n", idx);
            }else if (status.unpack_status != 0){
            LOG(LOG_DEBUG, "[ft8] candidate=%d unpack failed\n", idx);
            }
            clock_gettime(CLOCK_MONOTONIC, &t_ndec);
            LOG(LOG_DEBUG, "[ft8] candidate=%d decode failed elapsed=%.1fms\n",
                idx, elapsed_ms(&t_cand0, &t_ndec));
            continue;
        }
        
        message.freq_hz = freq_hz; // Save so we can sort on it and display it
        message.time_sec = time_sec; // Time offset of start from nominal UTC :00/:15/:30/:45 or :00/:07.5/:15/...
        message.score = cand->score;
        message.snr_raw_db = estimate_candidate_snr_db_2500(&mon.wf, cand, plain174, noise_power);
        {
            // Affine calibration from raw estimator to WSJT-like displayed SNR.
            float snr = FT8_SNR_RAW_SCALE * message.snr_raw_db + FT8_SNR_SCORE_SCALE * (float)message.score + FT8_SNR_OFFSET;
            if (snr < FT8_SNR_MIN_DB)
            snr = FT8_SNR_MIN_DB;
            if (snr > FT8_SNR_MAX_DB)
            snr = FT8_SNR_MAX_DB;
            message.snr_db = snr;
        }

        LOG(LOG_DEBUG, "[ft8] candidate=%d dedup score=%d time=%.3fs freq=%.1fHz\n",
            idx, cand->score, time_sec, freq_hz);
        const bool newly_decoded = store_decoded_message(message);
        if (newly_decoded){
            clock_gettime(CLOCK_MONOTONIC, &t_dec1);
            if (status.osd_score > 0.0f)
                LOG(LOG_DEBUG, "[ft8] OSD accepted candidate=%d metric=%.3f\n", idx, status.osd_score);
            LOG(LOG_DEBUG,"[ft8] decoded candidate=%d score=%d metric=%.3f snr=%+.1f dt=%+.3fs freq=%.1fHz text=%s decode=%.1fms\n",
                idx, cand->score, status.codeword_metric, message.snr_db, time_sec-0.14f, freq_hz,
                message.text, elapsed_ms(&t_cand0, &t_dec1));

            // Deliver each unique decode immediately. The cycle-complete
            // callback below still flushes batched consumers only once after
            // all decode/subtraction passes have finished.
            if (should_emit_decoded_callback(tmp, message.text)) {
              ft8_on_message_decoded(is_early_pass ? "early" : "final",
                                     tmp, sec, base_freq, &message);
            }

            // A caller-supplied position overrides fine synchronization for
            // one pass-0 candidate only.  It must not disable normal SIC for
            // all the other messages, nor be reapplied to an unrelated
            // candidate with the same index in a later residual pass.
            const bool manual_coordinates =
                cand_to_subtract >= 0 && pass == 0 && idx == cand_to_subtract;
            const bool subtract_this_candidate = newly_decoded || alternate_cancel;

            // On the alternate native branch, decode duplicates must still
            // be cancelled because this branch starts from the original slot.
            if (subtract_this_candidate &&
                (branch_pass < primary_passes - 1)){
                uint8_t tones[79];
                get_ft8_tones_from_plain174(plain174, tones);

                float subtract_freq;
                float subtract_delay;
                if (manual_coordinates) {
                    subtract_freq = freq_hz_subtract;
                    subtract_delay = time_delay_subtract;
                    LOG(LOG_DEBUG,
                        "[ft8] manual subtract candidate=%d freq=%.3fHz delay=%.6fs\n",
                        idx, subtract_freq, subtract_delay);
                } else {
                    // Refine the coarse waterfall hypothesis against the
                    // current residual before synthesizing and subtracting it.
                    refine_ft8_joint(samples, num_samples, sample_rate, tones,
                                     time_sec - 0.14f, freq_hz, idx,
                                     &subtract_delay, &subtract_freq);
                }

                subtract(tones,
                    subtract_freq,
                    subtract_freq,
                    subtract_delay,
                    samples,
                    num_samples,
                    sample_rate);
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &t_cand1);
        LOG(LOG_DEBUG, "[ft8] candidate=%d complete elapsed=%.1fms\n",
            idx, elapsed_ms(&t_cand0, &t_cand1));
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t_dec1);
    LOG(LOG_INFO,
        "[ft8] stage=%s pass=%d candidates=%d new=%d total=%d decode=%.1fms\n",
        is_early_pass ? "early" : "full", pass, num_candidates,
        num_decoded - decoded_before_pass, num_decoded,
        elapsed_ms(&t_dec0, &t_dec1));
    

  } // End of pass

#if !defined(ARDUINO_ARCH_ESP32)
  monitor_free(&mon);
#endif

  // The result table is shared by both passes, so print once after all
  // conventional and residual decodes have been deduplicated.
  qsort(decoded_hashtable, kMax_decoded_messages, sizeof *decoded_hashtable, mcompare);
  LOG(LOG_INFO, "[ft8] results stage=%s slot=%02d:%02d:%02d.%03d decoded=%d\n",
      is_early_pass ? "early" : "full",
      tmp->tm_hour, tmp->tm_min, tmp->tm_sec,
      (int)lround(sec * 1000.0), num_decoded);
  for (int i = 0; i < num_decoded; ++i) {
    message_t const* mp = decoded_hashtable[i];
    if (mp == NULL)
      continue;
    const int clock_seconds = decode_clock_seconds(tmp, mp->time_sec);
    LOG_DECODED(
        "[ft8] %02d%02d%02d snr=%+d dt=%+.2fs freq=%.1fHz text=%s%s\n",
        clock_seconds / 3600, (clock_seconds / 60) % 60,
        clock_seconds % 60, (int)lroundf(mp->snr_db),
        mp->time_sec - 0.14f, mp->freq_hz, mp->text, log_color_reset());
  }
  if (!is_early_pass)
    ft8_on_decode_cycle_complete();
  fflush(stderr);

  free(decoded);
  free(decoded_hashtable);
#if defined(NATIVE_BUILD) && FT8_NATIVE_ALTERNATE_CANCEL
  subtract_ramp = primary_subtract_ramp;
  free(native_original);
#endif
  return num_decoded; // Caller frees signal
}
