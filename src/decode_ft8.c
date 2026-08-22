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
#include <esp_heap_caps.h>

#include "ft8/decode.h"
#include "ft8/constants.h"
#include "ft8/ft8_config.h"

#include "common/wave.h"
#include "common/debug.h"
#include "fft/kiss_fftr.h"
#include "fft/kiss_fft.h"

#include "fft/esp-dsp.h"

#ifndef LOG_LEVEL
#define LOG_LEVEL LOG_DEBUG
#endif

const int kMin_score = FT8_MIN_SCORE; // Minimum sync score threshold for candidates
const int kMax_candidates = FT8_MAX_CANDIDATES; // for 12 kHz sample rate; scaled for other sample rates
const int kLDPC_iterations = FT8_LDPC_ITERATIONS;

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

static double elapsed_ms(const struct timespec *t0, const struct timespec *t1)
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
    LOG(LOG_DEBUG, "Waterfall size = %zu\n", mag_size);
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

#include "esp_heap_caps.h"

void monitor_init(monitor_t* me, const monitor_config_t* cfg)
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
        "Monitor: sample_rate=%d time_osr=%d freq_osr=%d\n",
        cfg->sample_rate,
        cfg->time_osr,
        cfg->freq_osr);

    LOG(LOG_INFO,
        "Block size = %d samples (%.3f ms)\n",
        me->block_size,
        1000.0f * symbol_period);

    LOG(LOG_INFO,
        "Subblock size = %d samples (%.3f ms)\n",
        me->subblock_size,
        1000.0f * symbol_period / cfg->time_osr);

    LOG(LOG_INFO,
        "N_FFT = %d samples (%.3f ms)\n",
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
        LOG(LOG_ERROR, "monitor_init: FFT buffer allocation failed\n");
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

    for (int i = 0; i < me->nfft; ++i)
        me->window[i] = hann_i(i, me->nfft);

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
        "FFT work area = %zu bytes\n",
        fft_work_size);

    me->fft_work =
        heap_caps_malloc(
            fft_work_size,
            MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT);

    if (!me->fft_work)
    {
        LOG(LOG_ERROR,
            "monitor_init: FFT work allocation failed\n");
        return;
    }

    me->fft_cfg =
        kiss_fftr_alloc(
            me->nfft,
            0,
            me->fft_work,
            &fft_work_size);

#if defined ARDUINO_ARCH_ESP32
    esp_dsp_fftr_init();
#endif

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

void monitor_init_ori(monitor_t* me, const monitor_config_t* cfg)
{
    float slot_time = (cfg->protocol == PROTO_FT4) ? FT4_SLOT_TIME : FT8_SLOT_TIME;
    float symbol_period = (cfg->protocol == PROTO_FT4) ? FT4_SYMBOL_PERIOD : FT8_SYMBOL_PERIOD;
    // Compute DSP parameters that depend on the sample rate
    me->block_size = (int)(cfg->sample_rate * symbol_period); // samples corresponding to one FSK symbol
    me->subblock_size = me->block_size / cfg->time_osr;
    me->nfft = me->block_size * cfg->freq_osr;
    me->fft_norm = 2.0f / me->nfft;
    // const int len_window = 1.8f * me->block_size; // hand-picked and optimized

    me->window = (float *)malloc(me->nfft * sizeof(me->window[0]));
    for (int i = 0; i < me->nfft; ++i)
    {
        // window[i] = 1;
        me->window[i] = hann_i(i, me->nfft);
        // me->window[i] = blackman_i(i, me->nfft);
        // me->window[i] = hamming_i(i, me->nfft);
        // me->window[i] = (i < len_window) ? hann_i(i, len_window) : 0;
    }
    me->last_frame = (float *)malloc(me->nfft * sizeof(me->last_frame[0]));
    me->timedata = (kiss_fft_scalar*)malloc(me->nfft * sizeof(me->timedata[0]));
    me->freqdata = (kiss_fft_cpx*)malloc((me->nfft / 2 + 1) * sizeof(me->freqdata[0]));

    size_t fft_work_size;
    kiss_fftr_alloc(me->nfft, 0, 0, &fft_work_size);

    LOG(LOG_INFO, "Block size = %d\n", me->block_size);
    LOG(LOG_INFO, "Subblock size = %d\n", me->subblock_size);
    LOG(LOG_INFO, "N_FFT = %d\n", me->nfft);
    LOG(LOG_DEBUG, "FFT work area = %zu\n", fft_work_size);

    me->fft_work = malloc(fft_work_size);
    LOG(LOG_DEBUG, "init FFT for %d points work area allocated at %p\n", me->nfft, me->fft_work);
    me->fft_cfg = kiss_fftr_alloc(me->nfft, 0, me->fft_work, &fft_work_size);

    #if defined ARDUINO_ARCH_ESP32
    dsps_fft2r_init_fc32(NULL, 2048); // Configura per 2048 punti
    #endif

    const int max_blocks = (int)(slot_time / symbol_period);
    const int num_bins = (int)(cfg->sample_rate * symbol_period / 2);
    waterfall_init(&me->wf, max_blocks, num_bins, cfg->time_osr, cfg->freq_osr);
    me->wf.protocol = cfg->protocol;
    me->symbol_period = symbol_period;

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
    struct timespec t0, t1;
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
    const float fft_norm   = me->fft_norm;

    /*
     * One monitor_process() call receives exactly one symbol
     * worth of new samples.
     */
    int frame_pos = 0;

    int offset =
        me->wf.num_blocks * me->wf.block_stride;

    float local_max_mag = me->max_mag;

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
         * [---------- 3072 ----------][--1024--]
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
         */
        for (int i = 0; i < nfft; ++i)
        {
            timedata[i] =
                fft_norm *
                window[i] *
                last_frame[i];
        }

        /*
         * FFT.
         */
#if defined ARDUINO_ARCH_ESP32

        esp_dsp_fftr(
            me->fft_cfg,
            timedata,
            (float *)me->freqdata,
            me->nfft
        );
        // kiss_fftr(
        //     me->fft_cfg,
        //     timedata,
        //     me->freqdata);


#else

        kiss_fftr(
            me->fft_cfg,
            timedata,
            me->freqdata);

#endif

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
        for (int freq_sub = 0;
             freq_sub < freq_osr;
             ++freq_sub)
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
    }

    me->max_mag = local_max_mag;

    /*
     * One complete 160-ms block has now been added.
     */
    ++me->wf.num_blocks;

    clock_gettime(CLOCK_MONOTONIC, &t1);

    LOG(LOG_DEBUG,
        "[ft8] monitor_process: %.3f ms\n",
        elapsed_ms(&t0, &t1));
}


void monitor_reset(monitor_t* me)
{
    me->wf.num_blocks = 0;
    me->max_mag = 0;
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


// ------------------------------------------------------------------------------------
// Process a buffer already loaded from a file
// Pass precise time of signal[0] (including fractional second) so we can reference to it
// ------------------------------------------------------------------------------------
int process_buffer(float const *signal,int sample_rate, int num_samples, bool is_ft8, float base_freq, struct tm const *tmp, double sec){
  assert(signal != NULL && tmp != NULL);

  struct timespec t_wf0 = {0};
  struct timespec t_wf1 = {0};
  struct timespec t_dec0 = {0};
  struct timespec t_dec1 = {0};

  LOG(LOG_INFO, "Sample rate %d Hz, %d samples, %.3f seconds\n", sample_rate, num_samples, (double)num_samples / sample_rate);

  clock_gettime(CLOCK_MONOTONIC, &t_wf0);

  // Compute Waterfall accumulation (FFT)
  monitor_t mon = {0};
  monitor_config_t const mon_cfg = {
    .f_min = 100,
    .f_max = sample_rate/2 - 500, // allow room for the receiver filter rolloff
    .sample_rate = sample_rate,
    .time_osr = kTime_osr,
    .freq_osr = kFreq_osr,
    .protocol = is_ft8 ? PROTO_FT8 : PROTO_FT4
  };

  monitor_init(&mon, &mon_cfg);
  LOG(LOG_DEBUG, "Waterfall allocated %d blocks of size %d\n", mon.wf.max_blocks, mon.block_size);
  for (int frame_pos = 0; frame_pos + mon.block_size <= num_samples; frame_pos += mon.block_size){
      // Process the waveform data frame by frame - you could have a live loop here with data from an audio device
      // (cool, now that we can get sample timings - KA9Q)
      monitor_process(&mon, signal + frame_pos);
  }

 clock_gettime(CLOCK_MONOTONIC, &t_wf1);
  LOG(LOG_INFO, "Waterfall accumulation: %d blocks in %.3f ms\n", mon.wf.num_blocks, elapsed_ms(&t_wf0, &t_wf1));
  LOG(LOG_INFO, "Max magnitude: %.1f dB\n", mon.max_mag);
  
  float const noise_power = estimate_global_noise_power(&mon.wf);
  
  // Find top candidates by Costas sync score and localize them in time and frequency
  int const candidate_size = (mon_cfg.f_max * kMax_candidates) / 3000; // Scale by bandwidth relative to the original 3 kHz
  candidate_t candidate_list[candidate_size];
  int num_candidates = ft8_find_sync(&mon.wf, candidate_size, candidate_list, kMin_score);

  // Hash table for decoded messages (to check for duplicates)
  int num_decoded = 0;
  // Pointer to kMax_decoded_messages-element array of message_t structures
  message_t *decoded = calloc(sizeof(message_t), kMax_decoded_messages);
  // Pointer to kMax_decoded_messsages-element array of pointers to message_t structures
  message_t **decoded_hashtable = calloc(sizeof(message_t *), kMax_decoded_messages);

  clock_gettime(CLOCK_MONOTONIC, &t_dec0);
  LOG(LOG_INFO, "Found %d candidates with score from %d in %.3f milliseconds\n", num_candidates, kMin_score, elapsed_ms(&t_wf1, &t_dec0));

  // Go over candidates and attempt to decode messages
  for (int idx = 0; idx < num_candidates; ++idx){
      const candidate_t* cand = &candidate_list[idx];
      if (cand->score < kMin_score)
	      continue;

      float const freq_hz = (cand->freq_offset + (float)cand->freq_sub / mon.wf.freq_osr) / mon.symbol_period;
      float const time_sec = (cand->time_offset + (float)cand->time_sub / mon.wf.time_osr) * mon.symbol_period;

      message_t message = {0}; // Written by ft8_decode()
      decode_status_t status = {0}; // ditto
      uint8_t plain174[FTX_LDPC_N];
      if (!ft8_decode(&mon.wf, cand, &message, kLDPC_iterations, &status, plain174)){
	      // printf("000000 %3d %+4.2f %4.0f ~  ---\n", cand->score, time_sec, freq_hz);
        if (status.ldpc_errors > 0){
          LOG(LOG_DEBUG, "LDPC decode: %d errors\n", status.ldpc_errors);
        }else if (status.crc_calculated != status.crc_extracted){
          LOG(LOG_DEBUG, "CRC mismatch!\n");
        }else if (status.unpack_status != 0){
          LOG(LOG_DEBUG, "Error while unpacking!\n");
        }
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

      LOG(LOG_DEBUG, "Checking hash table for %4.1fs / %4.1fHz [%d]...\n", time_sec, freq_hz, cand->score);
      int idx_hash = message.hash % kMax_decoded_messages;
      bool found_empty_slot = false;
      bool found_duplicate = false;
      do{
        if (decoded_hashtable[idx_hash] == NULL){
            LOG(LOG_DEBUG, "Found an empty slot\n");
            found_empty_slot = true;
        }else if ((decoded_hashtable[idx_hash]->hash == message.hash) && (0 == strcmp(decoded_hashtable[idx_hash]->text, message.text))){
            LOG(LOG_DEBUG, "Found a duplicate [%s]\n", message.text);
            found_duplicate = true;
        }else{
            LOG(LOG_DEBUG, "Hash table clash!\n");
            // Move on to check the next entry in hash table
            idx_hash = (idx_hash + 1) % kMax_decoded_messages;
        }
      } while (!found_empty_slot && !found_duplicate);

      if (found_empty_slot){
        // Fill the empty hashtable slot
        decoded[idx_hash] = message;
        decoded_hashtable[idx_hash] = &decoded[idx_hash];
  
        fprintf(stdout,"%4d/%02d/%02d %02d:%02d:%02d %3d %+4.3lf %'.1lf ~ %s\n",
            tmp->tm_year + 1900,
            tmp->tm_mon + 1,
            tmp->tm_mday,
            tmp->tm_hour,
            tmp->tm_min,
            tmp->tm_sec,
            (int)lroundf(message.snr_db),
            message.time_sec,
            message.freq_hz,
            message.text);
        
            ++num_decoded;
      }
  }
  
  clock_gettime(CLOCK_MONOTONIC, &t_dec1);
  LOG(LOG_INFO, "On %d candidates, decoded %d messages in %.3f ms\n", num_candidates, num_decoded, elapsed_ms(&t_dec0, &t_dec1));
  LOG(LOG_INFO, "Decoded %d messages\n", num_decoded);
  
  // Decoded messages are spread throughout hash table, so sort the whole thing including null entries
  qsort(decoded_hashtable, kMax_decoded_messages, sizeof *decoded_hashtable, mcompare);
  // Empty entries sorted to top, so first num_decoded elements of decoded_hashtable are valid
  double tbase = tmp->tm_sec; // Full seconds and fraction in minute, should be just above (not below) period multiple
  tbase = is_ft8 ? fmod(tbase,15.0) : fmod(tbase,7.5); // seconds after start of cycle (0/15/30/45 or 0/7.5/15/etc)
  tbase += sec; // sec could be negative, so add it only now

  for(int i=0; i < num_decoded; i++){
    message_t const *mp = decoded_hashtable[i];
    if(mp == NULL)
      continue; // Shouldn't happen

    // fprintf(stdout,"%4d/%02d/%02d %02d:%02d:%02d %3d %+4.3lf %'.1lf ~ %s\n",
	  //   tmp->tm_year + 1900,
	  //   tmp->tm_mon + 1,
	  //   tmp->tm_mday,
	  //   tmp->tm_hour,
	  //   tmp->tm_min,
	  //   tmp->tm_sec,
    //   (int)lroundf(mp->snr_db),
	  //   tbase + mp->time_sec,
	  //   mp->freq_hz,
    //   mp->text);
  }
  free(decoded);
  free(decoded_hashtable);

  monitor_free(&mon);
  return 0; // Caller frees signal
}
