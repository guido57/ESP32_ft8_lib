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
#include "ddc_ft8.h"

#include "common/wave.h"
#include "common/debug.h"
#include "fft/kiss_fftr.h"
#include "fft/kiss_fft.h"
#if defined(ARDUINO_ARCH_ESP32)
#include <esp_heap_caps.h>
#include "fft/esp-dsp.h"
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
#if defined ARDUINO_ARCH_ESP32
    bool esp_dsp_fft_ready; ///< ESP-DSP tables were allocated successfully
#endif
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
        "Monitor: sample_rate=%d num_samples=%d time_osr=%d freq_osr=%d\n",
        cfg->sample_rate,
        num_samples,
        cfg->time_osr,
        cfg->freq_osr);

    LOG(LOG_INFO,
        "Block size = %d samples (%.3f ms) Subblock size = %d samples (%.3f ms) N_FFT = %d samples (%.3f ms)\n",
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
    me->esp_dsp_fft_ready = esp_dsp_fftr_init();
    if (me->esp_dsp_fft_ready)
        LOG(LOG_INFO, "monitor_init: using ESP-DSP FFT\n");
    else
        LOG(LOG_WARN, "monitor_init: ESP-DSP FFT unavailable; using KISS FFT\n");
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
         */
        for (int i = 0; i < nfft; ++i)
        {
            timedata[i] = window[i] * last_frame[i];
        }

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

/*
 * ft8_decode() only needs 79 time rows of eight tone magnitudes.  Build that
 * minimal waterfall directly from an already-subtracted 250-Hz IQ stream so
 * a residual signal is not re-contaminated by the original wideband
 * waterfall.  The candidate is fine-synchronized first, then every local
 * waterfall bin is one exact 32-sample (160-ms) complex DFT.
 */
static bool decode_ddc_residual_candidate_200(
    const std::vector<IQ>& iq,
    int iq_first_input,
    const CostasCandidate& candidate,
    float ddc_lo_hz,
    message_t* message,
    decode_status_t* status,
    uint8_t* plain174,
    float* refined_delay_s,
    float* refined_freq_hz)
{
    constexpr int ddc_sample_rate = 250;
    constexpr int ddc_symbol_samples = 40;
    constexpr int ft8_symbols = 79;
    constexpr int ft8_tones = 8;
    constexpr float tone_spacing_hz = 6.25f;

    const float iq_start_s =
        (float)iq_first_input / 12000.0f;
    const CostasCandidate refined = refine_costas_candidate_200(
        iq, iq_first_input, ddc_lo_hz, candidate);
    const int start = (int)lroundf(
        (refined.delay_s - iq_start_s) * ddc_sample_rate);
    if (start < 0 || start >= (int)iq.size())
        return false;

    const int available_symbols = std::min(
        ft8_symbols,
        ((int)iq.size() - start) / ddc_symbol_samples);
    // All three Costas blocks must be present; a clipped data tail is okay.
    if (available_symbols < 75)
        return false;

    const float fine_freq = refined.freq_hz;
    const float fine_delay = refined.delay_s;

    uint8_t magnitudes[ft8_symbols * ft8_tones] = {};
    for (int symbol = 0; symbol < available_symbols; ++symbol)
    {
        for (int tone = 0; tone < ft8_tones; ++tone)
        {
            const float freq = (fine_freq - ddc_lo_hz) +
                tone_spacing_hz * (float)tone;
            const float phase_step =
                -2.0f * (float)M_PI * freq /
                (float)ddc_sample_rate;
            const float step_c = cosf(phase_step);
            const float step_s = sinf(phase_step);
            float c = 1.0f;
            float s = 0.0f;
            float re = 0.0f;
            float im = 0.0f;

            const int first = start + symbol * ddc_symbol_samples;
            for (int sample = 0; sample < ddc_symbol_samples; ++sample)
            {
                const IQ z = iq[first + sample];
                re += z.i * c - z.q * s;
                im += z.i * s + z.q * c;

                const float next_c = c * step_c - s * step_s;
                s = s * step_c + c * step_s;
                c = next_c;
            }

            const float power = re * re + im * im;
            int scaled = (int)(20.0f * log10f(1e-12f + power) + 240.0f);
            if (scaled < 0)
                scaled = 0;
            else if (scaled > 255)
                scaled = 255;
            magnitudes[symbol * ft8_tones + tone] = (uint8_t)scaled;
        }
    }

    waterfall_t residual_wf = {
        .max_blocks = ft8_symbols,
        .num_blocks = available_symbols,
        .num_bins = ft8_tones,
        .time_osr = 1,
        .freq_osr = 1,
        .mag = magnitudes,
        .block_stride = ft8_tones,
        .protocol = PROTO_FT8};
    candidate_t residual_cand = {
        .score = 0,
        .time_offset = 0,
        .freq_offset = 0,
        .time_sub = 0,
        .freq_sub = 0};

    if (refined_delay_s)
        *refined_delay_s = fine_delay;
    if (refined_freq_hz)
        *refined_freq_hz = fine_freq;
    return ft8_decode(&residual_wf,
                      &residual_cand,
                      message,
                      kLDPC_iterations,
                      status,
                      plain174);
}

// Estimate residual-message SINR without changing the active SIC stream.
// Measure the coherently integrated power in each decoded FT8 tone and use
// unoccupied 6.25-Hz bins as a robust, local noise estimate.  This avoids the
// old model/residual RMS ratio, which counted every other signal in the DDC
// passband as noise and made the result depend strongly on cancellation order.
static bool estimate_ddc_snr_2500(
    const std::vector<IQ>& iq,
    int iq_first_input,
    const uint8_t plain174[FTX_LDPC_N],
    float delay_s,
    float freq_hz,
    float ddc_lo_hz,
    float* snr_2500_db)
{
    constexpr int ddc_sample_rate = 250;
    constexpr int ddc_symbol_samples = 40;
    constexpr int ft8_symbols = 79;
    constexpr float tone_spacing_hz = 6.25f;
    constexpr float reference_bandwidth_hz = 2500.0f;

    uint8_t tones[79];
    get_ft8_tones_from_plain174(plain174, tones);

    const float iq_start_s = (float)iq_first_input / 12000.0f;
    const int start = (int)lroundf(
        (delay_s - iq_start_s) * (float)ddc_sample_rate);
    if (start < 0 || start >= (int)iq.size())
        return false;

    const int available_symbols = std::min(
        ft8_symbols,
        ((int)iq.size() - start) / ddc_symbol_samples);
    if (available_symbols < 75)
        return false;

    // Complex samples repeat every Fs in frequency.  Folding here affects
    // only this measurement and keeps aliased DDC coordinates equivalent.
    const float relative_tone0_hz =
        remainderf(freq_hz - ddc_lo_hz, (float)ddc_sample_rate);

    auto bin_power = [&](int first, float bin_hz) {
        const float phase_step =
            -2.0f * (float)M_PI * bin_hz / (float)ddc_sample_rate;
        const float step_c = cosf(phase_step);
        const float step_s = sinf(phase_step);
        float c = 1.0f;
        float s = 0.0f;
        float re = 0.0f;
        float im = 0.0f;
        for (int sample = 0; sample < ddc_symbol_samples; ++sample)
        {
            const IQ z = iq[first + sample];
            re += z.i * c - z.q * s;
            im += z.i * s + z.q * c;
            const float next_c = c * step_c - s * step_s;
            s = s * step_c + c * step_s;
            c = next_c;
        }
        return re * re + im * im;
    };

    std::vector<float> signal_powers;
    std::vector<float> noise_powers;
    signal_powers.reserve(available_symbols);
    noise_powers.reserve(available_symbols * 16);

    for (int symbol = 0; symbol < available_symbols; ++symbol)
    {
        const int first = start + symbol * ddc_symbol_samples;
        const float signal_hz = relative_tone0_hz +
            tone_spacing_hz * (float)tones[symbol];
        signal_powers.push_back(bin_power(first, signal_hz));

        // Bins 0..7 can contain this FT8 transmission.  Sample bins on both
        // sides, staying away from the transition region of the 80-Hz LPF.
        for (int bin = -12; bin <= 19; ++bin)
        {
            if (bin >= 0 && bin <= 7)
                continue;
            const float noise_hz = relative_tone0_hz +
                tone_spacing_hz * (float)bin;
            if (fabsf(noise_hz) > 70.0f)
                continue;
            noise_powers.push_back(bin_power(first, noise_hz));
        }
    }

    if (signal_powers.empty() || noise_powers.empty())
        return false;

    auto median = [](std::vector<float>& values) {
        const size_t middle = values.size() / 2;
        std::nth_element(values.begin(), values.begin() + middle, values.end());
        return values[middle];
    };
    const float measured_tone_power = median(signal_powers);
    // For complex Gaussian noise, FFT-bin power is exponential and its
    // median is ln(2) times its mean.  Correct the robust median back to the
    // mean power required by the SNR definition.
    const float noise_bin_power = median(noise_powers) / logf(2.0f);
    const float net_signal_power = measured_tone_power - noise_bin_power;
    if (noise_bin_power <= 0.0f || net_signal_power <= 0.0f)
        return false;

    const float snr_bin_db =
        10.0f * log10f(net_signal_power / noise_bin_power);
    *snr_2500_db = snr_bin_db -
        10.0f * log10f(reference_bandwidth_hz / tone_spacing_hz);
    return std::isfinite(*snr_2500_db);
}

/* The Costas peak is occasionally displaced by residual energy.  Do a small
 * decode-directed search for the strongest hypothesis only; this is much
 * cheaper than applying an LDPC grid to every blind candidate. */
static bool retry_ddc_residual_candidate_200(
    const std::vector<IQ>& iq,
    int iq_first_input,
    const CostasCandidate& candidate,
    float ddc_lo_hz,
    message_t* message,
    decode_status_t* status,
    uint8_t* plain174,
    float* refined_delay_s,
    float* refined_freq_hz)
{
    static constexpr int timing_offsets[] = {-4, 4, -8, 8};
    static constexpr float frequency_offsets[] = {0.0f, -3.0f, 3.0f};

    decode_status_t best_status = *status;
    for (int timing_offset : timing_offsets)
    {
        for (float frequency_offset : frequency_offsets)
        {
            CostasCandidate trial = candidate;
            trial.delay_s += (float)timing_offset / 250.0f;
            trial.freq_hz += frequency_offset;

            message_t trial_message = {0};
            decode_status_t trial_status = {0};
            uint8_t trial_plain174[FTX_LDPC_N];
            float trial_delay_s = 0.0f;
            float trial_freq_hz = 0.0f;
            if (decode_ddc_residual_candidate_200(
                    iq, iq_first_input, trial, ddc_lo_hz,
                    &trial_message, &trial_status, trial_plain174,
                    &trial_delay_s, &trial_freq_hz))
            {
                *message = trial_message;
                *status = trial_status;
                memcpy(plain174, trial_plain174, FTX_LDPC_N);
                *refined_delay_s = trial_delay_s;
                *refined_freq_hz = trial_freq_hz;
                return true;
            }
            if (trial_status.ldpc_errors < best_status.ldpc_errors)
                best_status = trial_status;
        }
    }
    *status = best_status;
    return false;
}

// ------------------------------------------------------------------------------------
// Process "num_samples" samples stored in "signal"
// ------------------------------------------------------------------------------------
int process_buffer(float *samples,int sample_rate, int num_samples, 
    bool is_ft8, int cand_to_subtract, float freq_hz_subtract, float time_delay_subtract, float base_freq, struct tm const *tmp, double sec){
  assert(samples != NULL && tmp != NULL);

  // This branch never mutates samples, so every stage can use this immutable
  // view of the original mixture without an ESP32-S3-prohibitive RAM copy.
  const float *const raw_samples = samples;


  // Compute Waterfall accumulation (FFT)
  monitor_config_t const mon_cfg = {
        .f_min = 100,
        .f_max = (sample_rate)/2.0f - 500.0f, // allow room for the receiver filter rolloff
        .sample_rate = sample_rate,
        .time_osr = FT8_TIME_OSR,
        .freq_osr = FT8_FREQ_OSR,
        .protocol = is_ft8 ? PROTO_FT8 : PROTO_FT4
  };

  monitor_t mon = {0};
  
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
    LOG(LOG_WARN, "Decoded-message table is full; dropping [%s]\\n",
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

  // A 250-Hz complex stream aliases signals separated by 250 Hz.  DDC bands
  // overlap, so the same already-known message can otherwise be sent through
  // the costly local DFT and LDPC decoder repeatedly under an alias.  The
  // blind candidates are only 16 ms / 1 Hz coarse estimates, hence the
  // deliberately wider-but-still-conservative matching tolerances below.
  auto is_known_ddc_candidate = [&](const CostasCandidate& candidate,
                                    float ddc_lo_hz) {
    const float candidate_relative = remainderf(
        candidate.freq_hz - ddc_lo_hz, 250.0f);
    for (int i = 0; i < kMax_decoded_messages; ++i) {
      const message_t* known = decoded_hashtable[i];
      if (!known)
        continue;
      const float known_relative = remainderf(
          known->freq_hz - ddc_lo_hz, 250.0f);
      if (fabsf(candidate_relative - known_relative) > 2.0f)
        continue;
      // Waterfall candidates are timestamped one symbol after the actual
      // FT8 start; DDC residual candidates use the start directly.
      const float waterfall_start = known->time_sec - FT8_SYMBOL_PERIOD;
      if (fabsf(candidate.delay_s - known->time_sec) < 0.060f ||
          fabsf(candidate.delay_s - waterfall_start) < 0.060f)
        return true;
    }
    return false;
  };
  
  // ---------------------------------------------------
  //  ONE PASS: standard waterfall decode followed by local DDC-SIC.
  //  Do not re-run a waterfall on a conventionally subtracted waveform.
  // ---------------------------------------------------
  for(int pass=0; pass<1; ++pass){  
    printf("========================================\n");
    printf("Starting pass %d\n", pass);
    printf("========================================\n");
    
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

    // Compute Waterfall accumulation (FFT)
    monitor_config_t const mon_cfg = {
            .f_min = 100,
            .f_max = (sample_rate)/2.0f - 500.0f, // allow room for the receiver filter rolloff
            .sample_rate = sample_rate,
            .time_osr = FT8_TIME_OSR,
            .freq_osr = FT8_FREQ_OSR,
            .protocol = is_ft8 ? PROTO_FT8 : PROTO_FT4
    };

    monitor_init(&mon, &mon_cfg, num_samples);
    LOG(LOG_DEBUG, "Waterfall allocated %d blocks of size %d\n", mon.wf.max_blocks, mon.block_size);

    for (int frame_pos = 0; frame_pos + mon.block_size <= num_samples; frame_pos += mon.block_size){
        // Process the waveform data frame by frame - you could have a live loop here with data from an audio device
        // (cool, now that we can get sample timings - KA9Q)
        monitor_process(&mon, raw_samples + frame_pos);
    }

    clock_gettime(CLOCK_MONOTONIC, &t_wf1);
    LOG(LOG_INFO, "Waterfall accumulation: %d blocks in %.3f ms  max magnitude %.1f dB\n", mon.wf.num_blocks, elapsed_ms(&t_wf0, &t_wf1), mon.max_mag);
    
    float const noise_power = estimate_global_noise_power(&mon.wf);
    
    // Find top candidates by Costas sync score and localize them in time and frequency
    int const candidate_size = (mon_cfg.f_max * kMax_candidates) / 3000; // Scale by bandwidth relative to the original 3 kHz
    candidate_t candidate_list[candidate_size];
    int num_candidates = ft8_find_sync(&mon.wf, candidate_size, candidate_list, kMin_score);

    clock_gettime(CLOCK_MONOTONIC, &t_dec0);
    LOG(LOG_INFO, "Found %d candidates with score from %d in %.3f milliseconds\n", num_candidates, kMin_score, elapsed_ms(&t_wf1, &t_dec0));

    // ==================================================
    // Go over candidates and attempt to decode messages
    // ==================================================
    for (int idx = 0; idx < num_candidates; ++idx){
    
        clock_gettime(CLOCK_MONOTONIC, &t_cand0);
        const candidate_t* cand = &candidate_list[idx];
        if (cand->score < kMin_score)
            continue;

        float const freq_hz = (cand->freq_offset + (float)cand->freq_sub / mon.wf.freq_osr) / mon.symbol_period;
        float const time_sec = (cand->time_offset + (float)cand->time_sub / mon.wf.time_osr) * mon.symbol_period;

        message_t message = {0};      // Written by ft8_decode()
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
            clock_gettime(CLOCK_MONOTONIC, &t_ndec);
            LOG(LOG_DEBUG, "Decoding failed in %.3f milliseconds\n", elapsed_ms(&t_cand0, &t_ndec));
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

        LOG(LOG_DEBUG, "Checking decoded-message table for %4.1fs / %4.1fHz [%d]...\n", time_sec, freq_hz, cand->score);
        if (store_decoded_message(message)){
            clock_gettime(CLOCK_MONOTONIC, &t_dec1);
            LOG(LOG_INFO,"Decoded cand=%d score=%3d snr=%+.1f %+4.3f %4.3f ~  %s took %.3f msec\r\n", 
                idx, cand->score, message.snr_db, time_sec-0.14f, freq_hz,
                message.text, elapsed_ms(&t_cand0, &t_dec1));
        
            if(pass == 0){    

                uint8_t tones[79];
                get_ft8_tones_from_plain174(plain174, tones);

                // DDC has its own Costas-based time/frequency refinement
                // below.  Do not run the legacy wideband refiners or raw
                // waveform RMS diagnostics here: their results were only
                // printed and never influenced DDC cancellation or decode.

                std::vector<IQ> ddc_output;

                struct timespec ddc_t0;
                struct timespec ddc_t1;

                clock_gettime(CLOCK_MONOTONIC, &ddc_t0);

                // Centre the complex DDC at this candidate's coarse frequency.
                ddc_process_fast(
                    freq_hz, raw_samples, num_samples, ddc_output);


                clock_gettime(
                    CLOCK_MONOTONIC,
                    &ddc_t1);


                const double ddc_time_ms =
                    elapsed_ms(
                        &ddc_t0,
                        &ddc_t1);

                LOG(LOG_DEBUG, "DDC processing took %.3f msec\r\n", ddc_time_ms);

                // ------------------------------------------------------------
                // DDC FINE FREQUENCY AND DELAY ESTIMATE
                // ------------------------------------------------------------
                float fine_freq_offset_hz = 0.0f;
                float fine_time_offset_samp = 0.0f;
                float coarse_delay_time_ms = 0.0f;
                float fine_frequency_time_ms = 0.0f;
                float fine_delay_time_ms = 0.0f;

                struct timespec fine_sync_t0;
                struct timespec fine_sync_t1;

                clock_gettime(
                    CLOCK_MONOTONIC,
                    &fine_sync_t0);

                /*
                 * With the DDC LO at freq_hz, the base FT8 tone is at DC.
                 * The Costas offsets are therefore the DFT bins directly.
                 */
                const int base_bin = 0;

                constexpr float ddc_input_sample_rate = 12000.0f;
                constexpr float ddc_output_sample_rate = 250.0f;
                constexpr float ddc_decimation =
                    ddc_input_sample_rate / ddc_output_sample_rate;
                constexpr float ddc_fir_group_delay_samples = 60.0f;

                /*
                 * Candidate time points one FT8 symbol after the waveform
                 * start.  Map the raw-input start to the causal FIR's output
                 * coordinate, whose samples are delayed by 60 input samples.
                 */
                const float ddc_coarse_start_s =
                    time_sec - mon.symbol_period;
                const int ddc_start_output =
                    static_cast<int>(std::lround(
                        (ddc_coarse_start_s * ddc_input_sample_rate +
                         ddc_fir_group_delay_samples) /
                        ddc_decimation));
                                
                ddc_estimate_fine_sync(
                    ddc_output,
                    ddc_start_output,
                    base_bin,
                    COSTAS,
                    COSTAS_LEN,
                    COSTAS_START,
                    3,
                    &fine_freq_offset_hz,
                    &fine_time_offset_samp,
                    &coarse_delay_time_ms,
                    &fine_frequency_time_ms,
                    &fine_delay_time_ms);

                clock_gettime(
                    CLOCK_MONOTONIC,
                    &fine_sync_t1);

                const float ddc_fine_freq_hz =
                    freq_hz + fine_freq_offset_hz;
                const float ddc_fine_delay_s =
                    ((static_cast<float>(ddc_start_output) +
                      fine_time_offset_samp) *
                     ddc_decimation - ddc_fir_group_delay_samples) /
                    ddc_input_sample_rate;

                LOG(LOG_INFO,
                    "DDC refine: freq=%.6f Hz (offset=%+.6f), "
                    "delay=%.6f s (offset=%+.3f samples)\r\n",
                    ddc_fine_freq_hz,
                    fine_freq_offset_hz,
                    ddc_fine_delay_s,
                    fine_time_offset_samp);
                LOG(LOG_DEBUG,
                    "DDC estimator timing: search=%.3f ms, frequency=%.3f ms, delay=%.3f ms\r\n",
                    coarse_delay_time_ms,
                    fine_frequency_time_ms,
                    fine_delay_time_ms);

                // ------------------------------------------------------------
                // DDC-IQ subtraction and blind Costas search
                // ------------------------------------------------------------
                // ddc_output[0] is the causal FIR output at raw sample zero,
                // which represents input time -60 samples after group-delay
                // correction.  Use that same coordinate for subtraction and
                // for the blind residual search.
                constexpr int ddc_iq_first_input = -60;
                float ddc_rms_before = 0.0f;
                float ddc_rms_after = 0.0f;
                float ddc_rms_model = 0.0f;
                float ddc_residual_projection_db = 0.0f;
                struct timespec ddc_subtract_t0;
                struct timespec ddc_subtract_t1;
                clock_gettime(CLOCK_MONOTONIC, &ddc_subtract_t0);

                const bool ddc_subtracted = subtract_ft8_message_200(
                    ddc_output,
                    ddc_iq_first_input,
                    tones,
                    ddc_fine_delay_s,
                    ddc_fine_freq_hz,
                    freq_hz,
                    &ddc_rms_before,
                    &ddc_rms_after,
                    &ddc_rms_model,
                    &ddc_residual_projection_db);

                clock_gettime(CLOCK_MONOTONIC, &ddc_subtract_t1);

                if (!ddc_subtracted)
                {
                    LOG(LOG_WARN,
                        "DDC IQ subtraction skipped: message is outside the IQ window\r\n");
                }
                else
                {

                    LOG(LOG_DEBUG,
                        "DDC IQ subtract: rms %.6f -> %.6f (model %.6f, residual projection %.2f dB) in %.3f ms\r\n",
                        ddc_rms_before,
                        ddc_rms_after,
                        ddc_rms_model,
                        ddc_residual_projection_db,
                        elapsed_ms(&ddc_subtract_t0, &ddc_subtract_t1));

                    // Search the complete nominal DDC passband around the
                    // coarse tone.  Edge candidates may be attenuated by the
                    // 80-Hz low-pass filter, but remain useful for testing.
                    struct timespec ddc_blind_t0;
                    struct timespec ddc_blind_t1;
                    clock_gettime(CLOCK_MONOTONIC, &ddc_blind_t0);
                    const std::vector<CostasCandidate> residual_candidates =
                        find_costas_candidates_200(
                            ddc_output,
                            ddc_iq_first_input,
                            freq_hz,
                            freq_hz - 80.0f,
                            freq_hz + 80.0f);
                    clock_gettime(CLOCK_MONOTONIC, &ddc_blind_t1);

                    LOG(LOG_INFO,
                        "DDC residual Costas scan: %u candidates in %.3f ms\r\n",
                        (unsigned)residual_candidates.size(),
                        elapsed_ms(&ddc_blind_t0, &ddc_blind_t1));
                    bool have_first_residual = false;
                    uint8_t first_residual_plain174[FTX_LDPC_N] = {};
                    float first_residual_delay_s = 0.0f;
                    float first_residual_freq_hz = 0.0f;
                    for (size_t candidate_index = 0;
                         candidate_index < residual_candidates.size();
                         ++candidate_index)
                    {
                        const CostasCandidate& residual =
                            residual_candidates[candidate_index];
                        if (is_known_ddc_candidate(residual, freq_hz))
                        {
                            LOG(LOG_DEBUG,
                                "  DDC residual %u skipped: known alias\r\n",
                                (unsigned)(candidate_index + 1));
                            continue;
                        }
                        LOG(LOG_DEBUG,
                            "  DDC residual %u: delay %.6f s, tone-0 %.3f Hz, score %.6g\r\n",
                            (unsigned)(candidate_index + 1),
                            residual.delay_s,
                            residual.freq_hz,
                            residual.score);

                        message_t residual_message = {0};
                        decode_status_t residual_status = {0};
                        uint8_t residual_plain174[FTX_LDPC_N];
                        float residual_fine_delay_s = 0.0f;
                        float residual_fine_freq_hz = 0.0f;
                        bool residual_decoded =
                            decode_ddc_residual_candidate_200(
                                ddc_output,
                                ddc_iq_first_input,
                                residual,
                                freq_hz,
                                &residual_message,
                                &residual_status,
                                residual_plain174,
                                &residual_fine_delay_s,
                                &residual_fine_freq_hz);
                        if (!residual_decoded && candidate_index == 0)
                        {
                            residual_decoded = retry_ddc_residual_candidate_200(
                                ddc_output, ddc_iq_first_input, residual,
                                freq_hz, &residual_message, &residual_status,
                                residual_plain174, &residual_fine_delay_s,
                                &residual_fine_freq_hz);
                        }
                        if (residual_decoded)
                        {
                            residual_message.freq_hz = residual_fine_freq_hz;
                            residual_message.time_sec = residual_fine_delay_s;
                            residual_message.score = residual.score;
                            message_t* previous =
                                find_decoded_message(residual_message);
                            float residual_snr_2500_db =
                                previous ? previous->snr_db : NAN;
                            if (!previous)
                            {
                                (void)estimate_ddc_snr_2500(
                                    ddc_output, ddc_iq_first_input,
                                    residual_plain174, residual_fine_delay_s,
                                    residual_fine_freq_hz, freq_hz,
                                    &residual_snr_2500_db);
                                residual_message.snr_db =
                                    residual_snr_2500_db;
                            }
                            const bool is_new_residual =
                                store_decoded_message(residual_message);
                            LOG(LOG_INFO,
                                "    DDC residual decoded%s: snr2500=%+.1f dB, delay %.6f s, tone-0 %.6f Hz ~  %s\r\n",
                                is_new_residual ? "" : " (duplicate)",
                                residual_snr_2500_db,
                                residual_fine_delay_s,
                                residual_fine_freq_hz,
                                residual_message.text);
                            // Preserve the strongest decodable residual for
                            // one cancellation/reacquisition iteration.
                            if (!have_first_residual)
                            {
                                memcpy(first_residual_plain174,
                                       residual_plain174,
                                       sizeof(first_residual_plain174));
                                first_residual_delay_s = residual_fine_delay_s;
                                first_residual_freq_hz = residual_fine_freq_hz;
                                have_first_residual = true;
                            }
                            // A stronger residual may be a second message or
                            // an alias.  Keep checking this short candidate
                            // list for weaker, distinct FT8 messages.
                        }
                        else
                        {
                            LOG(LOG_DEBUG,
                                "    DDC residual decode failed: delay %.6f s, tone-0 %.6f Hz, LDPC errors %d\r\n",
                                residual_fine_delay_s,
                                residual_fine_freq_hz,
                                residual_status.ldpc_errors);
                        }
                    }

                    // Successive-interference cancellation is essential for
                    // overlapping FT8 signals: e.g. after IK4LZH is removed,
                    // remove JA1FWS OK2BV before attempting JA1FWS HA7CH.
                    if (have_first_residual)
                    {
                        uint8_t first_residual_tones[79];
                        get_ft8_tones_from_plain174(
                            first_residual_plain174, first_residual_tones);
                        if (subtract_ft8_message_200(
                                ddc_output, ddc_iq_first_input,
                                first_residual_tones,
                                first_residual_delay_s,
                                first_residual_freq_hz, freq_hz,
                                nullptr, nullptr, nullptr, nullptr))
                        {
                            const std::vector<CostasCandidate> second_candidates =
                                find_costas_candidates_200(
                                    ddc_output, ddc_iq_first_input, freq_hz,
                                    freq_hz - 80.0f, freq_hz + 80.0f);
                            LOG(LOG_INFO,
                                "DDC residual scan after subtracting %.3f Hz / %.3f s: %u candidates\r\n",
                                first_residual_freq_hz,
                                first_residual_delay_s,
                                (unsigned)second_candidates.size());
                            for (size_t second_index = 0;
                                 second_index < second_candidates.size();
                                 ++second_index)
                            {
                                const CostasCandidate& second =
                                    second_candidates[second_index];
                                if (is_known_ddc_candidate(second, freq_hz))
                                {
                                    LOG(LOG_DEBUG,
                                        "  DDC residual SIC %u skipped: known alias\r\n",
                                        (unsigned)(second_index + 1));
                                    continue;
                                }
                                message_t second_message = {0};
                                decode_status_t second_status = {0};
                                uint8_t second_plain174[FTX_LDPC_N];
                                float second_delay_s = 0.0f;
                                float second_freq_hz = 0.0f;
                                if (decode_ddc_residual_candidate_200(
                                        ddc_output, ddc_iq_first_input,
                                        second, freq_hz, &second_message,
                                        &second_status, second_plain174,
                                        &second_delay_s, &second_freq_hz))
                                {
                                    second_message.freq_hz = second_freq_hz;
                                    second_message.time_sec = second_delay_s;
                                    second_message.score = second.score;
                                    message_t* previous =
                                        find_decoded_message(second_message);
                                    float second_snr_2500_db =
                                        previous ? previous->snr_db : NAN;
                                    if (!previous)
                                    {
                                        (void)estimate_ddc_snr_2500(
                                            ddc_output, ddc_iq_first_input,
                                            second_plain174, second_delay_s,
                                            second_freq_hz, freq_hz,
                                            &second_snr_2500_db);
                                        second_message.snr_db =
                                            second_snr_2500_db;
                                    }
                                    const bool is_new_residual =
                                        store_decoded_message(second_message);
                                    LOG(LOG_INFO,
                                        "    DDC residual after cancellation%s: snr2500=%+.1f dB, delay %.6f s, tone-0 %.6f Hz ~  %s\r\n",
                                        is_new_residual ? "" : " (duplicate)",
                                        second_snr_2500_db,
                                        second_delay_s, second_freq_hz,
                                        second_message.text);
                                }
                            }
                        }
                    }
                }
                LOG(LOG_DEBUG, "Fine sync processing took %.3f msec\r\n", elapsed_ms(&fine_sync_t0, &fine_sync_t1));

                // No wideband subtract() here.  Cancellation is confined to
                // ddc_output, so later candidates always start from the
                // immutable raw waveform above.
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &t_cand1);
        LOG(LOG_DEBUG, "End of processing cand %d. It took %.3f msec\n", idx, elapsed_ms(&t_cand0, &t_cand1));
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t_dec1);
    LOG(LOG_INFO, "Pass %d: decoded %d messages on %d candidates, in %.3f ms\n", pass, num_decoded, num_candidates, elapsed_ms(&t_dec0, &t_dec1));
    

    // Decoded messages are spread throughout hash table, so sort the whole thing including null entries
    // qsort(decoded_hashtable, kMax_decoded_messages, sizeof *decoded_hashtable, mcompare);
    // Empty entries sorted to top, so first num_decoded elements of decoded_hashtable are valid
    double tbase = tmp->tm_sec; // Full seconds and fraction in minute, should be just above (not below) period multiple
    tbase = is_ft8 ? fmod(tbase,15.0) : fmod(tbase,7.5); // seconds after start of cycle (0/15/30/45 or 0/7.5/15/etc)
    tbase += sec; // sec could be negative, so add it only now

    // for(int i=0; i < num_decoded; i++){
    //     message_t const *mp = decoded_hashtable[i];
    //     if(mp == NULL)
    //     continue; // Shouldn't happen
    // }
    monitor_free(&mon);
  } // End of pass

  free(decoded);
  free(decoded_hashtable);
  return 0; // Caller frees signal
}


int process_buffer_ori(float const *signal,int sample_rate, int num_samples, bool is_ft8, float time_delay, struct tm const *tmp, double sec){
  assert(signal != NULL && tmp != NULL);
    
  float * samples_ = (float *)signal;

  struct timespec t_wf0 = {0};
  struct timespec t_wf1 = {0};
  struct timespec t_dec0 = {0};
  struct timespec t_dec1 = {0};

  clock_gettime(CLOCK_MONOTONIC, &t_wf0);

  // Compute Waterfall accumulation (FFT)
  monitor_t mon = {0};
  monitor_config_t const mon_cfg = {
    .f_min = 100,
    .f_max = (sample_rate)/2.0f - 500.0f, // allow room for the receiver filter rolloff
    .sample_rate = sample_rate,
    .time_osr = kTime_osr,
    .freq_osr = kFreq_osr,
    .protocol = is_ft8 ? PROTO_FT8 : PROTO_FT4
  };

  monitor_init(&mon, &mon_cfg, num_samples);
  LOG(LOG_DEBUG, "Waterfall allocated %d blocks of size %d\n", mon.wf.max_blocks, mon.block_size);
  for (int frame_pos = 0; frame_pos + mon.block_size <= num_samples; frame_pos += mon.block_size){
      // Process the waveform data frame by frame - you could have a live loop here with data from an audio device
      // (cool, now that we can get sample timings - KA9Q)
      monitor_process(&mon, signal + frame_pos);
  }

 clock_gettime(CLOCK_MONOTONIC, &t_wf1);
  LOG(LOG_INFO, "Waterfall accumulation: %d blocks in %.3f ms  max magnitude %.1f dB\n", mon.wf.num_blocks, elapsed_ms(&t_wf0, &t_wf1), mon.max_mag);
  
  float const noise_power = estimate_global_noise_power(&mon.wf);
  
  // Find top candidates by Costas sync score and localize them in time and frequency
  int const candidate_size = (mon_cfg.f_max * kMax_candidates) / 3000; // Scale by bandwidth relative to the original 3 kHz
  candidate_t candidate_list[candidate_size];
  int num_candidates = ft8_find_sync(&mon.wf, candidate_size, candidate_list, kMin_score);

  // Hash table for decoded messages (to check for duplicates)
  int num_decoded = 0;
  // Pointer to kMax_decoded_messages-element array of message_t structures
  message_t *decoded = (message_t *) calloc(sizeof(message_t), kMax_decoded_messages);
  // Pointer to kMax_decoded_messsages-element array of pointers to message_t structures
  message_t **decoded_hashtable = (message_t **) calloc(sizeof(message_t *), kMax_decoded_messages);

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
    
            fprintf(stdout,"%4d/%02d/%02d %02d:%02d:%02d %3d %.3f %.6f %3d ~ %s\n",
                tmp->tm_year + 1900,
                tmp->tm_mon + 1,
                tmp->tm_mday,
                tmp->tm_hour,
                tmp->tm_min,
                tmp->tm_sec,
                idx,
                freq_hz,
                time_sec,
                (int)lroundf(message.snr_db),
                message.text);
            
            ++num_decoded;

            printf("DECODED: %3d %+4.2f %4.0f ~  %s\n", cand->score, time_sec, freq_hz, message.text);
            uint8_t tones[79];
            get_ft8_tones_from_plain174(plain174, tones);

            float rms_before = rms(samples_, num_samples);
    
            printf("subtracting tones from waterfall... message.freq_hz=%.3f message.time_sec=%.6f num_samples=%d\n", 
                    freq_hz, time_sec, num_samples);
            
            subtract(tones,
                freq_hz,
                freq_hz,
                time_delay,
                samples_,
                num_samples,
                sample_rate);
            
                printf("SUBTRACT: time=%.6f sec freq=%.3f Hz samples=%d sample_rate=%d\n",
                time_delay,
                freq_hz,
                num_samples,
                sample_rate);
                
            float rms_after = rms(samples_, num_samples);
            double cancellation_db =
                20.0 * std::log10(rms_before / rms_after);

            printf("RMS before subtract = %.6f  after subtract = %.6f  ratio = %.3f dB\n", rms_before, rms_after, cancellation_db);
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

  }
  free(decoded);
  free(decoded_hashtable);

  monitor_free(&mon);
  return 0; // Caller frees signal
}
