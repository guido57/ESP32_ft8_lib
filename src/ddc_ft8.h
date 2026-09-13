
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <time.h>

#ifdef ARDUINO
#include <Arduino.h>
#endif


// ============================================================
// PARAMETERS
// ============================================================

// static constexpr int INPUT_FS  = 12000;
// static constexpr int OUTPUT_FS = 250;
// static constexpr int DECIM     = 48;

// static constexpr float LO_FREQ    = 1450.0f;
// static constexpr float FIR_CUTOFF = 110.0f;
// static constexpr int FIR_TAPS     = 121;

// static constexpr int FFT_N = 40;

// static constexpr double SIGNAL_START_SEC = 2.000;

// static constexpr int SYMBOL_SAMPLES =
//     OUTPUT_FS * 160 / 1000;       // 40

// static constexpr int SIGNAL_START_INPUT =
//     static_cast<int>(
//         SIGNAL_START_SEC * INPUT_FS);

// static constexpr int FIR_DELAY =
//     (FIR_TAPS - 1) / 2;            // 60 input samples

/*
 * IMPORTANT:
 *
 * Timing of the FFT test is tied to the original signal
 * start and the decimation grid.
 *
 * 24000 / 60 = 400
 *
 * therefore:
 *
 * FFT #0 starts at 400 * 5 ms = 2.000 s
 */
// static constexpr int SIGNAL_START_OUTPUT =
//     SIGNAL_START_INPUT / DECIM;


// ============================================================
// COSTAS
// ============================================================

static constexpr int COSTAS_LEN = 7;

static constexpr int COSTAS[COSTAS_LEN] =
{
    3, 1, 4, 0, 6, 5, 2
};

static constexpr int COSTAS_START[3] =
{
    0, 36, 72
};


// ============================================================
// COMPLEX SAMPLE
// ============================================================

struct IQ
{
    float i;
    float q;
};

/// A blind FT8 Costas acquisition result in raw-input time coordinates.
struct CostasCandidate
{
    float delay_s;
    float freq_hz;
    float score;
};



// ============================================================
// FIR DESIGN
// ============================================================

// static void design_fir(float *h)
// {
//     const double fc =
//         static_cast<double>(FIR_CUTOFF) /
//         static_cast<double>(INPUT_FS);

//     const int M = FIR_TAPS - 1;
//     const int mid = M / 2;

//     double sum = 0.0;

//     for (int n = 0; n < FIR_TAPS; ++n)
//     {
//         const int k = n - mid;

//         double sinc;

//         if (k == 0)
//         {
//             sinc = 2.0 * fc;
//         }
//         else
//         {
//             const double x =
//                 2.0 * M_PI *
//                 fc *
//                 static_cast<double>(k);

//             sinc =
//                 std::sin(x) /
//                 (M_PI * static_cast<double>(k));
//         }

//         const double w =
//             0.54 -
//             0.46 *
//             std::cos(
//                 2.0 * M_PI *
//                 static_cast<double>(n) /
//                 static_cast<double>(M));

//         h[n] =
//             static_cast<float>(
//                 sinc * w);

//         sum += h[n];
//     }

//     for (int n = 0; n < FIR_TAPS; ++n)
//     {
//         h[n] =
//             static_cast<float>(
//                 static_cast<double>(h[n]) /
//                 sum);
//     }
// }


// ============================================================
// FAST DDC V2
// ============================================================
//
// Only decimated output samples are calculated.
//
// 3000 output samples
// 121 FIR taps
//
// = 363000 FIR coefficients
//
// No complete mixed input buffer is generated.
//
// ============================================================

struct ComplexCoeff
{
    float re;
    float im;
};

void ddc_process_fast(float lo_freq_hz,
    const float *__restrict input,
    int num_samples,
    std::vector<IQ> &output);

struct SyncEstimate
{
    float freq_offset_hz;    // fine correction, added to base_bin frequency
    float time_offset_samp;  // fine correction, added to signal_start_output (output samples)
};

// Twiddles for the exact 40-point DFT bins used by Costas sync.
// They are generated once with a phasor recurrence, avoiding 32
// sin()/cos() pairs every time a Costas symbol is evaluated.
struct DftTwiddle
{
    float c;
    float s;
};


static float costas_sync_energy(
    const std::vector<IQ> &x,
    int signal_start_output,
    int base_bin,
    const int *costas,
    int costas_len,
    const int *costas_start,
    int costas_arrays,
    int start_shift);
    
// ------------------------------------------------------------
// Public entry point.
//
// costas[costas_len]        : Costas bin offsets (e.g. {3,1,4,0,6,5,2})
// costas_start[costas_arrays]: first symbol index of each Costas array
// ------------------------------------------------------------

void ddc_estimate_fine_sync(
    const std::vector<IQ> &output,
    int signal_start_output,
    int base_bin,
    const int *costas,
    int costas_len,
    const int *costas_start,
    int costas_arrays,
    float *out_freq_offset_hz,
    float *out_time_offset_samples,
    float *out_coarse_delay_time_ms,
    float *out_fine_frequency_time_ms,
    float *out_fine_delay_time_ms);

/**
 * Remove one decoded FT8 message from the 250 Hz complex DDC stream.
 * iq_first_input is the raw-input sample time represented by iq[0]; for the
 * full-stream DDC this is -60 because of its 121-tap FIR group delay.
 */
bool subtract_ft8_message_200(
    std::vector<IQ> &iq,
    int iq_first_input,
    const uint8_t tones[79],
    float delay_s,
    float final_freq_hz,
    float ddc_lo_hz,
    float *rms_before,
    float *rms_after,
    float *rms_model,
    float *residual_projection_db);

/**
 * Find FT8 Costas candidates in an already-decimated 250 Hz complex stream.
 * min_freq_hz and max_freq_hz are absolute tone-0 frequencies.
 */
std::vector<CostasCandidate> find_costas_candidates_200(
    const std::vector<IQ> &iq,
    int iq_first_input,
    float ddc_lo_hz,
    float min_freq_hz,
    float max_freq_hz);

/**
 * Refine a blind Costas candidate on the 250-Hz IQ grid.  This searches the
 * full 21-symbol Costas pattern around the candidate, so it also works when
 * the tone-0 frequency lies between the 6.25-Hz FFT bins.
 */
CostasCandidate refine_costas_candidate_200(
    const std::vector<IQ> &iq,
    int iq_first_input,
    float ddc_lo_hz,
    const CostasCandidate &candidate);
    
