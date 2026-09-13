/*
 * ============================================================
 * FT8 DDC V2
 * ============================================================
 *
 * 12 kHz real
 *      |
 *      | complex mixer @ 1450 Hz
 *      v
 *  FIR low-pass 110 Hz / 121 taps
 *      |
 *      | decimate by 60
 *      v
 *  250 Hz complex IQ
 *      |
 *      v
 *  FFT40 -> 160 ms / symbol
 *
 * V2 TEST:
 *
 *   180000 input samples
 *        ->
 *   3750 complex IQ samples
 *
 *   FFT40 every 40 IQ samples
 *
 *   Costas verification:
 *
 *       symbols 0..6
 *       symbols 36..42
 *       symbols 72..78
 *
 *   Expected residual frequency:
 *
 *       1500 - 1450 = 50 Hz
 *
 *   FFT bin:
 *
 *       50 / (200/32) = 8
 *
 *   Costas bins:
 *
 *       8 + {3,1,4,0,6,5,2}
 *       =
 *       {11,9,12,8,14,13,10}
 *
 * IMPORTANT:
 *
 *   The FIR delay is reported separately.
 *
 *   The FT8 verification windows remain referenced to the
 *   original signal start and the decimation grid:
 *
 *       2.000 ... 2.160 s
 *       7.760 ... 7.920 s
 *      13.520 ... 13.680 s
 *
 * ============================================================
 */

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

extern double elapsed_ms(const struct timespec *t0, const struct timespec *t1);

// ============================================================
// PARAMETERS
// ============================================================

static constexpr int INPUT_FS  = 12000;
static constexpr int OUTPUT_FS = 250;
static constexpr int DECIM     = 48;

static constexpr float LO_FREQ    = 1450.0f;
static constexpr float FIR_CUTOFF = 110.0f;
static constexpr int FIR_TAPS     = 121;

static constexpr int FFT_N = 40;

// static constexpr double SIGNAL_START_SEC = 2.000;

static constexpr int SYMBOL_SAMPLES =
    OUTPUT_FS * 160 / 1000;       // 40

// static constexpr int SIGNAL_START_INPUT =
//     static_cast<int>(
//         SIGNAL_START_SEC * INPUT_FS);

static constexpr int FIR_DELAY =
    (FIR_TAPS - 1) / 2;            // 60 input samples

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


// ============================================================
// FIR DESIGN
// ============================================================

static void design_fir(float *h)
{
    const double fc =
        static_cast<double>(FIR_CUTOFF) /
        static_cast<double>(INPUT_FS);

    const int M = FIR_TAPS - 1;
    const int mid = M / 2;

    double sum = 0.0;

    for (int n = 0; n < FIR_TAPS; ++n)
    {
        const int k = n - mid;

        double sinc;

        if (k == 0)
        {
            sinc = 2.0 * fc;
        }
        else
        {
            const double x =
                2.0 * M_PI *
                fc *
                static_cast<double>(k);

            sinc =
                std::sin(x) /
                (M_PI * static_cast<double>(k));
        }

        const double w =
            0.54 -
            0.46 *
            std::cos(
                2.0 * M_PI *
                static_cast<double>(n) /
                static_cast<double>(M));

        h[n] =
            static_cast<float>(
                sinc * w);

        sum += h[n];
    }

    for (int n = 0; n < FIR_TAPS; ++n)
    {
        h[n] =
            static_cast<float>(
                static_cast<double>(h[n]) /
                sum);
    }
}


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
    std::vector<IQ> &output)
{
    static ComplexCoeff hc[FIR_TAPS];
    static bool initialized = false;
    static float coefficient_lo_hz = 0.0f;

    // ============================================================
    // PRECOMPUTE COMPLEX FIR COEFFICIENTS
    // ============================================================

    /*
     * The complex FIR includes the mixer phasor, so it must be rebuilt when
     * processing a candidate at a different coarse frequency.  The filter
     * itself has only 121 taps, making this negligible beside the DDC work.
     */
    if (!initialized || std::fabs(lo_freq_hz - coefficient_lo_hz) > 1e-6f)
    {
        float h[FIR_TAPS];

        design_fir(h);

        const float omega =
            2.0f * static_cast<float>(M_PI) *
            // LO_FREQ /
            lo_freq_hz /
            static_cast<float>(INPUT_FS);

        const float cs = std::cos(omega);
        const float sn = std::sin(omega);

        float ph_re = 1.0f;
        float ph_im = 0.0f;

        for (int k = 0; k < FIR_TAPS; ++k)
        {
            hc[k].re = h[k] * ph_re;
            hc[k].im = h[k] * ph_im;

            const float nr =
                ph_re * cs - ph_im * sn;

            const float ni =
                ph_re * sn + ph_im * cs;

            ph_re = nr;
            ph_im = ni;
        }

        coefficient_lo_hz = lo_freq_hz;
        initialized = true;
    }

    // ============================================================
    // OUTPUT SIZE
    // ============================================================

    const int count =
        (num_samples + DECIM - 1) / DECIM;

    output.resize(count);

    int out_index = 0;

    /*
     * Mix the decimated samples down by the same LO used in the FIR.  The old
     * four-case rotation was valid only for a 1450 Hz LO, where the decimated
     * phase increment happens to be exactly -pi/2.
     */
    const float output_phase_step =
        -2.0f * static_cast<float>(M_PI) * lo_freq_hz *
        static_cast<float>(DECIM) / static_cast<float>(INPUT_FS);
    const float output_step_c = std::cos(output_phase_step);
    const float output_step_s = std::sin(output_phase_step);
    float output_phase_c = 1.0f;
    float output_phase_s = 0.0f;

    // ============================================================
    // DECIMATED FIR
    // ============================================================

    for (int n = 0; n < num_samples; n += DECIM)
    {
        float acc_re = 0.0f;
        float acc_im = 0.0f;

        // --------------------------------------------------------
        // Startup section
        //
        // Only the first 60 output positions need boundary
        // checking. After that all 121 taps are valid.
        // --------------------------------------------------------

        if (n < FIR_TAPS)
        {
            int k = 0;

            for (; k + 3 < FIR_TAPS && k + 3 <= n; k += 4)
            {
                const float s0 = input[n - k];
                const float s1 = input[n - k - 1];
                const float s2 = input[n - k - 2];
                const float s3 = input[n - k - 3];

                acc_re +=
                    s0 * hc[k].re +
                    s1 * hc[k + 1].re +
                    s2 * hc[k + 2].re +
                    s3 * hc[k + 3].re;

                acc_im +=
                    s0 * hc[k].im +
                    s1 * hc[k + 1].im +
                    s2 * hc[k + 2].im +
                    s3 * hc[k + 3].im;
            }

            for (; k <= n && k < FIR_TAPS; ++k)
            {
                const float s = input[n - k];

                acc_re += s * hc[k].re;
                acc_im += s * hc[k].im;
            }
        }
        else
        {
            // ----------------------------------------------------
            // FULL FIR
            //
            // No boundary test.
            // 121 taps = 30 groups of 4 + 1 tap.
            // ----------------------------------------------------

            int k = 0;

            for (; k < 120; k += 4)
            {
                const float s0 = input[n - k];
                const float s1 = input[n - k - 1];
                const float s2 = input[n - k - 2];
                const float s3 = input[n - k - 3];

                acc_re +=
                    s0 * hc[k].re +
                    s1 * hc[k + 1].re +
                    s2 * hc[k + 2].re +
                    s3 * hc[k + 3].re;

                acc_im +=
                    s0 * hc[k].im +
                    s1 * hc[k + 1].im +
                    s2 * hc[k + 2].im +
                    s3 * hc[k + 3].im;
            }

            // Tap 120
            const float s = input[n - 120];

            acc_re += s * hc[120].re;
            acc_im += s * hc[120].im;
        }

        output[out_index].i =
            acc_re * output_phase_c - acc_im * output_phase_s;
        output[out_index].q =
            acc_re * output_phase_s + acc_im * output_phase_c;

        const float next_phase_c =
            output_phase_c * output_step_c -
            output_phase_s * output_step_s;
        output_phase_s =
            output_phase_s * output_step_c +
            output_phase_c * output_step_s;
        output_phase_c = next_phase_c;

        ++out_index;
    }
}


// ============================================================
// FINE SYNC ESTIMATION (frequency + timing)
// ============================================================
//
// Once an FT8 message has already been decoded, its coarse
// residual frequency (base_bin) and coarse symbol timing
// (signal_start_output) are known from the decoder. This
// estimator refines both, using only the 21 Costas sync
// symbols, so the decoded signal can later be regenerated and
// subtracted (in the time domain) from the raw IQ stream.
//
// Frequency:
//
//   A direct DFT (single-bin Goertzel-style) value is evaluated
//   at the *expected* Costas bin for each of the 21 sync
//   symbols. Consecutive Costas symbols are exactly FFT_N output
//   samples apart, so any residual (sub-bin) frequency offset
//   appears as a constant phase rotation between consecutive
//   symbols:
//
//       phase_step = 2*pi * df * (FFT_N / OUTPUT_FS)
//
//   df is recovered with a standard phase-difference ("Kay")
//   estimator: sum X[s] * conj(X[s-1]) over every consecutive
//   in-block symbol pair, then take the angle of the resulting
//   complex sum.
//
// Timing:
//
//   Using the corrected frequency, the total sync energy
//   (sum of |X|^2 at the expected bins, de-rotated by the fine
//   frequency estimate) is evaluated at the coarse start and at
//   +-1 output sample. A parabolic interpolation across these
//   three energy values gives the sub-sample timing offset that
//   maximizes coherence, i.e. the true message start time.
//
// ============================================================

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

static DftTwiddle dft_twiddles[FFT_N][FFT_N];
static bool dft_twiddles_ready = false;

static void initialize_dft_twiddles()
{
    if (dft_twiddles_ready)
        return;

    for (int bin = 0; bin < FFT_N; ++bin)
    {
        const float step =
            -2.0f * static_cast<float>(M_PI) *
            static_cast<float>(bin) /
            static_cast<float>(FFT_N);

        const float cs = std::cos(step);
        const float sn = std::sin(step);

        float ph_re = 1.0f;
        float ph_im = 0.0f;

        for (int n = 0; n < FFT_N; ++n)
        {
            dft_twiddles[bin][n] = DftTwiddle{ph_re, ph_im};

            const float next_re = ph_re * cs - ph_im * sn;
            const float next_im = ph_re * sn + ph_im * cs;

            ph_re = next_re;
            ph_im = next_im;
        }
    }

    dft_twiddles_ready = true;
}

// Direct DFT value at one exact (non-searched) bin.
static IQ dft_bin(
    const std::vector<IQ> &x,
    int start,
    int bin)
{
    initialize_dft_twiddles();

    float re = 0.0f;
    float im = 0.0f;

    // DFT bins are periodic modulo FFT_N.  The expected Costas bins
    // are normally 0..31, but this also handles an edge bin safely.
    bin %= FFT_N;
    if (bin < 0)
        bin += FFT_N;

    for (int n = 0; n < FFT_N; ++n)
    {
        const DftTwiddle &w = dft_twiddles[bin][n];

        re += x[start + n].i * w.c - x[start + n].q * w.s;
        im += x[start + n].i * w.s + x[start + n].q * w.c;
    }

    return IQ{re, im};
}

// Sum of |X|^2 at the expected Costas bins, for a given
// output-sample start shift.
//
// A common phase rotation does not affect |X|^2, so frequency is
// deliberately not an argument here.  Fine frequency is still
// estimated below from phase differences between Costas symbols.
static float costas_sync_energy(
    const std::vector<IQ> &x,
    int signal_start_output,
    int base_bin,
    const int *costas,
    int costas_len,
    const int *costas_start,
    int costas_arrays,
    int start_shift)
{
    float energy = 0.0f;

    for (int array = 0; array < costas_arrays; ++array)
    {
        const int first = costas_start[array];

        for (int s = 0; s < costas_len; ++s)
        {
            const int symbol = first + s;

            const int fft_start =
                signal_start_output + start_shift +
                symbol * FFT_N;

            if (fft_start < 0 ||
                fft_start + FFT_N > static_cast<int>(x.size()))
                continue;

            const int bin = base_bin + costas[s];

            const IQ X = dft_bin(x, fft_start, bin);

            energy += X.i * X.i + X.q * X.q;
        }
    }

    return energy;
}

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
    float *out_fine_delay_time_ms)
{
    struct timespec t0;
    struct timespec t1;

    // --------------------------------------------------------
    // 0) Coarse timing search.
    //
    // The caller's signal_start_output may be off by more than
    // one output sample (e.g. tens of ms), so first do an
    // integer-sample grid search for the shift that maximizes
    // Costas sync energy (frequency assumed 0 for this step;
    // a few Hz of residual error only mildly attenuates the
    // energy metric and does not bias the peak location).
    // --------------------------------------------------------

    static constexpr int COARSE_SEARCH_RANGE = 20;  // +-80 ms @ 250 Hz

    clock_gettime(CLOCK_MONOTONIC, &t0);

    int best_shift = 0;
    float best_energy = -1.0f;

    for (int shift = -COARSE_SEARCH_RANGE;
         shift <= COARSE_SEARCH_RANGE;
         ++shift)
    {
        const float e =
            costas_sync_energy(
                output, signal_start_output, base_bin,
                costas, costas_len, costas_start, costas_arrays,
                shift);

        if (e > best_energy)
        {
            best_energy = e;
            best_shift = shift;
        }
    }

    const int aligned_start =
        signal_start_output + best_shift;

    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (out_coarse_delay_time_ms)
        *out_coarse_delay_time_ms = static_cast<float>(elapsed_ms(&t0, &t1));

    // --------------------------------------------------------
    // 1) Fine frequency: averaged phase difference between
    //    consecutive Costas symbols at their expected bins,
    //    evaluated at the coarse-aligned start.
    // --------------------------------------------------------

    clock_gettime(CLOCK_MONOTONIC, &t0);

    float sum_re = 0.0f;
    float sum_im = 0.0f;

    for (int array = 0; array < costas_arrays; ++array)
    {
        const int first = costas_start[array];

        IQ prev{0.0f, 0.0f};
        bool have_prev = false;

        for (int s = 0; s < costas_len; ++s)
        {
            const int symbol = first + s;

            const int fft_start =
                aligned_start + symbol * FFT_N;

            if (fft_start < 0 ||
                fft_start + FFT_N > static_cast<int>(output.size()))
            {
                have_prev = false;
                continue;
            }

            const int bin = base_bin + costas[s];

            const IQ X = dft_bin(output, fft_start, bin);

            if (have_prev)
            {
                // X * conj(prev)
                sum_re += X.i * prev.i + X.q * prev.q;
                sum_im += X.q * prev.i - X.i * prev.q;
            }

            prev = X;
            have_prev = true;
        }
    }

    const float phase_step =
        std::atan2(sum_im, sum_re);

    const float symbol_period =
        static_cast<float>(FFT_N) /
        static_cast<float>(OUTPUT_FS);

    const float freq_offset_hz =
        phase_step /
        (2.0f * static_cast<float>(M_PI) * symbol_period);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (out_fine_frequency_time_ms)
        *out_fine_frequency_time_ms = static_cast<float>(elapsed_ms(&t0, &t1));

    // --------------------------------------------------------
    // 2) Fine timing: parabolic interpolation of sync energy
    //    around the coarse-aligned (integer output-sample)
    //    start, now using the recovered fine frequency.
    // --------------------------------------------------------

    clock_gettime(CLOCK_MONOTONIC, &t0);

    const float e_minus =
        costas_sync_energy(
            output, aligned_start, base_bin,
            costas, costas_len, costas_start, costas_arrays,
            -1);

    const float e_zero =
        costas_sync_energy(
            output, aligned_start, base_bin,
            costas, costas_len, costas_start, costas_arrays,
            0);

    const float e_plus =
        costas_sync_energy(
            output, aligned_start, base_bin,
            costas, costas_len, costas_start, costas_arrays,
            +1);

    const float denom =
        e_minus - 2.0f * e_zero + e_plus;

    float subsample_offset = 0.0f;

    if (std::fabs(denom) > 1e-9f)
    {
        subsample_offset =
            0.5f * (e_minus - e_plus) / denom;

        // Parabolic fit is only valid within +-1 sample.
        if (subsample_offset > 1.0f)  subsample_offset = 1.0f;
        if (subsample_offset < -1.0f) subsample_offset = -1.0f;
    }

    *out_freq_offset_hz = freq_offset_hz;
    *out_time_offset_samples =
        static_cast<float>(best_shift) + subsample_offset;

    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (out_fine_delay_time_ms)
        *out_fine_delay_time_ms = static_cast<float>(elapsed_ms(&t0, &t1));
}

// ============================================================
// BLIND COSTAS ACQUISITION AND COMPLEX-IQ SUBTRACTION
// ============================================================
//
// Adapted from ddc_ft8_blind.cpp for this decoder's existing DDC format:
// 250 complex samples/s, 40 samples per FT8 symbol, and a 60-sample
// input-side FIR delay.  Keeping this here lets the decoded-message DDC
// residual be searched without creating a second, incompatible IQ stream.

static constexpr int DDC_FT8_NTONES = 79;
static constexpr int DDC_FT8_SPS = 40;
// The final Costas symbol is symbol 74.  The remaining four FT8 symbols are
// data, so acquisition can accept a message that is clipped at the slot tail.
static constexpr int DDC_FT8_SYNC_SPAN_SAMPLES = 75 * DDC_FT8_SPS;

struct CostasCandidate
{
    float delay_s;
    float freq_hz;
    float score;
};

static float ddc_full_costas_score_200(
    const std::vector<IQ>& iq,
    int start,
    float freq_hz,
    float ddc_lo_hz)
{
    float total = 0.0f;

    for (int array = 0; array < 3; ++array)
    {
        const int first_symbol = COSTAS_START[array];
        for (int symbol = 0; symbol < 7; ++symbol)
        {
            const int p0 = start + (first_symbol + symbol) * DDC_FT8_SPS;
            if (p0 < 0 || p0 + DDC_FT8_SPS > (int)iq.size())
                continue;

            const float w =
                2.0f * static_cast<float>(M_PI) *
                ((freq_hz - ddc_lo_hz) +
                 6.25f * (float)COSTAS[symbol]) /
                (float)OUTPUT_FS;
            const float cd = cosf(w);
            const float sd = sinf(w);
            float c = 1.0f;
            float s = 0.0f;
            float ci = 0.0f;
            float cq = 0.0f;

            for (int n = 0; n < DDC_FT8_SPS; ++n)
            {
                const IQ z = iq[p0 + n];
                ci += z.i * c + z.q * s;
                cq += z.q * c - z.i * s;

                const float nc = c * cd - s * sd;
                s = s * cd + c * sd;
                c = nc;
            }

            total += ci * ci + cq * cq;
        }
    }

    return total / 21.0f;
}

CostasCandidate refine_costas_candidate_200(
    const std::vector<IQ>& iq,
    int iq_first_input,
    float ddc_lo_hz,
    const CostasCandidate& candidate)
{
    constexpr float frequency_range_hz = 2.5f;
    constexpr float frequency_step_hz = 0.25f;
    constexpr int time_range_samples = 2;

    const float iq_start_s = (float)iq_first_input / (float)INPUT_FS;
    const int coarse_start = (int)lroundf(
        (candidate.delay_s - iq_start_s) * (float)OUTPUT_FS);
    const int last_start = (int)iq.size() - DDC_FT8_SYNC_SPAN_SAMPLES;
    CostasCandidate best = {candidate.delay_s, candidate.freq_hz, -1.0f};

    for (int time_offset = -time_range_samples;
         time_offset <= time_range_samples;
         ++time_offset)
    {
        const int start = coarse_start + time_offset;
        if (start < 0 || start > last_start)
            continue;

        for (float freq_hz = candidate.freq_hz - frequency_range_hz;
             freq_hz <= candidate.freq_hz + frequency_range_hz + 0.001f;
             freq_hz += frequency_step_hz)
        {
            const float score = ddc_full_costas_score_200(
                iq, start, freq_hz, ddc_lo_hz);
            if (score > best.score)
            {
                best.delay_s = iq_start_s +
                    (float)start / (float)OUTPUT_FS;
                best.freq_hz = freq_hz;
                best.score = score;
            }
        }
    }

    return best;
}

std::vector<CostasCandidate> find_costas_candidates_200(
    const std::vector<IQ>& iq,
    int iq_first_input,
    float ddc_lo_hz,
    float min_freq_hz,
    float max_freq_hz)
{
    constexpr float freq_step_hz = 1.0f;
    // The later Costas refinement searches +/-2 samples, so a 4-sample
    // (16 ms) acquisition grid still covers every possible symbol start
    // while halving the dominant blind-search loop.
    constexpr int time_step_samples = 4;
    // Keep more cheap 9-symbol hypotheses than we will decode.  Their
    // ordering can change substantially after the full 21-symbol rescore.
    constexpr int coarse_max_candidates = 64;
    constexpr int final_max_candidates = 16;
    constexpr float min_freq_separation_hz = 3.0f;
    constexpr float min_time_separation_s = 0.080f;
    static constexpr int sync_symbols[9] = {0, 1, 2, 36, 37, 38, 72, 73, 74};
    static constexpr int sync_tones[9] = {3, 1, 4, 3, 1, 4, 3, 1, 4};

    std::vector<CostasCandidate> best;
    best.reserve(coarse_max_candidates);
    const int last_start = (int)iq.size() - DDC_FT8_SYNC_SPAN_SAMPLES;
    const float iq_start_s = (float)iq_first_input / (float)INPUT_FS;
    if (last_start < 0 || min_freq_hz > max_freq_hz)
        return best;

    for (float freq_hz = min_freq_hz;
         freq_hz <= max_freq_hz + 0.001f;
         freq_hz += freq_step_hz)
    {
        float osc_i[7][DDC_FT8_SPS];
        float osc_q[7][DDC_FT8_SPS];

        for (int tone = 0; tone < 7; ++tone)
        {
            const float w =
                2.0f * static_cast<float>(M_PI) *
                ((freq_hz - ddc_lo_hz) + 6.25f * (float)tone) /
                (float)OUTPUT_FS;
            const float cd = cosf(w);
            const float sd = sinf(w);
            float c = 1.0f;
            float s = 0.0f;

            for (int n = 0; n < DDC_FT8_SPS; ++n)
            {
                osc_i[tone][n] = c;
                osc_q[tone][n] = s;
                const float nc = c * cd - s * sd;
                s = s * cd + c * sd;
                c = nc;
            }
        }

        float ci[9];
        float cq[9];
        int pos[9];
        float pc[9];
        float ps[9];
        float ec[9];
        float es[9];

        for (int index = 0; index < 9; ++index)
        {
            const int p0 = sync_symbols[index] * DDC_FT8_SPS;
            const int tone = sync_tones[index];
            float sum_i = 0.0f;
            float sum_q = 0.0f;

            for (int n = 0; n < DDC_FT8_SPS; ++n)
            {
                const IQ z = iq[p0 + n];
                sum_i += z.i * osc_i[tone][n] + z.q * osc_q[tone][n];
                sum_q += z.q * osc_i[tone][n] - z.i * osc_q[tone][n];
            }

            ci[index] = sum_i;
            cq[index] = sum_q;
            pos[index] = p0;

            const float w =
                2.0f * static_cast<float>(M_PI) *
                ((freq_hz - ddc_lo_hz) +
                 6.25f * (float)tone) /
                (float)OUTPUT_FS;
            pc[index] = cosf(w);
            ps[index] = sinf(w);
            ec[index] = cosf(-w * (float)DDC_FT8_SPS);
            es[index] = sinf(-w * (float)DDC_FT8_SPS);
        }

        for (int start = 0; start <= last_start; start += time_step_samples)
        {
            CostasCandidate candidate = {
                iq_start_s + (float)start / (float)OUTPUT_FS,
                freq_hz,
                0.0f};

            for (int i = 0; i < 9; ++i)
                candidate.score += ci[i] * ci[i] + cq[i] * cq[i];
            candidate.score /= 9.0f;

            int overlap_index = -1;
            for (size_t i = 0; i < best.size(); ++i)
            {
                if (fabsf(candidate.freq_hz - best[i].freq_hz) <
                        min_freq_separation_hz &&
                    fabsf(candidate.delay_s - best[i].delay_s) <
                        min_time_separation_s)
                {
                    overlap_index = (int)i;
                    break;
                }
            }

            int changed_index = -1;
            if (overlap_index >= 0)
            {
                if (candidate.score > best[overlap_index].score)
                {
                    best[overlap_index] = candidate;
                    changed_index = overlap_index;
                }
            }
            else if ((int)best.size() < coarse_max_candidates)
            {
                best.push_back(candidate);
                changed_index = (int)best.size() - 1;
            }
            else if (candidate.score > best.back().score)
            {
                best.back() = candidate;
                changed_index = (int)best.size() - 1;
            }

            // The list is already sorted.  A grid point can change only one
            // entry, so bubbling that entry to its correct position is
            // equivalent to sorting the complete 64-entry list, but avoids
            // doing thousands of full std::sort() calls per blind scan.
            if (changed_index >= 0)
            {
                while (changed_index > 0 &&
                       best[changed_index].score >
                           best[changed_index - 1].score)
                {
                    std::swap(best[changed_index],
                              best[changed_index - 1]);
                    --changed_index;
                }
                while (changed_index + 1 < (int)best.size() &&
                       best[changed_index + 1].score >
                           best[changed_index].score)
                {
                    std::swap(best[changed_index],
                              best[changed_index + 1]);
                    ++changed_index;
                }
            }

            if (start + time_step_samples > last_start)
                continue;

            for (int step = 0; step < time_step_samples; ++step)
            {
                for (int i = 0; i < 9; ++i)
                {
                    const IQ old_sample = iq[pos[i]];
                    const IQ new_sample = iq[pos[i] + DDC_FT8_SPS];
                    const float ti = ci[i] - old_sample.i +
                        new_sample.i * ec[i] - new_sample.q * es[i];
                    const float tq = cq[i] - old_sample.q +
                        new_sample.i * es[i] + new_sample.q * ec[i];
                    ci[i] = ti * pc[i] - tq * ps[i];
                    cq[i] = ti * ps[i] + tq * pc[i];
                    ++pos[i];
                }
            }
        }
    }

    for (CostasCandidate& candidate : best)
    {
        const int start = (int)lroundf(
            (candidate.delay_s - iq_start_s) * (float)OUTPUT_FS);
        candidate.score = ddc_full_costas_score_200(
            iq, start, candidate.freq_hz, ddc_lo_hz);
    }

    std::sort(best.begin(), best.end(),
        [](const CostasCandidate& a, const CostasCandidate& b) {
            return a.score > b.score;
        });
    if ((int)best.size() > final_max_candidates)
        best.resize(final_max_candidates);
    return best;
}

struct DdcIqModelParams
{
    float amplitude;
    float phase;
};

static IQ ddc_estimate_symbol_coefficient_200(
    const std::vector<IQ>& iq,
    int start,
    int tone,
    float final_freq_hz,
    float ddc_lo_hz)
{
    IQ out = {0.0f, 0.0f};
    if (start < 0 || start >= (int)iq.size())
        return out;

    const int first = DDC_FT8_SPS / 10;
    // A message may extend a little beyond the captured 15 s slot.  The
    // available part of its final symbol still provides a useful estimate.
    const int last = std::min(
        DDC_FT8_SPS - first, (int)iq.size() - start);
    const int corr_len = last - first;
    if (corr_len <= 0)
        return out;
    const float w =
        2.0f * static_cast<float>(M_PI) *
        ((final_freq_hz - ddc_lo_hz) +
         6.25f * (float)(tone & 7)) /
        (float)OUTPUT_FS;
    const float cd = cosf(w);
    const float sd = sinf(w);
    float c = cosf(w * (float)first);
    float s = sinf(w * (float)first);

    for (int n = first; n < last; ++n)
    {
        const IQ z = iq[start + n];
        out.i += z.i * c + z.q * s;
        out.q += z.q * c - z.i * s;

        const float nc = c * cd - s * sd;
        s = s * cd + c * sd;
        c = nc;
    }

    const float inv = 1.0f / (float)corr_len;
    out.i *= inv;
    out.q *= inv;
    return out;
}

static void ddc_synthesize_ft8_message_200(
    std::vector<IQ>& dst,
    int start0,
    const uint8_t tones[DDC_FT8_NTONES],
    const DdcIqModelParams model[DDC_FT8_NTONES],
    float baseband_freq_hz)
{
    int ramp = (int)lroundf((float)DDC_FT8_SPS * 0.11f);
    if (ramp < 1)
        ramp = 1;

    for (int symbol = 0; symbol < DDC_FT8_NTONES; ++symbol)
    {
        const float freq = baseband_freq_hz +
            6.25f * (float)(tones[symbol] & 7);
        const float dtheta =
            2.0f * static_cast<float>(M_PI) * freq / (float)OUTPUT_FS;
        const float cd = cosf(dtheta);
        const float sd = sinf(dtheta);
        const float theta_start = model[symbol].phase + ramp * dtheta;
        float c = cosf(theta_start);
        float s = sinf(theta_start);

        for (int n = ramp; n < DDC_FT8_SPS - ramp; ++n)
        {
            const int index = start0 + symbol * DDC_FT8_SPS + n;
            if (index >= 0 && index < (int)dst.size())
            {
                dst[index].i += model[symbol].amplitude * c;
                dst[index].q += model[symbol].amplitude * s;
            }

            const float nc = c * cd - s * sd;
            s = s * cd + c * sd;
            c = nc;
        }

        float theta = model[symbol].phase +
            (DDC_FT8_SPS - ramp) * dtheta;
        const float next_freq = (symbol + 1 == DDC_FT8_NTONES) ? freq :
            baseband_freq_hz + 6.25f * (float)(tones[symbol + 1] & 7);
        const float next_phase = (symbol + 1 == DDC_FT8_NTONES) ?
            model[symbol].phase : model[symbol + 1].phase;
        const float dtheta_next =
            2.0f * static_cast<float>(M_PI) * next_freq / (float)OUTPUT_FS;
        const float increment = (dtheta_next - dtheta) / (2.0f * ramp);
        const float actual = theta + dtheta * 2.0f * ramp +
            increment * 2.0f * ramp * ramp;
        const float target = next_phase + dtheta_next * ramp;
        const float adjustment = remainderf(
            target - actual, 2.0f * static_cast<float>(M_PI)) /
            (2.0f * ramp);

        float delta = dtheta + adjustment;
        float ct = cosf(theta);
        float st = sinf(theta);
        float rc = cosf(delta);
        float rs = sinf(delta);
        const float ric = cosf(increment);
        const float ris = sinf(increment);
        const int transition_end = (symbol + 1 == DDC_FT8_NTONES) ?
            DDC_FT8_SPS : DDC_FT8_SPS + ramp;

        for (int n = DDC_FT8_SPS - ramp; n < transition_end; ++n)
        {
            const int index = start0 + symbol * DDC_FT8_SPS + n;
            if (index >= 0 && index < (int)dst.size())
            {
                float gain = model[symbol].amplitude;
                if (symbol + 1 == DDC_FT8_NTONES)
                    gain *= 1.0f -
                        (float)(n - (DDC_FT8_SPS - ramp)) / (float)ramp;
                dst[index].i += gain * ct;
                dst[index].q += gain * st;
            }

            const float nct = ct * rc - st * rs;
            st = st * rc + ct * rs;
            ct = nct;

            const float nrc = rc * ric - rs * ris;
            rs = rs * ric + rc * ris;
            rc = nrc;
        }
    }
}

bool subtract_ft8_message_200(
    std::vector<IQ>& iq,
    int iq_first_input,
    const uint8_t tones[DDC_FT8_NTONES],
    float delay_s,
    float final_freq_hz,
    float ddc_lo_hz,
    float* rms_before,
    float* rms_after,
    float* rms_model,
    float* residual_projection_db)
{
    const float iq_start_s = (float)iq_first_input / (float)INPUT_FS;
    const int start0 = (int)lroundf(
        (delay_s - iq_start_s) * (float)OUTPUT_FS);
    if (start0 < 0 || start0 >= (int)iq.size())
        return false;

    // The three Costas sequences end at symbol 75.  Permit a clipped data
    // tail, but only subtract when all synchronization symbols are present.
    if ((int)iq.size() - start0 < DDC_FT8_SYNC_SPAN_SAMPLES)
        return false;

    DdcIqModelParams model[DDC_FT8_NTONES];
    for (int symbol = 0; symbol < DDC_FT8_NTONES; ++symbol)
    {
        const IQ coefficient = ddc_estimate_symbol_coefficient_200(
            iq, start0 + symbol * DDC_FT8_SPS, tones[symbol],
            final_freq_hz, ddc_lo_hz);
        model[symbol].amplitude = sqrtf(
            coefficient.i * coefficient.i + coefficient.q * coefficient.q);
        model[symbol].phase = atan2f(coefficient.q, coefficient.i);
    }

    static std::vector<IQ> reference;
    reference.assign(iq.size(), IQ{0.0f, 0.0f});
    ddc_synthesize_ft8_message_200(
        reference, start0, tones, model, final_freq_hz - ddc_lo_hz);

    double energy_before = 0.0;
    double energy_model = 0.0;
    double energy_after = 0.0;
    double cross_re = 0.0;
    double cross_im = 0.0;

    const int last = std::min(
        (int)iq.size(), start0 + DDC_FT8_NTONES * DDC_FT8_SPS);
    for (int index = start0; index < last; ++index)
    {
        const IQ x = iq[index];
        const IQ model_sample = reference[index];
        const double mr = model_sample.i;
        const double mq = model_sample.q;
        const double xr = x.i;
        const double xq = x.q;

        cross_re += mr * xr + mq * xq;
        cross_im += mr * xq - mq * xr;
        energy_model += mr * mr + mq * mq;
        energy_before += xr * xr + xq * xq;

        iq[index] = {
            x.i - model_sample.i,
            x.q - model_sample.q};
        energy_after += (double)iq[index].i * iq[index].i +
                        (double)iq[index].q * iq[index].q;
    }

    if (energy_model <= 0.0)
        return false;

    const double count = (double)(last - start0);
    if (rms_before)
        *rms_before = sqrtf((float)(energy_before / count));
    if (rms_model)
        *rms_model = sqrtf((float)(energy_model / count));
    if (rms_after)
        *rms_after = sqrtf((float)(energy_after / count));
    if (residual_projection_db)
    {
        const double residual_re = cross_re - energy_model;
        const double residual = sqrt(residual_re * residual_re +
                                     cross_im * cross_im) / energy_model;
        *residual_projection_db = residual > 0.0 ?
            (float)(20.0 * log10(residual)) : -INFINITY;
    }

    return true;
}



// ============================================================
// FFT40
// ============================================================

static int fft_peak_bin(
    const std::vector<IQ> &x,
    int start)
{
    float mag[FFT_N];


    for (int k = 0;
         k < FFT_N;
         ++k)
    {
        float re = 0.0f;
        float im = 0.0f;


        for (int n = 0;
             n < FFT_N;
             ++n)
        {
            const float angle =
                -2.0f *
                static_cast<float>(M_PI) *
                static_cast<float>(k*n) /
                static_cast<float>(FFT_N);


            const float c =
                std::cos(angle);

            const float s =
                std::sin(angle);


            re +=
                x[start+n].i * c -
                x[start+n].q * s;

            im +=
                x[start+n].i * s +
                x[start+n].q * c;
        }


        mag[k] =
            re*re + im*im;
    }


    int peak = 0;


    for (int k = 1;
         k < FFT_N;
         ++k)
    {
        if (mag[k] > mag[peak])
            peak = k;
    }


    return peak;
}
