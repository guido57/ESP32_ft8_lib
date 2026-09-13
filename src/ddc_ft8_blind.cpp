/*
 * ============================================================
 * FT8 V9.2 - ONE DDC -> 200 Hz + FINE FREQ/TIME SYNC
 * ============================================================
 *
 * Linux / Ubuntu + ESP32-S3 Arduino / LittleFS
 *
 * V9.4 architecture:
 *
 *   12 kHz real WAV
 *       |
 *       | ONE DDC at the message coarse frequency
 *       | windowed to the message time region
 *       | FIR LP + decimate 1:60
 *       v
 *   200 Hz complex IQ
 *       |
 *       | fine frequency + delay from the four waterfall hypotheses
 *       | (TEST HARNESS ONLY: same IQ is reused for all 4 hypotheses)
 *       v
 *   final timing refinement on all 79 symbols
 *
 * IMPORTANT:
 *   In real operation the DDC is performed ONCE per message.
 *   The four corners below are only alternative waterfall hypotheses
 *   fed to the fine synchronizer; they do NOT trigger another DDC.
 *
 * TEST WAV:
 *   12 kHz mono PCM16, known FT8 message at 1500 Hz, true delay 2.000 s.
 *
 * DDC test configuration:
 *   DDC_LO = 1498 Hz (one single DDC pass).
 *   The four waterfall hypotheses are then tested against the same 200-Hz IQ:
 *      1.950 s / 1498 Hz
 *      2.050 s / 1498 Hz
 *      1.950 s / 1502 Hz
 *      2.050 s / 1502 Hz
 *
 * This deliberately decouples the DDC LO from the hypothesis metadata so
 * that one DDC can validate all four fine-sync corner hypotheses.
 *
 * Fine synchronization:
 *   - 9 Costas symbols
 *   - absolute frequency +/-2 Hz around each waterfall frequency
 *   - 0.1 Hz grid + parabolic interpolation
 *   - delay +/-50 ms at 200 Hz (5 ms/sample)
 *   - final timing: all 79 symbols, +/-2 samples, parabolic interpolation
 *
 * DDC optimization:
 *   Only the required message window is processed.  For the 4 tests,
 *   the window is chosen to cover the full +/-50 ms timing uncertainty,
 *   all 79 symbols, and FIR boundary samples.
 * ============================================================
 */

#if defined(ARDUINO)
#include <Arduino.h>
#include <LittleFS.h>
#else
#include <chrono>
#include <fstream>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifndef PI
#define PI 3.14159265358979323846f
#endif
#ifndef TWO_PI
#define TWO_PI 6.28318530717958647692f
#endif

static constexpr int INPUT_FS  = 12000;
// 200-Hz RF subband, sampled at 250 Hz complex IQ to retain transition-band
// guard at the edges.  Adjacent subbands are stepped by 150 Hz.
static constexpr int OUTPUT_FS = 250;
static constexpr int DECIM = 48;
static constexpr int NTONES = 79;
static constexpr int SPS_IN = 1920;
static constexpr int SPS_200 = 40;
static constexpr float TONE_SPACING = 6.25f;
static constexpr float TEST_FREQ = 1500.0f;
static constexpr float DDC_LO_TEST = 1521.875f;
static constexpr float WEAK_TEST_DELAY_S = 1.500f;
static constexpr float WEAK_TEST_FREQ_HZ = 1515.0f;
static constexpr int FIR_TAPS = 121;
static constexpr int FIR_GD = (FIR_TAPS - 1) / 2;
static constexpr float FIR_CUTOFF_HZ = 110.0f;

static constexpr uint8_t SYNC_SYMBOLS[9] = {
    0, 1, 2, 36, 37, 38, 72, 73, 74
};

static constexpr uint8_t COSTAS_TONES[7] = {
    3, 1, 4, 0, 6, 5, 2
};

// Three leading symbols from each 7-symbol Costas block.  This is the
// inexpensive coarse-acquisition pattern; candidates are later verified
// with all 21 Costas symbols.
static constexpr uint8_t COARSE_COSTAS_TONES[9] = {
    3, 1, 4, 3, 1, 4, 3, 1, 4
};

static inline uint64_t bench_now_us()
{
#if defined(ARDUINO)
    return (uint64_t)micros();
#else
    return (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
#endif
}

struct IQ {
    float i;
    float q;
};

struct FineEstimate {
    float delay_s = 0.0f;
    float freq_hz = TEST_FREQ;
    float score = -1.0f;
    float integer_delay_s = 0.0f;
    float freq_grid_hz = TEST_FREQ;
    float delay_fraction_samples = 0.0f;
    uint64_t joint_us = 0;
    uint64_t final_us = 0;
};

struct CostasCandidate {
    float delay_s = 0.0f;
    float freq_hz = 0.0f;
    float score = -1.0f;
};

// ============================================================
// Fixed 121-tap Hamming-windowed low-pass FIR, cutoff 110 Hz at 12 kHz.
// Store its first half plus centre; the filter is exactly symmetric.
// Keeping these as constants removes FIR design (and its sinf/cosf calls)
// from the one-message DDC timing path.
// ============================================================
static constexpr float FIR_HALF[FIR_GD + 1] = {
    -1.4802712208e-04f, -1.2459101419e-04f, -1.0103317240e-04f, -7.6077328093e-05f, -4.8343324071e-05f, -1.6359536568e-05f,
     2.1423383255e-05f,  6.6617567669e-05f,  1.2088192388e-04f,  1.8590520483e-04f,  2.6338831700e-04f,  3.5502603274e-04f,
     4.6248828337e-04f,  5.8740121550e-04f,  7.3132819832e-04f,  8.9575097240e-04f,  1.0820511317e-03f,  1.2914921295e-03f,
     1.5252019967e-03f,  1.7841569555e-03f,  2.0691661053e-03f,  2.3808573495e-03f,  2.7196647205e-03f,  3.0858172501e-03f,
     3.4793295174e-03f,  3.8999939919e-03f,  4.3473752731e-03f,  4.8208063122e-03f,  5.3193866795e-03f,  5.8419829252e-03f,
     6.3872310605e-03f,  6.9535411649e-03f,  7.5391041062e-03f,  8.1419003379e-03f,  8.7597107204e-03f,  9.3901292909e-03f,
     1.0030577887e-02f,  1.0678322513e-02f,  1.1330491319e-02f,  1.1984094045e-02f,  1.2636042766e-02f,  1.3283173775e-02f,
     1.3922270399e-02f,  1.4550086561e-02f,  1.5163370891e-02f,  1.5758891150e-02f,  1.6333458781e-02f,  1.6883953337e-02f,
     1.7407346596e-02f,  1.7900726127e-02f,  1.8361318101e-02f,  1.8786509142e-02f,  1.9173867016e-02f,  1.9521159969e-02f,
     1.9826374543e-02f,  2.0087731701e-02f,  2.0303701111e-02f,  2.0473013473e-02f,  2.0594670751e-02f,  2.0667954238e-02f,
     2.0692430372e-02f
};

// ============================================================
// ONE DDC + FIR + decimate, limited to the required window
// ============================================================
static std::vector<IQ> ddc_decimate_200_window(const float* samples,
                                               int num_samples,
                                               float ddc_freq_hz,
                                               float window_start_s,
                                               float window_end_s,
                                               int* out_first_input)
{

    // Only process the input interval actually needed by the message,
    // including the FIR half-length on both sides.
    int first = (int)floorf(window_start_s * INPUT_FS) - FIR_GD;
    int last  = (int)ceilf(window_end_s * INPUT_FS) + FIR_GD;
    first = std::max(0, first);
    last = std::min(num_samples - 1, last);

    if (last < first) {
        if (out_first_input) *out_first_input = first - FIR_GD;
        return {};
    }

    // First decimated output: the FIR center is at n = FIR_GD + m*60.
    // IQ sample 0 therefore represents the corrected time n-FIR_GD.
    const int first_out_global = std::max(FIR_GD,
        ((first - FIR_GD + DECIM - 1) / DECIM) * DECIM + FIR_GD);
    const int nout = (last >= first_out_global) ?
        ((last - first_out_global) / DECIM + 1) : 0;

    std::vector<IQ> out;
    out.reserve((size_t)nout);

    // ------------------------------------------------------------
    // Direct decimated DDC/FIR.
    //
    // The previous V9.4 tried to use only 2-3 FIR taps per output.
    // That is not a valid polyphase decimator implementation: each
    // decimated output still contains all 121 FIR taps.  Here we keep
    // all taps, but evaluate ONLY the ~200-Hz output instants.
    //
    // For y[n] = sum_k h[k] x[n-k] exp(-j*w*(n-k)):
    //   y[n] = exp(-j*w*n) * sum_k h[k] x[n-k] exp(+j*w*k)
    //
    // Thus sin/cos for k are precomputed once and the 12-kHz mixer
    // does not have to run over every input sample.
    // ------------------------------------------------------------
    const float w = TWO_PI * ddc_freq_hz / (float)INPUT_FS;
    static float hc[FIR_TAPS];
    static float hs[FIR_TAPS];
    static float ck = 1.0f;
    static float sk = 0.0f;

    // Global-index independent FIR phase: exp(+j*w*k).
    const float cd = cosf(w);
    const float sd = sinf(w);
    float c = 1.0f;
    float s = 0.0f;
    for (int k = 0; k < FIR_TAPS; ++k) {
        const float h = FIR_HALF[(k <= FIR_GD) ? k : FIR_TAPS - 1 - k];
        hc[k] = h * c;
        hs[k] = h * s;
        const float nc = c * cd - s * sd;
        const float ns = s * cd + c * sd;
        c = nc;
        s = ns;
    }
    (void)ck;
    (void)sk;

    // Mixer phase for the first output's global input index n.
    float phase = w * (float)first_out_global;
    float co = cosf(phase);
    float si = sinf(phase);
    const float out_cd = cosf(w * (float)DECIM);
    const float out_sd = sinf(w * (float)DECIM);

    for (int m = 0; m < nout; ++m) {
        const int n = first_out_global + m * DECIM;

        float a = 0.0f; // real part of sum h*x*exp(+jwk)
        float b = 0.0f; // imag part of sum h*x*exp(+jwk)

        // 121 taps.  This is still much cheaper than mixing every
        // 12-kHz sample and then filtering every sample.
        for (int k = 0; k < FIR_TAPS; ++k) {
            const int input_index = n - k;
            const float x = (input_index >= 0 && input_index < num_samples) ?
                samples[input_index] : 0.0f;
            a += x * hc[k];
            b += x * hs[k];
        }

        // Multiply by exp(-j*w*n):
        //   I = a*cos + b*sin
        //   Q = b*cos - a*sin
        const float yi = a * co + b * si;
        const float yq = b * co - a * si;
        out.push_back({yi, yq});

        const float nco = co * out_cd - si * out_sd;
        const float nsi = si * out_cd + co * out_sd;
        co = nco;
        si = nsi;
    }

    if (out_first_input) *out_first_input = first_out_global - FIR_GD;
    return out;
}

// ============================================================
// Test helper: recover known FT8 tones from original WAV
// ============================================================
static uint8_t estimate_tone(const float* samples,
                             int num_samples,
                             int start,
                             float freq)
{
    float best_power = -1.0f;
    int best_tone = 0;
    if (start < 0 || start >= num_samples) return 0;

    const int count = std::min(SPS_IN, num_samples - start);

    for (int tone = 0; tone < 8; ++tone) {
        const float f = freq + tone * TONE_SPACING;
        const float w = TWO_PI * f / (float)INPUT_FS;
        const float cd = cosf(w);
        const float sd = sinf(w);
        float c = 1.0f, s = 0.0f;
        float ci = 0.0f, cq = 0.0f;

        for (int n = 0; n < count; ++n) {
            const float x = samples[start + n];
            ci += x * c;
            cq -= x * s;
            const float nc = c * cd - s * sd;
            const float ns = s * cd + c * sd;
            c = nc;
            s = ns;
        }

        const float p = ci * ci + cq * cq;
        if (p > best_power) {
            best_power = p;
            best_tone = tone;
        }
    }
    return (uint8_t)best_tone;
}

static bool extract_ft8_tones(const float* samples,
                              int num_samples,
                              float delay_s,
                              float freq_hz,
                              uint8_t tones[NTONES])
{
    const int start = (int)lroundf(delay_s * INPUT_FS);
    if (start < 0 || start >= num_samples) return false;

    for (int k = 0; k < NTONES; ++k)
        tones[k] = estimate_tone(samples, num_samples,
                                 start + k * SPS_IN, freq_hz);
    return true;
}

// ============================================================
// 200 Hz oscillator templates: three Costas tones, relative to DDC LO
// ============================================================
static void build_costas_templates(float absolute_freq_hz,
                                   float ddc_lo_hz,
                                   float ci[3][SPS_200],
                                   float si[3][SPS_200])
{
    const float base = absolute_freq_hz - ddc_lo_hz;
    const float tone_abs[3] = {
        base + TONE_SPACING * 3.0f,
        base + TONE_SPACING * 1.0f,
        base + TONE_SPACING * 4.0f
    };

    for (int t = 0; t < 3; ++t) {
        const float w = TWO_PI * tone_abs[t] / (float)OUTPUT_FS;
        const float cd = cosf(w);
        const float sd = sinf(w);
        float c = 1.0f;
        float s = 0.0f;
        for (int n = 0; n < SPS_200; ++n) {
            ci[t][n] = c;
            si[t][n] = s;
            const float nc = c * cd - s * sd;
            const float ns = s * cd + c * sd;
            c = nc;
            s = ns;
        }
    }
}

// ============================================================
// Fine joint frequency + delay at 200 Hz
// ============================================================
static FineEstimate fine_joint_200(const std::vector<IQ>& iq,
                                   int iq_first_input,
                                   const uint8_t tones[NTONES],
                                   float delay_guess_s,
                                   float coarse_freq_hz,
                                   float ddc_lo_hz)
{
    constexpr int DELAY_RADIUS = 10; // +/-50 ms at 200 Hz
    constexpr float FREQ_RADIUS = 2.0f;
    constexpr float FREQ_STEP = 0.1f;
    constexpr int NFREQ = 41;

    const uint64_t t0 = bench_now_us();

    FineEstimate r;
    r.delay_s = delay_guess_s;
    r.freq_hz = coarse_freq_hz;
    r.freq_grid_hz = coarse_freq_hz;

    // Convert the absolute coarse delay into the local 200-Hz IQ index.
    // FIR group delay is removed from the IQ time origin: the first output
    // corresponds to input time (first_out_global - GD), not first.
    const float window_start_s = (float)iq_first_input / (float)INPUT_FS;
    const int center = (int)lroundf((delay_guess_s - window_start_s) * OUTPUT_FS);
    const int first = std::max(0, center - DELAY_RADIUS);
    const int last = std::min((int)iq.size() - SPS_200 * NTONES,
                              center + DELAY_RADIUS);
    if (last < first) {
        r.joint_us = bench_now_us() - t0;
        return r;
    }

    const int ndelay = last - first + 1;
    float best_grid_scores[NFREQ] = {};

    for (int fi = 0; fi < NFREQ; ++fi) {
        const float delta = -FREQ_RADIUS + fi * FREQ_STEP;
        const float abs_freq = coarse_freq_hz + delta;

        float osc_i[3][SPS_200];
        float osc_q[3][SPS_200];
        build_costas_templates(abs_freq, ddc_lo_hz, osc_i, osc_q);

        float ci[9] = {};
        float cq[9] = {};
        int pos[9] = {};
        float pc[9] = {};
        float ps[9] = {};
        float ec[9] = {};
        float es[9] = {};
        bool valid[9] = {};

        for (int si = 0; si < 9; ++si) {
            const int k = SYNC_SYMBOLS[si];
            const int p0 = first + k * SPS_200;
            if (p0 < 0 || p0 + SPS_200 > (int)iq.size()) continue;

            const uint8_t tone = tones[k];
            int t = -1;
            if (tone == 3) t = 0;
            else if (tone == 1) t = 1;
            else if (tone == 4) t = 2;
            if (t < 0) continue;

            const float fbase = abs_freq - ddc_lo_hz;
            const float f = fbase + TONE_SPACING * (float)tone;
            const float w = TWO_PI * f / (float)OUTPUT_FS;
            pc[si] = cosf(w);
            ps[si] = sinf(w);
            ec[si] = cosf(-w * (float)SPS_200);
            es[si] = sinf(-w * (float)SPS_200);
            pos[si] = p0;
            valid[si] = true;

            float ii = 0.0f;
            float qq = 0.0f;
            for (int n = 0; n < SPS_200; ++n) {
                const float x = iq[p0 + n].i;
                const float y = iq[p0 + n].q;
                const float c = osc_i[t][n];
                const float s = osc_q[t][n];
                ii += x * c + y * s;
                qq += y * c - x * s;
            }
            ci[si] = ii;
            cq[si] = qq;
        }

        float best_score_f = -1.0f;
        int best_delay_idx_f = 0;

        for (int d = 0; d < ndelay; ++d) {
            float total = 0.0f;
            int used = 0;
            for (int si = 0; si < 9; ++si) {
                if (!valid[si]) continue;
                total += ci[si] * ci[si] + cq[si] * cq[si];
                ++used;
            }

            const float score = used ? total / (float)used : -1.0f;
            if (score > best_score_f) {
                best_score_f = score;
                best_delay_idx_f = d;
            }
            if (d == ndelay - 1) break;

            for (int si = 0; si < 9; ++si) {
                if (!valid[si]) continue;
                const int old_pos = pos[si];
                const int new_pos = old_pos + 1;
                if (new_pos + SPS_200 - 1 >= (int)iq.size()) {
                    valid[si] = false;
                    continue;
                }

                const IQ oldx = iq[old_pos];
                const IQ newx = iq[new_pos + SPS_200 - 1];

                const float ti = ci[si] -
                                 (oldx.i * 1.0f) +
                                 (newx.i * ec[si] - newx.q * es[si]);
                const float tq = cq[si] -
                                 (oldx.q * 1.0f) +
                                 (newx.i * es[si] + newx.q * ec[si]);
                ci[si] = ti * pc[si] - tq * ps[si];
                cq[si] = ti * ps[si] + tq * pc[si];
                pos[si] = new_pos;
            }
        }

        best_grid_scores[fi] = best_score_f;
        if (best_score_f > r.score) {
            r.score = best_score_f;
            const int best_delay = first + best_delay_idx_f;
            r.delay_s = window_start_s + (float)best_delay / (float)OUTPUT_FS;
            r.integer_delay_s = r.delay_s;
            r.freq_hz = abs_freq;
            r.freq_grid_hz = abs_freq;
        }
    }

    int best_fi = (int)lroundf((r.freq_grid_hz - coarse_freq_hz +
                                FREQ_RADIUS) / FREQ_STEP);
    best_fi = std::max(1, std::min(NFREQ - 2, best_fi));

    {
        const float ym = best_grid_scores[best_fi - 1];
        const float y0 = best_grid_scores[best_fi];
        const float yp = best_grid_scores[best_fi + 1];
        const float den = ym - 2.0f * y0 + yp;
        float frac = 0.0f;
        if (fabsf(den) > 1e-20f)
            frac = 0.5f * (ym - yp) / den;
        frac = std::max(-0.5f, std::min(0.5f, frac));
        r.freq_hz += frac * FREQ_STEP;
    }

    r.joint_us = bench_now_us() - t0;
    return r;
}

// ============================================================
// Final timing refinement at 200 Hz
// ============================================================
static FineEstimate final_time_200(const std::vector<IQ>& iq,
                                   int iq_first_input,
                                   const uint8_t tones[NTONES],
                                   float delay0_s,
                                   float ddc_lo_hz,
                                   float final_freq_hz)
{
    constexpr int RADIUS = 2;
    constexpr int COUNT = 2 * RADIUS + 1;

    const uint64_t t0 = bench_now_us();
    FineEstimate r;
    r.delay_s = delay0_s;
    r.freq_hz = final_freq_hz;

    const float window_start_s = (float)iq_first_input / (float)INPUT_FS;
    const int center = (int)lroundf((delay0_s - window_start_s) * OUTPUT_FS);
    const int first = std::max(0, center - RADIUS);
    const int last = std::min((int)iq.size() - SPS_200 * NTONES,
                              center + RADIUS);
    if (last < first) {
        r.final_us = bench_now_us() - t0;
        return r;
    }

    const int nscan = std::min(COUNT, last - first + 1);

    static float osc_i[8][SPS_200];
    static float osc_q[8][SPS_200];

    for (int tone = 0; tone < 8; ++tone) {
        const float f = (final_freq_hz - ddc_lo_hz) +
                        TONE_SPACING * (float)tone;
        const float w = TWO_PI * f / (float)OUTPUT_FS;
        const float cd = cosf(w);
        const float sd = sinf(w);
        float c = 1.0f, s = 0.0f;
        for (int n = 0; n < SPS_200; ++n) {
            osc_i[tone][n] = c;
            osc_q[tone][n] = s;
            const float nc = c * cd - s * sd;
            const float ns = s * cd + c * sd;
            c = nc;
            s = ns;
        }
    }

    float score[COUNT] = {};
    bool valid[COUNT] = {};

    for (int d = 0; d < nscan; ++d) {
        const int start = first + d;
        float total = 0.0f;
        int used = 0;

        for (int k = 0; k < NTONES; ++k) {
            const int p0 = start + k * SPS_200;
            if (p0 < 0 || p0 + SPS_200 > (int)iq.size()) continue;

            const int tone = tones[k] & 7;
            float ci = 0.0f;
            float cq = 0.0f;
            for (int n = 0; n < SPS_200; ++n) {
                const float x = iq[p0 + n].i;
                const float y = iq[p0 + n].q;
                const float c = osc_i[tone][n];
                const float s = osc_q[tone][n];
                ci += x * c + y * s;
                cq += y * c - x * s;
            }
            total += ci * ci + cq * cq;
            ++used;
        }

        if (used) {
            score[d] = total / (float)used;
            valid[d] = true;
        }
    }

    float best = -1.0f;
    int bi = 0;
    for (int d = 0; d < nscan; ++d) {
        if (valid[d] && score[d] > best) {
            best = score[d];
            bi = d;
        }
    }

    float frac = 0.0f;
    if (bi > 0 && bi + 1 < nscan) {
        const float ym = score[bi - 1];
        const float y0 = score[bi];
        const float yp = score[bi + 1];
        const float den = ym - 2.0f * y0 + yp;
        if (fabsf(den) > 1e-20f)
            frac = 0.5f * (ym - yp) / den;
        frac = std::max(-0.5f, std::min(0.5f, frac));
    }

    r.score = best;
    r.integer_delay_s = window_start_s +
                         (float)(first + bi) / (float)OUTPUT_FS;
    r.delay_fraction_samples = frac;
    r.delay_s = window_start_s +
                ((float)(first + bi) + frac) / (float)OUTPUT_FS;
    r.final_us = bench_now_us() - t0;
    return r;
}

// ============================================================
// Blind Costas scan for another FT8 signal in the same 200-Hz IQ band.
// It uses only the three known 7-symbol Costas blocks; payload tones are
// deliberately not required, because they are unknown before decoding.
// ============================================================
static void initialize_costas_correlations(
    const std::vector<IQ>& iq,
    const float osc_i[7][SPS_200],
    const float osc_q[7][SPS_200],
    float ci[9], float cq[9], int pos[9])
{
    for (int index = 0; index < 9; ++index) {
        const int p0 = SYNC_SYMBOLS[index] * SPS_200;
        float sum_i = 0.0f;
        float sum_q = 0.0f;
        const int tone = COARSE_COSTAS_TONES[index];

        for (int n = 0; n < SPS_200; ++n) {
            const IQ z = iq[p0 + n];
            const float c = osc_i[tone][n];
            const float s = osc_q[tone][n];
            sum_i += z.i * c + z.q * s;
            sum_q += z.q * c - z.i * s;
        }
        ci[index] = sum_i;
        cq[index] = sum_q;
        pos[index] = p0;
    }
}

static float full_costas_score_200(const std::vector<IQ>& iq, int start,
                                   float freq_hz, float ddc_lo_hz)
{
    static constexpr int COSTAS_BLOCKS[3] = { 0, 36, 72 };
    float total = 0.0f;
    for (int block : COSTAS_BLOCKS) {
        for (int symbol = 0; symbol < 7; ++symbol) {
            const float w = TWO_PI * ((freq_hz - ddc_lo_hz) +
                TONE_SPACING * COSTAS_TONES[symbol]) / (float)OUTPUT_FS;
            const float cd = cosf(w), sd = sinf(w);
            float c = 1.0f, s = 0.0f, ci = 0.0f, cq = 0.0f;
            const int p0 = start + (block + symbol) * SPS_200;
            for (int n = 0; n < SPS_200; ++n) {
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

static std::vector<CostasCandidate> find_costas_candidates_200(
    const std::vector<IQ>& iq,
    int iq_first_input,
    float ddc_lo_hz,
    float min_freq_hz,
    float max_freq_hz)
{
    // Coarse full-slot acquisition: 1 Hz x 8 ms grid.  Fine synchronization
    // is performed only after this stage identifies a Costas peak.
    constexpr float FREQ_STEP_HZ = 1.0f;
    constexpr int TIME_STEP_SAMPLES = 4;
    constexpr int MAX_CANDIDATES = 8;
    constexpr float MIN_FREQ_SEPARATION_HZ = 3.0f;
    constexpr float MIN_TIME_SEPARATION_S = 0.080f;

    std::vector<CostasCandidate> best;
    const int last_start = (int)iq.size() - NTONES * SPS_200;
    const float iq_start_s = (float)iq_first_input / (float)INPUT_FS;
    if (last_start < 0) return best;

    for (float freq = min_freq_hz; freq <= max_freq_hz + 0.001f;
         freq += FREQ_STEP_HZ) {
        // The oscillator templates depend on frequency, but not on time.
        // Build all seven once per frequency instead of once per candidate.
        float osc_i[7][SPS_200];
        float osc_q[7][SPS_200];
        for (int tone = 0; tone < 7; ++tone) {
            const float w = TWO_PI * ((freq - ddc_lo_hz) +
                                      TONE_SPACING * (float)tone) /
                            (float)OUTPUT_FS;
            const float cd = cosf(w);
            const float sd = sinf(w);
            float c = 1.0f;
            float s = 0.0f;
            for (int n = 0; n < SPS_200; ++n) {
                osc_i[tone][n] = c;
                osc_q[tone][n] = s;
                const float nc = c * cd - s * sd;
                const float ns = s * cd + c * sd;
                c = nc;
                s = ns;
            }
        }

        // One 40-sample correlation per Costas symbol initializes the scan.
        // Subsequent timing hypotheses slide each correlation by two samples.
        float ci[9];
        float cq[9];
        int pos[9];
        float pc[9];
        float ps[9];
        float ec[9];
        float es[9];
        initialize_costas_correlations(iq, osc_i, osc_q, ci, cq, pos);

        for (int index = 0; index < 9; ++index) {
            const float w = TWO_PI * ((freq - ddc_lo_hz) +
                TONE_SPACING * COARSE_COSTAS_TONES[index]) / (float)OUTPUT_FS;
            pc[index] = cosf(w);
            ps[index] = sinf(w);
            ec[index] = cosf(-w * SPS_200);
            es[index] = sinf(-w * SPS_200);
        }

        for (int start = 0; start <= last_start; start += TIME_STEP_SAMPLES) {
            CostasCandidate candidate;
            candidate.delay_s = iq_start_s + (float)start / (float)OUTPUT_FS;
            candidate.freq_hz = freq;
            float total = 0.0f;
            for (int i = 0; i < 9; ++i)
                total += ci[i] * ci[i] + cq[i] * cq[i];
            candidate.score = total / 9.0f;

            int overlap_index = -1;
            for (size_t i = 0; i < best.size(); ++i) {
                const CostasCandidate& old = best[i];
                if (fabsf(candidate.freq_hz - old.freq_hz) <
                        MIN_FREQ_SEPARATION_HZ &&
                    fabsf(candidate.delay_s - old.delay_s) <
                        MIN_TIME_SEPARATION_S) {
                    overlap_index = (int)i;
                    break;
                }
            }
            bool changed = false;
            if (overlap_index >= 0) {
                if (candidate.score > best[overlap_index].score) {
                    best[overlap_index] = candidate;
                    changed = true;
                }
            } else if ((int)best.size() < MAX_CANDIDATES) {
                best.push_back(candidate);
                changed = true;
            } else if (candidate.score > best.back().score) {
                best.back() = candidate;
                changed = true;
            }
            if (changed) {
                std::sort(best.begin(), best.end(),
                          [](const CostasCandidate& a, const CostasCandidate& b) {
                              return a.score > b.score;
                          });
            }

            if (start + TIME_STEP_SAMPLES > last_start) continue;
            for (int step = 0; step < TIME_STEP_SAMPLES; ++step) {
                for (int i = 0; i < 9; ++i) {
                    const IQ oldx = iq[pos[i]];
                    const IQ newx = iq[pos[i] + SPS_200];
                    const float ti = ci[i] - oldx.i +
                        newx.i * ec[i] - newx.q * es[i];
                    const float tq = cq[i] - oldx.q +
                        newx.i * es[i] + newx.q * ec[i];
                    ci[i] = ti * pc[i] - tq * ps[i];
                    cq[i] = ti * ps[i] + tq * pc[i];
                    ++pos[i];
                }
            }
        }
    }
    for (CostasCandidate& candidate : best) {
        const int start = (int)lroundf((candidate.delay_s - iq_start_s) *
                                       OUTPUT_FS);
        candidate.score = full_costas_score_200(iq, start,
                                                candidate.freq_hz, ddc_lo_hz);
    }
    std::sort(best.begin(), best.end(),
              [](const CostasCandidate& a, const CostasCandidate& b) {
                  return a.score > b.score;
              });
    return best;
}

static void print_costas_candidates(const char* title,
                                    const std::vector<CostasCandidate>& candidates)
{
    printf("%s\n", title);
    if (candidates.empty()) {
        printf("  none in the available IQ window\n");
        return;
    }
    for (size_t i = 0; i < candidates.size(); ++i) {
        const CostasCandidate& c = candidates[i];
        printf("  %zu: delay %.6f s, tone-0 %.2f Hz, Costas score %.6g\n",
               i + 1, c.delay_s, c.freq_hz, c.score);
    }
}


// ============================================================
// IQ subtraction based on the user's proven subtract.cpp logic.
//
// The original subtract.cpp operates on a REAL waveform at 12 kHz.
// Here the DDC has already produced COMPLEX IQ at 200 Hz, so we use
// the same per-symbol amplitude/phase estimation and the same FT8
// transition law, but synthesize the complex baseband waveform directly.
// ============================================================
struct IQModelParams {
    float amp = 0.0f;
    float phase = 0.0f;
};

static IQ estimate_symbol_coefficient_200(
    const std::vector<IQ>& iq,
    int start,
    int tone,
    float final_freq_hz,
    float ddc_lo_hz)
{
    IQ out{0.0f, 0.0f};
    if (start < 0 || start + SPS_200 > (int)iq.size())
        return out;

    const int block = SPS_200;
    const int first = block / 10;
    const int last  = block - block / 10;
    const int corr_len = last - first;

    const float f = (final_freq_hz - ddc_lo_hz) +
                    TONE_SPACING * (float)(tone & 7);
    const float w = TWO_PI * f / (float)OUTPUT_FS;
    const float cd = cosf(w);
    const float sd = sinf(w);

    // Start at the first sample actually used in the correlation.
    float c = cosf(w * (float)first);
    float s = sinf(w * (float)first);

    float sum_re = 0.0f;
    float sum_im = 0.0f;

    for (int n = first; n < last; ++n) {
        const IQ z = iq[start + n];
        sum_re += z.i * c + z.q * s;
        sum_im += z.q * c - z.i * s;

        const float nc = c * cd - s * sd;
        const float ns = s * cd + c * sd;
        c = nc;
        s = ns;
    }

    // Complex IQ: unlike the original real waveform, there is NO factor 2.
    const float inv = 1.0f / (float)corr_len;
    out.i = sum_re * inv;
    out.q = sum_im * inv;
    return out;
}

static inline void iq_add_rotated(
    IQ& dst, const IQ& a, float c, float s, float sign)
{
    const float ri = a.i * c - a.q * s;
    const float rq = a.i * s + a.q * c;
    dst.i += sign * ri;
    dst.q += sign * rq;
}

// Synthesize one complete 79-symbol FT8 message directly in complex IQ.
// This is the user's synthesize_float() transition law translated from
// a real cosine waveform to an analytic/complex waveform.
static void synthesize_iq_200(
    std::vector<IQ>& dst,
    int start0,
    const uint8_t* tones,
    const IQModelParams* p,
    float hz0,
    float sign)
{
    const int block = SPS_200;
    const int sample_rate = OUTPUT_FS;
    int ramp = (int)lroundf((float)block * 0.11f);
    if (ramp < 1) ramp = 1;

    for (int si = 0; si < NTONES; ++si) {
        const float amp = p[si].amp;
        const float phase = p[si].phase;
        const float freq = hz0 + TONE_SPACING * (float)(tones[si] & 7);
        float dtheta = TWO_PI * freq / (float)sample_rate;

        // -----------------------------------------------------
        // Steady part: amp * exp(j*(phase + n*dtheta))
        // -----------------------------------------------------
        const float theta_start = phase + (float)ramp * dtheta;
        float c = cosf(theta_start);
        float s = sinf(theta_start);
        const float cd = cosf(dtheta);
        const float sd = sinf(dtheta);

        for (int jj = ramp; jj < block - ramp; ++jj) {
            const int idx = start0 + si * block + jj;
            if (idx >= 0 && idx < (int)dst.size()) {
                dst[idx].i += sign * amp * c;
                dst[idx].q += sign * amp * s;
            }
            const float nc = c * cd - s * sd;
            const float ns = s * cd + c * sd;
            c = nc;
            s = ns;
        }

        // -----------------------------------------------------
        // Transition, same phase/frequency law as subtract.cpp
        // -----------------------------------------------------
        float theta = phase + (float)(block - ramp) * dtheta;
        float freq1;
        float phase1;

        if (si == 78) {
            freq1 = freq;
            phase1 = phase;
        } else {
            freq1 = hz0 + TONE_SPACING * (float)(tones[si + 1] & 7);
            phase1 = p[si + 1].phase;
        }

        const float dtheta1 = TWO_PI * freq1 / (float)sample_rate;
        const float inc = (dtheta1 - dtheta) /
                          (2.0f * (float)ramp);

        const float actual =
            theta +
            dtheta * 2.0f * (float)ramp +
            inc * 4.0f * (float)ramp * (float)ramp / 2.0f;

        float target = phase1 + dtheta1 * (float)ramp;
        float adj = remainderf(target - actual, TWO_PI);
        const float adj_step = adj / (2.0f * (float)ramp);

        int end = block + ramp;
        if (si == 78) end = block;
        const int transition_samples = end - (block - ramp);

        float delta = dtheta + adj_step;
        float ct = cosf(theta);
        float st = sinf(theta);
        float rc = cosf(delta);
        float rs = sinf(delta);
        const float ric = cosf(inc);
        const float ris = sinf(inc);

        for (int jj = 0; jj < transition_samples; ++jj) {
            const int idx = start0 + si * block + (block - ramp) + jj;
            if (idx >= 0 && idx < (int)dst.size()) {
                float gain = amp;
                if (si == 78)
                    gain *= 1.0f - (float)jj / (float)ramp;
                dst[idx].i += sign * gain * ct;
                dst[idx].q += sign * gain * st;
            }

            const float nc = ct * rc - st * rs;
            const float ns = st * rc + ct * rs;
            ct = nc;
            st = ns;

            const float nrc = rc * ric - rs * ris;
            const float nrs = rs * ric + rc * ris;
            rc = nrc;
            rs = nrs;
        }
    }
}

static bool subtract_ft8_message_200(
    std::vector<IQ>& iq,
    int iq_first_input,
    const uint8_t tones[NTONES],
    float delay_s,
    float final_freq_hz,
    float ddc_lo_hz,
    float* rms_before,
    float* rms_after,
    float* rms_model,
    float* residual_projection_db)
{
    const float window_start_s = (float)iq_first_input / (float)INPUT_FS;
    const int start0 = (int)lroundf((delay_s - window_start_s) * OUTPUT_FS);
    if (start0 < 0 || start0 + NTONES * SPS_200 > (int)iq.size())
        return false;

    // Estimate amp/phase independently for the 79 symbols, exactly as the
    // supplied subtract.cpp does for the real waveform.
    static IQModelParams model[NTONES];
    for (int k = 0; k < NTONES; ++k) {
        const IQ a = estimate_symbol_coefficient_200(
            iq, start0 + k * SPS_200, tones[k], final_freq_hz, ddc_lo_hz);
        model[k].amp = sqrtf(a.i * a.i + a.q * a.q);
        model[k].phase = atan2f(a.q, a.i);
    }

    // Build the full complex model in a static buffer to avoid loopTask stack.
    static std::vector<IQ> ref;
    ref.assign(iq.size(), IQ{0.0f, 0.0f});
    synthesize_iq_200(ref, start0, tones, model,
                      final_freq_hz - ddc_lo_hz, +1.0f);

    double e_before = 0.0;
    double e_model = 0.0;
    double e_after = 0.0;
    double xs_re = 0.0;
    double xs_im = 0.0;
    double ss = 0.0;
    uint64_t count = 0;

    const int first = std::max(0, start0);
    const int last = std::min((int)iq.size(), start0 + NTONES * SPS_200);

    for (int n = first; n < last; ++n) {
        const IQ x = iq[n];
        const IQ m = ref[n];

        const double mr = m.i;
        const double mq = m.q;
        const double xr = x.i;
        const double xq = x.q;

        // <x,m> for complex vectors: sum conj(m)*x
        xs_re += mr * xr + mq * xq;
        xs_im += mr * xq - mq * xr;
        ss += mr * mr + mq * mq;
        e_before += xr * xr + xq * xq;
        e_model += mr * mr + mq * mq;

        const IQ r{ x.i - m.i, x.q - m.q };
        e_after += (double)r.i * r.i + (double)r.q * r.q;
        iq[n] = r;
        ++count;
    }

    if (!count || ss <= 0.0) return false;

    if (rms_before) *rms_before = sqrtf((float)(e_before / (double)count));
    if (rms_model) *rms_model = sqrtf((float)(e_model / (double)count));
    if (rms_after) *rms_after = sqrtf((float)(e_after / (double)count));

    // <x - m, m> = <x, m> - <m, m>.  This reports how much of the
    // subtracted message remains in the residual, rather than its strength
    // before subtraction.
    const double residual_re = xs_re - ss;
    const double residual_im = xs_im;
    const double proj_mag = sqrt(residual_re * residual_re +
                                 residual_im * residual_im);
    const double alpha_after = proj_mag / ss;
    if (residual_projection_db) {
        *residual_projection_db =
            (alpha_after > 0.0) ? (float)(20.0 * log10(alpha_after)) : -INFINITY;
    }
    return true;
}

// ============================================================
// Save the 200-Hz complex residual so the real FT8 decoder can be
// connected in the next step.
// Format: float32 I,Q interleaved, native little-endian.
// ============================================================
static bool save_residual_200(const std::vector<IQ>& iq, const char* path)
{
#if defined(ARDUINO)
    File f = LittleFS.open(path, "w");
    if (!f) return false;
    for (const IQ& z : iq) {
        if (f.write((const uint8_t*)&z.i, sizeof(float)) != sizeof(float) ||
            f.write((const uint8_t*)&z.q, sizeof(float)) != sizeof(float)) {
            f.close();
            return false;
        }
    }
    f.close();
    return true;
#else
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    for (const IQ& z : iq) {
        f.write((const char*)&z.i, sizeof(float));
        f.write((const char*)&z.q, sizeof(float));
    }
    return (bool)f;
#endif
}

// ============================================================
// WAV loader - ESP32 / Linux
// ============================================================
#if defined(ARDUINO)
static float* load_wav(const char* path,
                       int* out_num_samples,
                       int* out_num_channels,
                       int* out_sample_rate)
{
    File f = LittleFS.open(path, "r");
    if (!f) return nullptr;

    uint8_t riff[12];
    if (f.read(riff, 12) != 12 ||
        memcmp(riff, "RIFF", 4) != 0 ||
        memcmp(riff + 8, "WAVE", 4) != 0) {
        f.close();
        return nullptr;
    }

    uint16_t audio_format = 0, num_channels = 0, bits_per_sample = 0;
    uint32_t sample_rate = 0, data_pos = 0, data_size = 0;
    bool have_fmt = false, have_data = false;

    while (f.available()) {
        uint8_t h[8];
        if (f.read(h, 8) != 8) break;
        const uint32_t chunk_size = h[4] | ((uint32_t)h[5] << 8) |
                                    ((uint32_t)h[6] << 16) | ((uint32_t)h[7] << 24);
        const uint32_t chunk_start = f.position();

        if (memcmp(h, "fmt ", 4) == 0) {
            if (chunk_size < 16) { f.close(); return nullptr; }
            uint8_t fmt[16];
            if (f.read(fmt, 16) != 16) { f.close(); return nullptr; }
            audio_format = fmt[0] | ((uint16_t)fmt[1] << 8);
            num_channels = fmt[2] | ((uint16_t)fmt[3] << 8);
            sample_rate = fmt[4] | ((uint32_t)fmt[5] << 8) |
                          ((uint32_t)fmt[6] << 16) | ((uint32_t)fmt[7] << 24);
            bits_per_sample = fmt[14] | ((uint16_t)fmt[15] << 8);
            have_fmt = true;
        } else if (memcmp(h, "data", 4) == 0) {
            data_pos = f.position();
            data_size = chunk_size;
            have_data = true;
        }

        f.seek(chunk_start + chunk_size);
        if (chunk_size & 1) f.seek(f.position() + 1);
        if (have_fmt && have_data) break;
    }

    if (!have_fmt || !have_data || audio_format != 1 ||
        bits_per_sample != 16 || num_channels < 1) {
        f.close();
        return nullptr;
    }

    const int frame_bytes = 2 * num_channels;
    const int num_samples = data_size / frame_bytes;
    float* samples = (float*)malloc((size_t)num_samples * sizeof(float));
    if (!samples) { f.close(); return nullptr; }

    f.seek(data_pos);
    for (int n = 0; n < num_samples; ++n) {
        uint8_t b[2];
        if (f.read(b, 2) != 2) { free(samples); f.close(); return nullptr; }
        const int16_t v = (int16_t)(b[0] | ((uint16_t)b[1] << 8));
        samples[n] = (float)v / 32768.0f;
        for (int ch = 1; ch < num_channels; ++ch)
            f.seek(f.position() + 2);
    }

    f.close();
    *out_num_samples = num_samples;
    *out_num_channels = (int)num_channels;
    *out_sample_rate = (int)sample_rate;
    return samples;
}
#else
static uint16_t rd_u16(const uint8_t* p)
{
    return (uint16_t)p[0] | ((uint16_t)p[1] << 8);
}
static uint32_t rd_u32(const uint8_t* p)
{
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static float* load_wav(const char* path,
                       int* out_num_samples,
                       int* out_num_channels,
                       int* out_sample_rate)
{
    std::ifstream f(path, std::ios::binary);
    if (!f && path[0] == '/') {
        f.clear();
        f.open(path + 1, std::ios::binary);
    }
    if (!f) return nullptr;

    uint8_t riff[12];
    if (!f.read((char*)riff, 12)) return nullptr;
    if (memcmp(riff, "RIFF", 4) != 0 || memcmp(riff + 8, "WAVE", 4) != 0)
        return nullptr;

    uint16_t audio_format = 0, num_channels = 0;
    uint16_t bits_per_sample = 0, block_align = 0;
    uint32_t sample_rate = 0, data_size = 0;
    std::streamoff data_pos = 0;
    bool have_fmt = false, have_data = false;

    while (f) {
        uint8_t h[8];
        if (!f.read((char*)h, 8)) break;
        const uint32_t chunk_size = rd_u32(h + 4);
        const std::streamoff chunk_start = f.tellg();

        if (memcmp(h, "fmt ", 4) == 0) {
            if (chunk_size < 16) return nullptr;
            std::vector<uint8_t> fmt(chunk_size);
            if (!f.read((char*)fmt.data(), chunk_size)) return nullptr;
            audio_format = rd_u16(fmt.data());
            num_channels = rd_u16(fmt.data() + 2);
            sample_rate = rd_u32(fmt.data() + 4);
            block_align = rd_u16(fmt.data() + 12);
            bits_per_sample = rd_u16(fmt.data() + 14);
            have_fmt = true;
        } else if (memcmp(h, "data", 4) == 0) {
            data_pos = chunk_start;
            data_size = chunk_size;
            have_data = true;
        } else {
            f.seekg((std::streamoff)chunk_size, std::ios::cur);
        }

        if (chunk_size & 1) f.seekg(1, std::ios::cur);
        if (have_fmt && have_data) break;
    }

    if (!have_fmt || !have_data || audio_format != 1 ||
        bits_per_sample != 16 || num_channels < 1) return nullptr;

    const int num_samples = data_size / block_align;
    float* samples = (float*)malloc((size_t)num_samples * sizeof(float));
    if (!samples) return nullptr;

    f.clear();
    f.seekg(data_pos);
    for (int n = 0; n < num_samples; ++n) {
        uint8_t b[2];
        if (!f.read((char*)b, 2)) { free(samples); return nullptr; }
        const int16_t v = (int16_t)(b[0] | ((uint16_t)b[1] << 8));
        samples[n] = (float)v / 32768.0f;
        const std::streamoff skip = (std::streamoff)block_align - 2;
        if (skip > 0) f.seekg(skip, std::ios::cur);
    }

    *out_num_samples = num_samples;
    *out_num_channels = (int)num_channels;
    *out_sample_rate = (int)sample_rate;
    return samples;
}
#endif

// ============================================================
// V10 subtraction test
// ============================================================
static void run_all_tests()
{
    printf("\r\n\r\n");
    printf("############################################################\n");
    printf(" FT8 V10.1 - 200 Hz IQ SUBTRACTION USING subtract.cpp MODEL\n");
    printf("############################################################\n");
#if defined(ESP32)
    printf("Platform         : ESP32-S3 / Arduino\n");
    printf("CPU frequency    : %u MHz\n", ESP.getCpuFreqMHz());
#else
    printf("Platform         : Linux / Ubuntu\n");
#endif
    printf("Input sample rate: %d Hz\n", INPUT_FS);
    printf("Output sample rate: %d Hz\n", OUTPUT_FS);
    printf("DDC                : ONE pass at %.1f Hz\n", DDC_LO_TEST);
    printf("Purpose            : synchronize -> reconstruct A -> subtract A\n");
#if defined(ARDUINO)
    const char* residual_path = "/residual_200hz_iq.bin";
#else
    const char* residual_path = "residual_200hz_iq.bin";
#endif
    printf("Residual output    : %s\n", residual_path);

#if defined(ESP32)
    printf("Free heap        : %u bytes\n", ESP.getFreeHeap());
    printf("Free PSRAM       : %u bytes\n", ESP.getFreePsram());
    printf("\r\nMounting LittleFS...\n");
    if (!LittleFS.begin(true)) {
        printf("ERROR: LittleFS.begin() failed\n");
        return;
    }
    printf("LittleFS mounted.\n");
#endif

#if defined(ARDUINO)
    const char* filename = "/test_real_ft8.wav";
#else
    const char* filename = "test_real_ft8.wav";
#endif
    int num_samples = 0, num_channels = 0, sample_rate = 0;
    float* samples = load_wav(filename, &num_samples, &num_channels, &sample_rate);
    if (!samples) {
        printf("ERROR: WAV load failed\n");
        return;
    }
    if (sample_rate != INPUT_FS || num_channels != 1) {
        printf("ERROR: WAV format mismatch\n");
        free(samples);
        return;
    }

    uint8_t tones[NTONES];
    if (!extract_ft8_tones(samples, num_samples, 2.0f, TEST_FREQ, tones)) {
        printf("ERROR: cannot extract known FT8 tones\n");
        free(samples);
        return;
    }

    // The waterfall supplies one coarse hypothesis. In this controlled
    // test use the corner -50 ms / -2 Hz, exactly as in V9.x.
    const float coarse_delay = 1.950f;
    const float coarse_freq  = 1498.0f;
    // Retain the complete FT8 slot.  This permits blind candidates from
    // T = 0 through the last complete 12.64-second message in the WAV.
    const float WIN_START = 0.0f;
    const float WIN_END = (float)num_samples / (float)INPUT_FS;
    int iq_first_input = 0;

    const uint64_t d0 = bench_now_us();
    std::vector<IQ> iq200 = ddc_decimate_200_window(
        samples, num_samples, DDC_LO_TEST,
        WIN_START, WIN_END, &iq_first_input);
    const uint64_t ddc_us = bench_now_us() - d0;

    printf("\n------------------------------------------------------------\n");
    printf("DDC\n");
    printf("SEARCH       : T 0..2360 ms, F0 1425..1575 Hz\n");
    printf("DDC LO       : %.1f Hz\n", DDC_LO_TEST);
    printf("IQ samples   : %zu\n", iq200.size());
    printf("DDC TIME     : %.3f ms\n", (double)ddc_us / 1000.0);

    const uint64_t coarse_scan_start = bench_now_us();
    const std::vector<CostasCandidate> before_candidates =
        find_costas_candidates_200(iq200, iq_first_input, DDC_LO_TEST,
                                   1425.0f, 1575.0f);
    const uint64_t coarse_scan_us = bench_now_us() - coarse_scan_start;
    printf("BLIND SCAN TIME: %.3f ms\n", (double)coarse_scan_us / 1000.0);
    print_costas_candidates("COSTAS CANDIDATES BEFORE SUBTRACTION",
                            before_candidates);
    if (before_candidates.empty()) {
        printf("ERROR: no primary FT8 Costas candidate found\n");
        free(samples);
        return;
    }

    FineEstimate j = fine_joint_200(
        iq200, iq_first_input, tones,
        coarse_delay, coarse_freq, DDC_LO_TEST);
    FineEstimate f = final_time_200(
        iq200, iq_first_input, tones,
        j.delay_s, DDC_LO_TEST, j.freq_hz);

    printf("\n------------------------------------------------------------\n");
    printf("FINE SYNC\n");
    printf("FREQUENCY   : %.6f Hz\n", f.freq_hz);
    printf("DELAY       : %.6f s\n", f.delay_s);
    printf("DELAY GRID  : %.6f s\n", f.integer_delay_s);
    printf("TIME ERROR  : %+0.3f ms\n", (f.delay_s - 2.000f) * 1000.0f);
    printf("FREQ ERROR  : %+0.6f Hz\n", f.freq_hz - TEST_FREQ);
    printf("SYNC TIME   : %.3f ms\n",
           (double)(j.joint_us + f.final_us) / 1000.0);

    float rms_before = 0.0f, rms_after = 0.0f, rms_model = 0.0f, residual_projection_db = 0.0f;
    const uint64_t s0 = bench_now_us();
    // The blind Costas peak supplies the exact sample-aligned start used by
    // the current 200-Hz synthesizer.  Fractional timing needs a fractional
    // delay model before it can safely be used for subtraction.
    const CostasCandidate& primary = before_candidates.front();
    const bool ok = subtract_ft8_message_200(
        iq200, iq_first_input, tones,
        primary.delay_s, primary.freq_hz,
        DDC_LO_TEST,
        &rms_before, &rms_after, &rms_model, &residual_projection_db);
    const uint64_t sus = bench_now_us() - s0;

    printf("\n------------------------------------------------------------\n");
    printf("SUBTRACTION OF DECODED MESSAGE A\n");
    if (!ok) {
        printf("ERROR: subtraction failed\n");
        free(samples);
        return;
    }
    printf("RMS BEFORE  : %.8f\n", rms_before);
    printf("RMS MODEL   : %.8f\n", rms_model);
    printf("RMS RESIDUAL : %.8f\n", rms_after);
    if (rms_after > 0.0f) {
        const float db = 20.0f * log10f(rms_before / rms_after);
        printf("CANCELLATION (RMS): %.2f dB\n", db);
        printf("RESIDUAL MODEL PROJECTION: %.2f dB\n", residual_projection_db);
    }
    printf("SUBTRACT TIME: %.3f ms\n", (double)sus / 1000.0);

    const uint64_t residual_scan_start = bench_now_us();
    const std::vector<CostasCandidate> after_candidates =
        find_costas_candidates_200(iq200, iq_first_input, DDC_LO_TEST,
                                   1425.0f, 1575.0f);
    const uint64_t residual_scan_us = bench_now_us() - residual_scan_start;
    printf("RESIDUAL SCAN TIME: %.3f ms\n",
           (double)residual_scan_us / 1000.0);
    print_costas_candidates("COSTAS CANDIDATES AFTER SUBTRACTION",
                            after_candidates);
    if (!after_candidates.empty()) {
        const CostasCandidate& weak = after_candidates.front();
        const bool weak_verified =
            fabsf(weak.delay_s - WEAK_TEST_DELAY_S) <= 0.020f &&
            fabsf(weak.freq_hz - WEAK_TEST_FREQ_HZ) <= 2.0f;
        printf("WEAK MESSAGE  : expected %.3f s / %.1f Hz -> %s\n",
               WEAK_TEST_DELAY_S, WEAK_TEST_FREQ_HZ,
               weak_verified ? "VERIFIED" : "NOT VERIFIED");
    }

    const bool saved = save_residual_200(iq200, residual_path);
    printf("RESIDUAL SAVE: %s\n", saved ? "OK" : "FAIL");
    printf("NEXT STEP    : feed residual IQ to the FT8 decoder / waterfall\n");

    printf("\n------------------------------------------------------------\n");
    printf("TOTAL REAL PATH\n");
    printf("DDC          : %.3f ms\n", (double)ddc_us / 1000.0);
    printf("FINE SYNC    : %.3f ms\n",
           (double)(j.joint_us + f.final_us) / 1000.0);
    printf("SUBTRACTION  : %.3f ms\n", (double)sus / 1000.0);
    printf("BLIND SCAN   : %.3f ms\n", (double)coarse_scan_us / 1000.0);
    printf("RESIDUAL SCAN: %.3f ms\n", (double)residual_scan_us / 1000.0);
    printf("TOTAL        : %.3f ms\n", (double)(ddc_us + j.joint_us +
           f.final_us + sus + coarse_scan_us + residual_scan_us) / 1000.0);

    free(samples);
    printf("\n\n############################################################\n");
    printf("V10.1 IQ SUBTRACTION TEST COMPLETE\n");
#if defined(ESP32)
    printf("Final free heap : %u bytes\n", ESP.getFreeHeap());
    printf("Final free PSRAM: %u bytes\n", ESP.getFreePsram());
#endif
    printf("############################################################\n");
}

#if defined(ARDUINO)
void setup()
{
    Serial.begin(115200);
    delay(2000);
    run_all_tests();
}
void loop()
{
    delay(1000);
}
#else
int main()
{
    run_all_tests();
    return 0;
}
#endif
