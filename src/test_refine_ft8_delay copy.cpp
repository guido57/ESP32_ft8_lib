// ============================================================
// test_refine_ft8_delay_original_wav.cpp
//
// BENCHMARK DELLA refine_ft8_delay() ORIGINALE
//
// ESP32-S3 + LittleFS
//
// Test:
//   /test_real_ft8_D1500.wav
//   /test_real_ft8_D2000.wav
//   /test_real_ft8_D2500.wav
//
// Ogni WAV viene:
//   1. caricato da LittleFS
//   2. analizzato per ricavare i 79 FT8 tones
//   3. passato direttamente a refine_ft8_delay()
//   4. liberato
//   5. si passa al file successivo
//
// IMPORTANTE:
//   refine_ft8_delay() NON viene modificata.
//
// ============================================================


#if defined(ARDUINO)
#include <Arduino.h>
#include <LittleFS.h>
#else
#include <chrono>
#include <thread>
#include <fstream>
#include <cstdint>
#include <cstring>
#endif

#include <math.h>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cstring>

#ifndef PI
static constexpr float PI = 3.14159265358979323846f;
#endif

#ifndef TWO_PI
static constexpr float TWO_PI = 6.28318530717958647692f;
#endif

// Portable timing helper.
// Arduino: micros()
// Linux:   std::chrono::steady_clock
static inline uint64_t portable_micros()
{
#if defined(ARDUINO)
    return (uint64_t)micros();
#else
    static const auto t0 = std::chrono::steady_clock::now();
    const auto now = std::chrono::steady_clock::now();
    return (uint64_t)std::chrono::duration_cast<
        std::chrono::microseconds>(now - t0).count();
#endif
}

// ============================================================
// CONFIG
// ============================================================

static constexpr int SAMPLE_RATE = 12000;
static constexpr int NTONES      = 79;
static constexpr int SPS         = 1920;       // 160 ms
static constexpr float TONE_SPACING = 6.25f;

static constexpr float TEST_FREQ = 1500.0f;

static inline int blocksize(int sample_rate)
{
    // FT8: 160 ms
    return (sample_rate * 160) / 1000;
}



// ============================================================
// refine_ft8_delay_v4()
//
// FAST FT8 DELAY ESTIMATOR
//
// Pipeline:
//
//   RF 12 kHz
//      |
//      v
//   complex downconversion
//      |
//      v
//   FIR LPF 200 Hz
//      |
//      v
//   decimate x4
//      |
//      v
//   complex baseband @ 3 kHz
//      |
//      v
//   SLIDING correlation, coarse delay
//      |
//      v
//   fine correlation @ 12 kHz
//      |
//      v
//   1-sample final resolution
//
// Important:
//   - coarse search does NOT recalculate 480 samples for
//     every delay
//   - one initial correlation + sliding update
//   - fine search is only +/- 3 samples
//
// ============================================================

float refine_ft8_delay_v4(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)cand_to_print;

    constexpr int FS       = 12000;
    constexpr int NTONES   = 79;
    constexpr int SPS      = 1920;

    constexpr int DECIM    = 4;
    constexpr int FS_DEC   = FS / DECIM;
    constexpr int SPS_DEC  = SPS / DECIM;

    constexpr int FIR_TAPS = 129;
    constexpr int FIR_DELAY = (FIR_TAPS - 1) / 2;

    constexpr float CUTOFF_HZ = 200.0f;

    
    // ========================================================
    // DELAY RANGE
    //
    // Coarse search is performed in units of 4 RF samples.
    // ========================================================

    const int first_delay =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.05f) * FS));

    const int last_delay =
        std::min(
            num_samples - 1,
            (int)std::lround(
                (delay0 + 0.05f) * FS));

    if (last_delay < first_delay)
        return delay0;

    // ========================================================
    // FIR coefficients
    //
    // Hann-windowed sinc LPF
    // ========================================================

    float fir[FIR_TAPS];

    const float fc =
        CUTOFF_HZ / (float)FS;

    float fir_sum = 0.0f;

    for (int n = 0; n < FIR_TAPS; ++n)
    {
        const float m =
            (float)n -
            (float)FIR_DELAY;

        float h;

        if (fabsf(m) < 1.0e-12f)
        {
            h = 2.0f * fc;
        }
        else
        {
            h =
                sinf(TWO_PI * fc * m) /
                (PI * m);
        }

        const float w =
            0.5f -
            0.5f *
            cosf(
                TWO_PI *
                (float)n /
                (float)(FIR_TAPS - 1));

        fir[n] = h * w;

        fir_sum += fir[n];
    }

    const float inv_sum =
        1.0f / fir_sum;

    for (int n = 0; n < FIR_TAPS; ++n)
        fir[n] *= inv_sum;

    // ========================================================
    // RF processing range
    //
    // FIR output corresponding to first_delay occurs at:
    //
    //     first_delay + FIR_DELAY
    //
    // ========================================================

    const int rf_start =
        std::max(
            0,
            first_delay - (FIR_TAPS - 1));

    const int rf_end =
        std::min(
            num_samples - 1,
            last_delay +
            NTONES * SPS +
            FIR_DELAY);

    if (rf_end < rf_start)
        return delay0;

    const int rf_len =
        rf_end - rf_start + 1;

    // ========================================================
    // Decimated baseband buffer
    //
    // One complex sample every 4 RF samples.
    // ========================================================

    const int dec_len =
        (rf_len + DECIM - 1) / DECIM;

    std::vector<float> bi(dec_len);
    std::vector<float> bq(dec_len);

    // ========================================================
    // Circular FIR history
    // ========================================================

    float hist_i[FIR_TAPS] = {};
    float hist_q[FIR_TAPS] = {};

    int hist_pos = 0;
    int out_count = 0;

    // ========================================================
    // LO
    //
    // RF -> complex baseband:
    //
    // I = x*cos()
    // Q = -x*sin()
    // ========================================================

    const float lo_w =
        TWO_PI *
        freq /
        (float)FS;

    const float lo_cd =
        cosf(lo_w);

    const float lo_sd =
        sinf(lo_w);

    const float initial_phase =
        TWO_PI *
        freq *
        (float)rf_start /
        (float)FS;

    float lo_c =
        cosf(initial_phase);

    float lo_s =
        sinf(initial_phase);

    // ========================================================
    // Generate FIR + decimation
    //
    // The decimation phase is chosen so that the first output
    // corresponds to first_delay after FIR group-delay
    // compensation.
    //
    // ========================================================

    const int first_output_abs =
        first_delay + FIR_DELAY;

    for (int n = 0; n < rf_len; ++n)
    {
        const int abs_n =
            rf_start + n;

        const float x =
            samples[abs_n];

        const float xi =
            x * lo_c;

        const float xq =
            -x * lo_s;

        hist_i[hist_pos] = xi;
        hist_q[hist_pos] = xq;

        // ----------------------------------------------------
        // Only calculate FIR output at the required decimation
        // phase.
        // ----------------------------------------------------

        if (abs_n >= first_output_abs &&
            ((abs_n - first_output_abs) % DECIM) == 0)
        {
            float yi = 0.0f;
            float yq = 0.0f;

            int p = hist_pos;

            for (int k = 0;
                 k < FIR_TAPS;
                 ++k)
            {
                const float h = fir[k];

                yi += hist_i[p] * h;
                yq += hist_q[p] * h;

                --p;

                if (p < 0)
                    p = FIR_TAPS - 1;
            }

            if (out_count < dec_len)
            {
                bi[out_count] = yi;
                bq[out_count] = yq;
                ++out_count;
            }
        }

        ++hist_pos;

        if (hist_pos >= FIR_TAPS)
            hist_pos = 0;

        // ----------------------------------------------------
        // oscillator recurrence
        // ----------------------------------------------------

        const float nc =
            lo_c * lo_cd -
            lo_s * lo_sd;

        const float ns =
            lo_s * lo_cd +
            lo_c * lo_sd;

        lo_c = nc;
        lo_s = ns;

        // Occasionally renormalize oscillator.
        // Prevents long-term amplitude drift.
        if ((n & 0x0FFF) == 0x0FFF)
        {
            const float norm =
                1.0f /
                sqrtf(
                    lo_c * lo_c +
                    lo_s * lo_s);

            lo_c *= norm;
            lo_s *= norm;
        }
    }

    const int actual_dec_len =
        out_count;

    if (actual_dec_len <= 0)
        return delay0;

    // ========================================================
    // COARSE SEARCH
    //
    // Delay step = 4 RF samples = 1 decimated sample.
    //
    // Number of candidates:
    //
    //     ~301
    //
    // For every tone:
    //
    //     initial correlation = 480 samples
    //
    // then each following delay requires only:
    //
    //     remove old sample
    //     rotate
    //     add new sample
    //
    // ========================================================

    const int coarse_first =
        first_delay / DECIM;

    const int coarse_last =
        last_delay / DECIM;

    const int ncoarse =
        coarse_last -
        coarse_first +
        1;

    std::vector<float>
        coarse_score(ncoarse, 0.0f);

    std::vector<int>
        coarse_used(ncoarse, 0);

    for (int k = 0; k < NTONES; ++k)
    {
        const float tone =
            6.25f *
            (float)tones[k];

        const float w =
            TWO_PI *
            tone /
            (float)FS_DEC;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);

        // ----------------------------------------------------
        // First delay
        // ----------------------------------------------------

        const int pos0 =
            k * SPS_DEC;

        if (pos0 >= actual_dec_len)
            continue;

        const int count0 =
            std::min(
                SPS_DEC,
                actual_dec_len - pos0);

        if (count0 <= 0)
            continue;

        float ci = 0.0f;
        float cq = 0.0f;

        float c = 1.0f;
        float s = 0.0f;

        for (int n = 0;
             n < count0;
             ++n)
        {
            const float ai =
                bi[pos0 + n];

            const float aq =
                bq[pos0 + n];

            // Multiply baseband by exp(-j tone*n)
            ci +=
                ai * c +
                aq * s;

            cq +=
                aq * c -
                ai * s;

            const float nc =
                c * cd -
                s * sd;

            const float ns =
                s * cd +
                c * sd;

            c = nc;
            s = ns;
        }

        coarse_score[0] +=
            ci * ci +
            cq * cq;

        coarse_used[0]++;

        // ----------------------------------------------------
        // Sliding through delay candidates
        // ----------------------------------------------------

        for (int d = 1;
             d < ncoarse;
             ++d)
        {
            const int pos =
                d +
                k * SPS_DEC;

            const int new_index =
                pos +
                SPS_DEC -
                1;

            if (pos < 0 ||
                new_index >= actual_dec_len)
                break;

            const float old_i =
                bi[pos - 1];

            const float old_q =
                bq[pos - 1];

            const float new_i =
                bi[new_index];

            const float new_q =
                bq[new_index];

            // Remove oldest sample.
            //
            // The correlation accumulator corresponds to the
            // current phase c/s.
            //
            // old sample enters with phase 0 after the
            // appropriate rotation below.

            const float ri =
                ci -
                old_i;

            const float rq =
                cq -
                old_q;

            // Rotate accumulator by one sample.
            const float rotated_i =
                ri * cd -
                rq * sd;

            const float rotated_q =
                ri * sd +
                rq * cd;

            // End phase of the window.
            //
            // c/s currently correspond to the phase following
            // the last sample of the initial window.
            //
            float end_c = 1.0f;
            float end_s = 0.0f;

            // We need exp(-j tone*(SPS_DEC-1)).
            //
            // Compute it once for this tone.
            //
            // This loop is intentionally outside the expensive
            // delay loop in the optimized implementation below.
            //
            // Replaced immediately after this block.
            (void)end_c;
            (void)end_s;

            // ------------------------------------------------
            // For correctness use the phase of the new sample.
            // ------------------------------------------------

            // At this point c/s have been advanced to the phase
            // after count0 samples. For a complete symbol:
            //
            // phase(new_index) = exp(-j*2*pi*tone*(SPS_DEC-1))
            //
            // Calculate using the known final phase.
            // ------------------------------------------------

            static float dummy = 0.0f;
            (void)dummy;

            // The first window is always 480 samples except
            // the final truncated symbol. For the normal case,
            // calculate end phase analytically.

            float ec;
            float es;

            const float end_phase =
                -w *
                (float)(SPS_DEC - 1);

            ec = cosf(end_phase);
            es = sinf(end_phase);

            ci =
                rotated_i +
                new_i * ec -
                new_q * es;

            cq =
                rotated_q +
                new_i * es +
                new_q * ec;

            coarse_score[d] +=
                ci * ci +
                cq * cq;

            coarse_used[d]++;
        }
    }

    // ========================================================
    // Find coarse maximum
    // ========================================================

    int best_coarse_index = 0;

    float best_coarse_score =
        -1.0f;

    for (int d = 0;
         d < ncoarse;
         ++d)
    {
        if (coarse_used[d] == 0)
            continue;

        const float score =
            coarse_score[d] /
            (float)coarse_used[d];

        if (score > best_coarse_score)
        {
            best_coarse_score =
                score;

            best_coarse_index = d;
        }
    }

    // --------------------------------------------------------
    // Coarse delay in original 12 kHz samples.
    // --------------------------------------------------------

    const int coarse_delay =
        (coarse_first +
         best_coarse_index) *
        DECIM;

    // ========================================================
    // FINE SEARCH
    //
    // Only +/- 3 samples.
    //
    // IMPORTANT:
    // The fine search uses the known-good original RF
    // correlation, preserving exact 1-sample resolution.
    // ========================================================

    constexpr int FINE_RADIUS = 3;

    const int fine_first =
        std::max(
            first_delay,
            coarse_delay -
            FINE_RADIUS);

    const int fine_last =
        std::min(
            last_delay,
            coarse_delay +
            FINE_RADIUS);

    const int nfine =
        fine_last -
        fine_first +
        1;

    std::vector<float>
        fine_score(nfine, 0.0f);

    std::vector<int>
        fine_used(nfine, 0);

    for (int k = 0; k < NTONES; ++k)
    {
        const float tone_freq =
            6.25f *
            (float)tones[k];

        const float w =
            TWO_PI *
            tone_freq /
            (float)FS;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);

        // ----------------------------------------------------
        // Initial fine delay
        // ----------------------------------------------------

        const int start0 =
            fine_first +
            k * SPS;

        if (start0 < 0 ||
            start0 >= num_samples)
            continue;

        const int count =
            std::min(
                SPS,
                num_samples - start0);

        if (count <= 0)
            continue;

        float ci = 0.0f;
        float cq = 0.0f;

        float c = 1.0f;
        float s = 0.0f;

        for (int n = 0;
             n < count;
             ++n)
        {
            const float x =
                samples[start0 + n];

            ci += x * c;
            cq -= x * s;

            const float nc =
                c * cd -
                s * sd;

            const float ns =
                s * cd +
                c * sd;

            c = nc;
            s = ns;
        }

        fine_score[0] +=
            ci * ci +
            cq * cq;

        fine_used[0]++;

        // ----------------------------------------------------
        // Slide only a few samples
        // ----------------------------------------------------

        for (int d = 1;
             d < nfine;
             ++d)
        {
            const int start =
                fine_first +
                d +
                k * SPS;

            const int new_index =
                start +
                count -
                1;

            if (start <= 0 ||
                new_index >= num_samples)
                break;

            const float old_x =
                samples[start - 1];

            const float new_x =
                samples[new_index];

            const float ri =
                ci -
                old_x;

            const float rq =
                cq;

            const float rotated_re =
                ri * cd -
                rq * sd;

            const float rotated_im =
                ri * sd +
                rq * cd;

            // Phase of the new sample.
            const float phase =
                -w *
                (float)(count - 1);

            const float ec =
                cosf(phase);

            const float es =
                sinf(phase);

            ci =
                rotated_re +
                new_x * ec;

            cq =
                rotated_im +
                new_x * es;

            fine_score[d] +=
                ci * ci +
                cq * cq;

            fine_used[d]++;
        }
    }

    // ========================================================
    // FINAL MAXIMUM
    // ========================================================

    int best_delay =
        fine_first;

    float best_score =
        -1.0f;

    for (int d = 0;
         d < nfine;
         ++d)
    {
        if (fine_used[d] == 0)
            continue;

        const float score =
            fine_score[d] /
            (float)fine_used[d];

        if (score > best_score)
        {
            best_score =
                score;

            best_delay =
                fine_first + d;
        }
    }

    return
        best_delay /
        (float)FS;
}
// ============================================================
// refine_ft8_delay_v3()
//
// FAST FT8 DELAY ESTIMATOR
//
// Pipeline:
//
//   RF @ freq
//       |
//       v
//   complex downconversion
//       |
//       v
//   FIR LPF 200 Hz
//       |
//       v
//   decimate x4
//       |
//       v
//   Fs = 3000 Hz
//       |
//       v
//   coarse delay search
//       |
//       v
//   local refinement @ 12 kHz
//       |
//       v
//   final delay, 1-sample resolution
//
// Target:
//   ESP32-S3 <= 15 ms
//
// Reference:
//   original estimator: 0 ... -2 samples on real WAVs
//
// ============================================================

float refine_ft8_delay_v3(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int FS = 12000;
    constexpr int NTONES = 79;
    constexpr int SPS = 1920;

    // --------------------------------------------------------
    // Decimation
    // --------------------------------------------------------

    constexpr int DECIM = 4;
    constexpr int FS_DEC = FS / DECIM;
    constexpr int SPS_DEC = SPS / DECIM;

    // 1920 / 4 = 480
    static_assert(SPS_DEC == 480);

    // --------------------------------------------------------
    // FIR
    //
    // 129 taps
    // Fc = 200 Hz
    // Group delay = 64 samples @ 12 kHz
    // --------------------------------------------------------

    constexpr int FIR_TAPS = 129;
    constexpr int FIR_DELAY = 64;
    constexpr float CUTOFF_HZ = 200.0f;

    // constexpr float TWO_PI =
    //     6.2831853071795864769f;

    // --------------------------------------------------------
    // Coarse search range
    //
    // At 3 kHz:
    //
    // +/- 50 ms = +/- 150 samples
    //
    // Therefore only 301 candidates.
    // --------------------------------------------------------

    const int coarse_first =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.05f) *
                FS_DEC));

    const int coarse_last =
        std::min(
            num_samples / DECIM,
            (int)std::lround(
                (delay0 + 0.05f) *
                FS_DEC));

    const int ncoarse =
        coarse_last -
        coarse_first +
        1;

    if (ncoarse <= 0)
        return delay0;

    // --------------------------------------------------------
    // FIR coefficients
    // --------------------------------------------------------

    float fir[FIR_TAPS];

    const float fc =
        CUTOFF_HZ /
        (float)FS;

    float fir_sum = 0.0f;

    for (int n = 0;
         n < FIR_TAPS;
         ++n)
    {
        const float m =
            (float)n -
            (float)FIR_DELAY;

        float h;

        if (fabsf(m) < 1.0e-12f)
        {
            h = 2.0f * fc;
        }
        else
        {
            h =
                sinf(TWO_PI * fc * m) /
                (float)(M_PI * m);
        }

        const float win =
            0.5f -
            0.5f *
            cosf(
                TWO_PI *
                (float)n /
                (float)(FIR_TAPS - 1));

        fir[n] = h * win;

        fir_sum += fir[n];
    }

    const float inv_sum =
        1.0f / fir_sum;

    for (int n = 0;
         n < FIR_TAPS;
         ++n)
    {
        fir[n] *= inv_sum;
    }

    // --------------------------------------------------------
    // We process only the FT8 message region.
    //
    // Need enough RF samples for:
    //
    //   delay0 +/- 50 ms
    //   79 symbols
    //   FIR history
    // --------------------------------------------------------

    const int first_rf =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.05f) *
                FS) -
            FIR_DELAY);

    const int last_rf =
        std::min(
            num_samples - 1,
            (int)std::lround(
                (delay0 + 0.05f) *
                FS) +
            NTONES * SPS +
            FIR_DELAY);

    if (last_rf <= first_rf)
        return delay0;

    const int region_len =
        last_rf -
        first_rf +
        1;

    // --------------------------------------------------------
    // Decimated complex baseband buffer.
    //
    // We only keep every fourth FIR output.
    // --------------------------------------------------------

    const int dec_len =
        (region_len + DECIM - 1) /
        DECIM;

    std::vector<float> dec_i(dec_len);
    std::vector<float> dec_q(dec_len);

    // --------------------------------------------------------
    // RF -> complex baseband
    //
    // Instead of storing all 12-kHz baseband samples, use a
    // circular FIR history and produce only every 4th output.
    // --------------------------------------------------------

    float hist_i[FIR_TAPS] = {};
    float hist_q[FIR_TAPS] = {};

    int hist_pos = 0;
    int dec_index = 0;

    const float lo_w =
        TWO_PI *
        freq /
        (float)FS;

    const float lo_cd =
        cosf(lo_w);

    const float lo_sd =
        sinf(lo_w);

    const float start_phase =
        TWO_PI *
        freq *
        (float)first_rf /
        (float)FS;

    float lo_c =
        cosf(start_phase);

    float lo_s =
        sinf(start_phase);

    for (int n = 0;
         n < region_len;
         ++n)
    {
        const float x =
            samples[first_rf + n];

        const float xi =
            x * lo_c;

        const float xq =
            -x * lo_s;

        // ----------------------------------------------------
        // Insert into circular FIR buffer
        // ----------------------------------------------------

        hist_i[hist_pos] = xi;
        hist_q[hist_pos] = xq;

        // ----------------------------------------------------
        // FIR only every DECIM samples.
        //
        // We still need the input samples continuously for
        // correct filtering, but calculate an output only
        // every 4 samples.
        // ----------------------------------------------------

        if ((n % DECIM) == 0)
        {
            float yi = 0.0f;
            float yq = 0.0f;

            int p = hist_pos;

            for (int k = 0;
                 k < FIR_TAPS;
                 ++k)
            {
                yi +=
                    hist_i[p] *
                    fir[k];

                yq +=
                    hist_q[p] *
                    fir[k];

                --p;

                if (p < 0)
                    p = FIR_TAPS - 1;
            }

            dec_i[dec_index] = yi;
            dec_q[dec_index] = yq;

            ++dec_index;
        }

        ++hist_pos;

        if (hist_pos >= FIR_TAPS)
            hist_pos = 0;

        // ----------------------------------------------------
        // Oscillator update
        // ----------------------------------------------------

        const float nc =
            lo_c * lo_cd -
            lo_s * lo_sd;

        const float ns =
            lo_s * lo_cd +
            lo_c * lo_sd;

        lo_c = nc;
        lo_s = ns;
    }

    const int actual_dec_len =
        dec_index;

    if (actual_dec_len <= 0)
        return delay0;

    // --------------------------------------------------------
    // COARSE SEARCH @ 3 kHz
    // --------------------------------------------------------

    std::vector<float>
        coarse_score(
            ncoarse,
            0.0f);

    std::vector<int>
        coarse_used(
            ncoarse,
            0);

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        const float tone =
            6.25f *
            (float)tones[k];

        const float w =
            TWO_PI *
            tone /
            (float)FS_DEC;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);

        for (int d = 0;
             d < ncoarse;
             ++d)
        {
            const int delay_dec =
                coarse_first + d;

            // Convert original RF delay to the decimated
            // buffer coordinate.
            //
            // FIR output has 64-sample RF delay.
            //

            const int rf_start =
                delay_dec * DECIM +
                k * SPS;

            const int pos =
                (rf_start -
                 first_rf +
                 FIR_DELAY) /
                DECIM;

            if (pos < 0 ||
                pos >= actual_dec_len)
                continue;

            const int available =
                actual_dec_len - pos;

            const int count =
                std::min(
                    SPS_DEC,
                    available);

            if (count <= 0)
                continue;

            float ci = 0.0f;
            float cq = 0.0f;

            float c = 1.0f;
            float s = 0.0f;

            for (int n = 0;
                 n < count;
                 ++n)
            {
                const float ai =
                    dec_i[pos + n];

                const float aq =
                    dec_q[pos + n];

                ci +=
                    ai * c +
                    aq * s;

                cq +=
                    aq * c -
                    ai * s;

                const float nc =
                    c * cd -
                    s * sd;

                const float ns =
                    s * cd +
                    c * sd;

                c = nc;
                s = ns;
            }

            coarse_score[d] +=
                ci * ci +
                cq * cq;

            ++coarse_used[d];
        }
    }

    // --------------------------------------------------------
    // Find coarse maximum
    // --------------------------------------------------------

    int best_coarse =
        coarse_first;

    float best_coarse_score =
        -1.0f;

    for (int d = 0;
         d < ncoarse;
         ++d)
    {
        if (coarse_used[d] == 0)
            continue;

        const float score =
            coarse_score[d] /
            (float)coarse_used[d];

        if (score > best_coarse_score)
        {
            best_coarse_score =
                score;

            best_coarse =
                coarse_first + d;
        }
    }

    // --------------------------------------------------------
    // Convert coarse result back to 12-kHz samples.
    //
    // Coarse resolution = 4 samples.
    // --------------------------------------------------------

    const int coarse_delay_rf =
        best_coarse * DECIM;

    // --------------------------------------------------------
    // FINAL REFINEMENT
    //
    // Only +/- 3 samples around coarse result.
    //
    // This is deliberately the ORIGINAL correlation
    // algorithm, so that we retain 1-sample resolution and
    // can directly compare the result with the known-good
    // reference.
    //
    // Maximum candidates:
    //
    //     7 delays × 79 symbols
    //
    // instead of:
    //
    //     1201 delays × 79 symbols
    // --------------------------------------------------------

    constexpr int REFINE_RADIUS = 3;

    const int refine_first =
        std::max(
            0,
            coarse_delay_rf -
            REFINE_RADIUS);

    const int refine_last =
        std::min(
            num_samples - 1,
            coarse_delay_rf +
            REFINE_RADIUS);

    const int nrefine =
        refine_last -
        refine_first +
        1;

    std::vector<float>
        refine_score(
            nrefine,
            0.0f);

    std::vector<int>
        refine_used(
            nrefine,
            0);

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        const float tone_freq =
            6.25f *
            (float)tones[k];

        const float w =
            TWO_PI *
            tone_freq /
            (float)FS;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);

        for (int d = 0;
             d < nrefine;
             ++d)
        {
            const int delay =
                refine_first + d;

            const int start =
                delay +
                k * SPS;

            if (start < 0 ||
                start >= num_samples)
                continue;

            const int count =
                std::min(
                    SPS,
                    num_samples - start);

            if (count <= 0)
                continue;

            float ci = 0.0f;
            float cq = 0.0f;

            float c = 1.0f;
            float s = 0.0f;

            for (int n = 0;
                 n < count;
                 ++n)
            {
                const float x =
                    samples[start + n];

                ci +=
                    x * c;

                cq -=
                    x * s;

                const float nc =
                    c * cd -
                    s * sd;

                const float ns =
                    s * cd +
                    c * sd;

                c = nc;
                s = ns;
            }

            refine_score[d] +=
                ci * ci +
                cq * cq;

            ++refine_used[d];
        }
    }

    // --------------------------------------------------------
    // Final maximum
    // --------------------------------------------------------

    int best_delay =
        refine_first;

    float best_score =
        -1.0f;

    for (int d = 0;
         d < nrefine;
         ++d)
    {
        if (refine_used[d] == 0)
            continue;

        const float score =
            refine_score[d] /
            (float)refine_used[d];

        if (score > best_score)
        {
            best_score =
                score;

            best_delay =
                refine_first + d;
        }
    }

    return
        best_delay /
        (float)FS;
}

// ============================================================
// FT8 DELAY FUNCTION
//
// INCOLLARE QUI LA refine_ft8_delay() ORIGINALE
// ESATTAMENTE COME FORNITA
// ============================================================

// ============================================================
// refine_ft8_delay()
// ============================================================
//
// Pipeline:
//
//     RF
//      |
//      v
//   Downconvert @ freq
//      |
//      v
//   Complex FIR LPF 200 Hz
//      |
//      v
//   FT8 tone correlation
//      |
//      v
//   Delay estimate, 1 sample resolution
//
// FIR:
//   129 taps
//   cutoff = 200 Hz
//   group delay = 64 samples
//
// ============================================================

// ============================================================
// refine_ft8_delay_fir()
//
// Original delay estimator +:
//
//   1. RF -> complex baseband downconversion at 'freq'
//   2. Complex FIR LPF, 200 Hz
//   3. 129 taps -> group delay = 64 samples
//   4. Same +/- 50 ms delay search
//   5. Same 1-sample resolution
//
// FIR:
//   Fs       = 12000 Hz
//   cutoff   = 200 Hz
//   taps     = 129
//   delay    = 64 samples
//
// The FIR group delay is compensated explicitly.
// The final FT8 symbol may be partial, as in the original.
//
// ============================================================

float refine_ft8_delay_fir(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int sample_rate = 12000;
    constexpr int ntones      = 79;
    constexpr int block       = 1920;       // 160 ms @ 12 kHz

    // --------------------------------------------------------
    // FIR parameters
    // --------------------------------------------------------

    constexpr int FIR_TAPS  = 129;
    constexpr int FIR_DELAY = FIR_TAPS / 2; // 64 samples
    constexpr float FIR_CUTOFF = 200.0f;

    // --------------------------------------------------------
    // Delay search
    // --------------------------------------------------------

    const int first_delay =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.05f) *
                sample_rate));

    const int last_delay =
        std::min(
            num_samples - 1,
            (int)std::lround(
                (delay0 + 0.05f) *
                sample_rate));

    const int ndelays =
        last_delay - first_delay + 1;

    // --------------------------------------------------------
    // FIR coefficients
    //
    // Windowed-sinc LPF.
    //
    // 129 taps
    // fc = 200 Hz
    // Fs = 12000 Hz
    //
    // Hann window.
    // --------------------------------------------------------

    float fir[FIR_TAPS];

    constexpr float MY_TWO_PI =
        6.2831853071795864769f;

    const float fc =
        FIR_CUTOFF /
        (float)sample_rate;

    float fir_sum = 0.0f;

    for (int n = 0; n < FIR_TAPS; ++n)
    {
        const float m =
            (float)n -
            (float)FIR_DELAY;

        float h;

        if (fabsf(m) < 1.0e-12f)
        {
            h = 2.0f * fc;
        }
        else
        {
            h =
                sinf(
                    MY_TWO_PI *
                    fc *
                    m) /
                (float)(M_PI * m);
        }

        // Hann window
        const float w =
            0.5f -
            0.5f *
            cosf(
                MY_TWO_PI *
                (float)n /
                (float)(FIR_TAPS - 1));

        fir[n] = h * w;
        fir_sum += fir[n];
    }

    // Normalize DC gain to 1.
    if (fabsf(fir_sum) > 1.0e-20f)
    {
        const float inv_sum =
            1.0f / fir_sum;

        for (int n = 0; n < FIR_TAPS; ++n)
            fir[n] *= inv_sum;
    }

    // --------------------------------------------------------
    // Downconversion oscillator
    //
    // RF signal:
    //
    //   cos((freq + tone) * t)
    //
    // becomes approximately:
    //
    //   exp(+j tone t)
    //
    // after multiplication by exp(-j freq t).
    // --------------------------------------------------------

    const float lo_w =
        MY_TWO_PI *
        freq /
        (float)sample_rate;

    const float lo_cd =
        cosf(lo_w);

    const float lo_sd =
        sinf(lo_w);

    // --------------------------------------------------------
    // We need samples around the complete region searched.
    //
    // Start FIR input 64 samples before first candidate so
    // that the FIR output corresponding to first_delay is
    // available at index 0.
    //
    // We do NOT require the complete 79th symbol to be present.
    // --------------------------------------------------------

    const int input_start =
        std::max(
            0,
            first_delay - FIR_DELAY);

    const int input_end =
        num_samples - 1;

    const int input_count =
        input_end -
        input_start + 1;

    if (input_count <= 0)
        return delay0;

    // --------------------------------------------------------
    // Baseband buffers.
    //
    // bb_i/bb_q:
    //   RF downconverted complex signal
    //
    // fir_i/fir_q:
    //   FIR filtered complex signal
    // --------------------------------------------------------

    std::vector<float> bb_i(input_count);
    std::vector<float> bb_q(input_count);

    std::vector<float> fir_i(input_count);
    std::vector<float> fir_q(input_count);

    // --------------------------------------------------------
    // Downconvert RF -> complex baseband
    // --------------------------------------------------------

    float lo_c = 1.0f;
    float lo_s = 0.0f;

    // Advance oscillator to input_start.
    //
    // This avoids having to calculate sin/cos for every sample
    // during downconversion.
    //
    // First advance to the required absolute sample position.

    for (int n = 0; n < input_start; ++n)
    {
        const float nc =
            lo_c * lo_cd -
            lo_s * lo_sd;

        const float ns =
            lo_s * lo_cd +
            lo_c * lo_sd;

        lo_c = nc;
        lo_s = ns;
    }

    for (int n = 0; n < input_count; ++n)
    {
        const float x =
            samples[input_start + n];

        // Multiply by exp(-j*w*n)
        bb_i[n] =
            x * lo_c;

        bb_q[n] =
            -x * lo_s;

        const float nc =
            lo_c * lo_cd -
            lo_s * lo_sd;

        const float ns =
            lo_s * lo_cd +
            lo_c * lo_sd;

        lo_c = nc;
        lo_s = ns;
    }

    // --------------------------------------------------------
    // Complex FIR
    //
    // y[n] = sum h[k] x[n-k]
    //
    // For the first samples, unavailable samples are simply
    // treated as zero.
    // --------------------------------------------------------

    for (int n = 0; n < input_count; ++n)
    {
        float yi = 0.0f;
        float yq = 0.0f;

        const int kmax =
            std::min(
                FIR_TAPS - 1,
                n);

        for (int k = 0; k <= kmax; ++k)
        {
            const float h =
                fir[k];

            yi +=
                bb_i[n - k] * h;

            yq +=
                bb_q[n - k] * h;
        }

        fir_i[n] = yi;
        fir_q[n] = yq;
    }

    // --------------------------------------------------------
    // Delay score arrays
    // --------------------------------------------------------

    std::vector<float>
        score_total(
            ndelays,
            0.0f);

    std::vector<int>
        used_count(
            ndelays,
            0);

    // --------------------------------------------------------
    // Process all 79 FT8 symbols.
    // --------------------------------------------------------

    for (int k = 0; k < ntones; ++k)
    {
        // Residual baseband tone.
        //
        // The RF carrier 'freq' has already been removed.
        //

        const float tone_freq =
            6.25f *
            (float)tones[k];

        const float w =
            MY_TWO_PI *
            tone_freq /
            (float)sample_rate;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);

        // ----------------------------------------------------
        // exp(-j*w*(block-1))
        //
        // Same quantity used by the original sliding
        // implementation.
        // ----------------------------------------------------

        float end_c = 1.0f;
        float end_s = 0.0f;

        for (int n = 0;
             n < block - 1;
             ++n)
        {
            const float nc =
                end_c * cd +
                end_s * sd;

            const float ns =
                end_s * cd -
                end_c * sd;

            end_c = nc;
            end_s = ns;
        }

        // ----------------------------------------------------
        // First delay
        // ----------------------------------------------------

        float ci = 0.0f;
        float cq = 0.0f;

        const int symbol_rf_start =
            first_delay +
            k * block;

        // Convert RF sample index to FIR-output index.
        //
        // FIR output n corresponds to RF sample:
        //
        //     input_start + n - FIR_DELAY
        //
        // Therefore:
        //
        //     output_index =
        //       symbol_rf_start
        //       - input_start
        //       + FIR_DELAY
        //
        // Wait: because the causal FIR output y[n] contains
        // input samples ending at n, the filtered signal
        // represents input time n-FIR_DELAY.
        //
        // To represent RF sample symbol_rf_start:
        //
        //     n = symbol_rf_start-input_start+FIR_DELAY
        //
        // This explicitly compensates the FIR group delay.
        // ----------------------------------------------------

        const int symbol_start =
            symbol_rf_start -
            input_start +
            FIR_DELAY;

        if (symbol_start >= 0 &&
            symbol_start < input_count)
        {
            const int count =
                std::min(
                    block,
                    input_count -
                    symbol_start);

            float c = 1.0f;
            float s = 0.0f;

            for (int n = 0;
                 n < count;
                 ++n)
            {
                const float ai =
                    fir_i[
                        symbol_start + n];

                const float aq =
                    fir_q[
                        symbol_start + n];

                // Multiply complex baseband by
                // exp(-j*w*n)
                //
                // (ai + j aq)
                // *
                // (c - j s)
                //
                ci +=
                    ai * c +
                    aq * s;

                cq +=
                    aq * c -
                    ai * s;

                const float nc =
                    c * cd -
                    s * sd;

                const float ns =
                    s * cd +
                    c * sd;

                c = nc;
                s = ns;
            }

            score_total[0] +=
                ci * ci +
                cq * cq;

            ++used_count[0];
        }

        // ----------------------------------------------------
        // Sliding delay search
        // ----------------------------------------------------

        for (int d = 1;
             d < ndelays;
             ++d)
        {
            const int rf_start =
                first_delay +
                d +
                k * block;

            // Corresponding FIR output.
            const int start =
                rf_start -
                input_start +
                FIR_DELAY;

            if (start >= input_count)
                break;

            // Determine how many samples are actually
            // available. This preserves the original
            // handling of the final partial symbol.
            const int available =
                input_count -
                start;

            const int count =
                std::min(
                    block,
                    available);

            if (count <= 0)
                break;

            // For the first delay the correlation state
            // contains 'block' samples. For a partial symbol
            // there is no useful sliding state beyond the
            // available data, so recompute it.
            //
            // Normally all symbols except the final one use
            // the fast sliding path.

            if (count < block)
            {
                ci = 0.0f;
                cq = 0.0f;

                float c = 1.0f;
                float s = 0.0f;

                for (int n = 0;
                     n < count;
                     ++n)
                {
                    const float ai =
                        fir_i[start + n];

                    const float aq =
                        fir_q[start + n];

                    ci +=
                        ai * c +
                        aq * s;

                    cq +=
                        aq * c -
                        ai * s;

                    const float nc =
                        c * cd -
                        s * sd;

                    const float ns =
                        s * cd +
                        c * sd;

                    c = nc;
                    s = ns;
                }

                score_total[d] +=
                    ci * ci +
                    cq * cq;

                ++used_count[d];

                continue;
            }

            // ------------------------------------------------
            // Normal full-block sliding update.
            // ------------------------------------------------

            const float ai_old =
                fir_i[start - 1];

            const float aq_old =
                fir_q[start - 1];

            const float ai_new =
                fir_i[start + block - 1];

            const float aq_new =
                fir_q[start + block - 1];

            // Remove oldest sample.
            //
            // The old sample had phase exp(-j*w*0)
            // relative to the old window.
            //
            // First rotate the old correlation by exp(+j*w)
            // so that the new window has the correct phase
            // reference.
            //
            // Then add the new sample at phase
            // exp(-j*w*(block-1)).
            //

            const float rotated_re =
                (ci - ai_old) * cd -
                cq * sd;

            const float rotated_im =
                (ci - ai_old) * sd +
                cq * cd;

            ci =
                rotated_re +
                ai_new * end_c +
                aq_new * (-end_s);

            cq =
                rotated_im +
                aq_new * end_c -
                ai_new * end_s;

            score_total[d] +=
                ci * ci +
                cq * cq;

            ++used_count[d];
        }
    }

    // --------------------------------------------------------
    // Find best delay
    // --------------------------------------------------------

    int best_delay =
        first_delay;

    float best_score =
        -1.0f;

    for (int d = 0;
         d < ndelays;
         ++d)
    {
        if (used_count[d] == 0)
            continue;

        const float score =
            score_total[d] /
            (float)used_count[d];

        if (score > best_score)
        {
            best_score =
                score;

            best_delay =
                first_delay + d;
        }
    }

    return best_delay /
           (float)sample_rate;
}


float refine_ft8_delay_original(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int sample_rate = 12000;
    constexpr int ntones = 79;

    const int block = blocksize(sample_rate);

    const int first_delay =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.05f) * sample_rate));

    const int last_delay =
        std::min(
            num_samples - 1,
            (int)std::lround(
                (delay0 + 0.05f) * sample_rate));

    const int ndelays =
        last_delay - first_delay + 1;

    // constexpr float TWO_PI =
    //     6.2831853071795864769f;

    std::vector<float> score_total(ndelays, 0.0f);
    std::vector<int> used_count(ndelays, 0);

    for (int k = 0; k < ntones; ++k)
    {
        const float freq_hz =
            freq + 6.25f * (float)tones[k];

        const float w =
            TWO_PI * freq_hz /
            (float)sample_rate;

        const float cd = cosf(w);
        const float sd = sinf(w);

        float end_c = 1.0f;
        float end_s = 0.0f;

        for (int n = 0; n < block - 1; ++n)
        {
            const float nc =
                end_c * cd + end_s * sd;

            const float ns =
                end_s * cd - end_c * sd;

            end_c = nc;
            end_s = ns;
        }

        float ci = 0.0f;
        float cq = 0.0f;

        const int start0 =
            first_delay + k * block;

        if (start0 >= 0 &&
            start0 < num_samples)
        {
            const int count =
                std::min(
                    block,
                    num_samples - start0);

            float c = 1.0f;
            float s = 0.0f;

            for (int n = 0; n < count; ++n)
            {
                const float x =
                    samples[start0 + n];

                ci += x * c;
                cq -= x * s;

                const float nc =
                    c * cd -
                    s * sd;

                const float ns =
                    s * cd +
                    c * sd;

                c = nc;
                s = ns;
            }

            score_total[0] +=
                ci * ci + cq * cq;

            ++used_count[0];
        }

        for (int d = 1; d < ndelays; ++d)
        {
            const int start =
                first_delay + d + k * block;

            if (start >= num_samples)
                break;

            const int new_index =
                start + block - 1;

            if (new_index >= num_samples)
                break;

            const float x_old =
                samples[start - 1];

            const float x_new =
                samples[new_index];

            const float r =
                ci - x_old;

            const float im =
                cq;

            const float rotated_re =
                r * cd - im * sd;

            const float rotated_im =
                r * sd + im * cd;

            ci =
                rotated_re + x_new * end_c;

            cq =
                rotated_im + x_new * end_s;

            score_total[d] +=
                ci * ci + cq * cq;

            ++used_count[d];
        }
    }

    int best_delay = first_delay;
    float best_score = -1.0f;

    for (int d = 0; d < ndelays; ++d)
    {
        if (used_count[d] == 0)
            continue;

        const float score =
            score_total[d] /
            (float)used_count[d];

        if (score > best_score)
        {
            best_score = score;
            best_delay =
                first_delay + d;
        }
    }

    return best_delay /
           (float)sample_rate;
}

// ============================================================
// WAV LOADER
//
// Same interface on both platforms:
//
//   ESP32-S3 Arduino : LittleFS
//   Linux / Ubuntu   : std::ifstream
//
// Supported format:
//   RIFF/WAVE
//   PCM
//   16-bit
//   mono or multi-channel (first channel is used)
// ============================================================

#if defined(ARDUINO)

static float* load_wav(
    const char* path,
    int* out_num_samples,
    int* out_num_channels,
    int* out_sample_rate)
{
    File f = LittleFS.open(path, "r");

    if (!f)
    {
        printf("ERROR: cannot open %s\n", path);
        return nullptr;
    }

    uint8_t riff[12];

    if (f.read(riff, 12) != 12)
    {
        printf("ERROR: WAV header too short\n");
        f.close();
        return nullptr;
    }

    if (memcmp(riff, "RIFF", 4) != 0 ||
        memcmp(riff + 8, "WAVE", 4) != 0)
    {
        printf("ERROR: not a RIFF/WAVE file\n");
        f.close();
        return nullptr;
    }

    uint16_t audio_format = 0;
    uint16_t num_channels = 0;
    uint32_t sample_rate = 0;
    uint16_t bits_per_sample = 0;

    uint32_t data_pos = 0;
    uint32_t data_size = 0;

    bool have_fmt = false;
    bool have_data = false;

    while (f.available())
    {
        uint8_t chunk_header[8];

        if (f.read(chunk_header, 8) != 8)
            break;

        const uint32_t chunk_size =
            chunk_header[4] |
            ((uint32_t)chunk_header[5] << 8) |
            ((uint32_t)chunk_header[6] << 16) |
            ((uint32_t)chunk_header[7] << 24);

        const uint32_t chunk_start = f.position();

        if (memcmp(chunk_header, "fmt ", 4) == 0)
        {
            if (chunk_size < 16)
            {
                printf("ERROR: invalid fmt chunk\n");
                f.close();
                return nullptr;
            }

            uint8_t fmt[16];

            if (f.read(fmt, 16) != 16)
            {
                printf("ERROR: cannot read fmt\n");
                f.close();
                return nullptr;
            }

            audio_format =
                fmt[0] |
                ((uint16_t)fmt[1] << 8);

            num_channels =
                fmt[2] |
                ((uint16_t)fmt[3] << 8);

            sample_rate =
                fmt[4] |
                ((uint32_t)fmt[5] << 8) |
                ((uint32_t)fmt[6] << 16) |
                ((uint32_t)fmt[7] << 24);

            bits_per_sample =
                fmt[14] |
                ((uint16_t)fmt[15] << 8);

            have_fmt = true;
        }
        else if (memcmp(chunk_header, "data", 4) == 0)
        {
            data_pos = f.position();
            data_size = chunk_size;
            have_data = true;
        }

        f.seek(chunk_start + chunk_size);

        if (chunk_size & 1)
            f.seek(f.position() + 1);

        if (have_fmt && have_data)
            break;
    }

    if (!have_fmt || !have_data)
    {
        printf("ERROR: missing fmt/data chunk\n");
        f.close();
        return nullptr;
    }

    if (audio_format != 1)
    {
        printf("ERROR: unsupported WAV format %u\n",
               audio_format);
        f.close();
        return nullptr;
    }

    if (bits_per_sample != 16)
    {
        printf("ERROR: only 16-bit PCM supported, got %u\n",
               bits_per_sample);
        f.close();
        return nullptr;
    }

    if (num_channels < 1)
    {
        printf("ERROR: invalid channel count\n");
        f.close();
        return nullptr;
    }

    const int bytes_per_sample =
        bits_per_sample / 8;

    const int frame_bytes =
        bytes_per_sample * num_channels;

    const int num_samples =
        data_size / frame_bytes;

    float* samples =
        (float*)malloc(
            (size_t)num_samples * sizeof(float));

    if (!samples)
    {
        printf("ERROR: malloc failed for %d samples\n",
               num_samples);
        f.close();
        return nullptr;
    }

    f.seek(data_pos);

    for (int n = 0; n < num_samples; ++n)
    {
        uint8_t b[2];

        if (f.read(b, 2) != 2)
        {
            printf("ERROR: unexpected EOF\n");
            free(samples);
            f.close();
            return nullptr;
        }

        const int16_t v =
            (int16_t)(
                b[0] |
                ((uint16_t)b[1] << 8));

        samples[n] =
            (float)v / 32768.0f;

        // Skip remaining channels.
        for (int ch = 1;
             ch < num_channels;
             ++ch)
        {
            f.seek(
                f.position() +
                bytes_per_sample);
        }
    }

    f.close();

    *out_num_samples  = num_samples;
    *out_num_channels = num_channels;
    *out_sample_rate  = (int)sample_rate;

    return samples;
}

#else

static uint32_t read_u32_le(
    const uint8_t* p)
{
    return
        (uint32_t)p[0] |
        ((uint32_t)p[1] << 8) |
        ((uint32_t)p[2] << 16) |
        ((uint32_t)p[3] << 24);
}

static uint16_t read_u16_le(
    const uint8_t* p)
{
    return
        (uint16_t)p[0] |
        ((uint16_t)p[1] << 8);
}

static float* load_wav(
    const char* path,
    int* out_num_samples,
    int* out_num_channels,
    int* out_sample_rate)
{
    // On Linux the Arduino paths "/file.wav" are interpreted
    // as files in the current directory. This lets the exact
    // same test names be used on both platforms.
    const char* open_path = path;

    std::ifstream f(
        open_path,
        std::ios::binary);

    if (!f && path[0] == '/')
    {
        open_path = path + 1;
        f.clear();
        f.open(
            open_path,
            std::ios::binary);
    }

    if (!f)
    {
        printf("ERROR: cannot open %s\n", path);
        return nullptr;
    }

    uint8_t riff[12];

    f.read(
        (char*)riff,
        sizeof(riff));

    if (f.gcount() != (std::streamsize)sizeof(riff))
    {
        printf("ERROR: WAV header too short\n");
        return nullptr;
    }

    if (memcmp(riff, "RIFF", 4) != 0 ||
        memcmp(riff + 8, "WAVE", 4) != 0)
    {
        printf("ERROR: not a RIFF/WAVE file\n");
        return nullptr;
    }

    uint16_t audio_format = 0;
    uint16_t num_channels = 0;
    uint32_t sample_rate = 0;
    uint16_t bits_per_sample = 0;

    std::streamoff data_pos = 0;
    uint32_t data_size = 0;

    bool have_fmt = false;
    bool have_data = false;

    while (f)
    {
        uint8_t chunk_header[8];

        f.read(
            (char*)chunk_header,
            sizeof(chunk_header));

        if (f.gcount() !=
            (std::streamsize)sizeof(chunk_header))
            break;

        const uint32_t chunk_size =
            read_u32_le(chunk_header + 4);

        const std::streamoff chunk_start =
            f.tellg();

        if (memcmp(chunk_header, "fmt ", 4) == 0)
        {
            if (chunk_size < 16)
            {
                printf("ERROR: invalid fmt chunk\n");
                return nullptr;
            }

            uint8_t fmt[16];

            f.read(
                (char*)fmt,
                sizeof(fmt));

            if (f.gcount() !=
                (std::streamsize)sizeof(fmt))
            {
                printf("ERROR: cannot read fmt\n");
                return nullptr;
            }

            audio_format =
                read_u16_le(fmt + 0);

            num_channels =
                read_u16_le(fmt + 2);

            sample_rate =
                read_u32_le(fmt + 4);

            bits_per_sample =
                read_u16_le(fmt + 14);

            have_fmt = true;
        }
        else if (memcmp(chunk_header, "data", 4) == 0)
        {
            data_pos = f.tellg();
            data_size = chunk_size;
            have_data = true;
        }

        f.seekg(
            chunk_start +
            (std::streamoff)chunk_size);

        if (chunk_size & 1)
            f.seekg(
                f.tellg() + std::streamoff(1));

        if (have_fmt && have_data)
            break;
    }

    if (!have_fmt || !have_data)
    {
        printf("ERROR: missing fmt/data chunk\n");
        return nullptr;
    }

    if (audio_format != 1)
    {
        printf("ERROR: unsupported WAV format %u\n",
               audio_format);
        return nullptr;
    }

    if (bits_per_sample != 16)
    {
        printf("ERROR: only 16-bit PCM supported, got %u\n",
               bits_per_sample);
        return nullptr;
    }

    if (num_channels < 1)
    {
        printf("ERROR: invalid channel count\n");
        return nullptr;
    }

    const int bytes_per_sample =
        bits_per_sample / 8;

    const int frame_bytes =
        bytes_per_sample * num_channels;

    const int num_samples =
        (int)(data_size / frame_bytes);

    float* samples =
        (float*)malloc(
            (size_t)num_samples * sizeof(float));

    if (!samples)
    {
        printf("ERROR: malloc failed for %d samples\n",
               num_samples);
        return nullptr;
    }

    f.clear();
    f.seekg(data_pos);

    for (int n = 0; n < num_samples; ++n)
    {
        uint8_t b[2];

        f.read(
            (char*)b,
            2);

        if (f.gcount() != 2)
        {
            printf("ERROR: unexpected EOF\n");
            free(samples);
            return nullptr;
        }

        const int16_t v =
            (int16_t)(
                b[0] |
                ((uint16_t)b[1] << 8));

        samples[n] =
            (float)v / 32768.0f;

        // Skip remaining channels.
        if (num_channels > 1)
        {
            f.seekg(
                f.tellg() +
                std::streamoff(
                    (num_channels - 1) *
                    bytes_per_sample));
        }
    }

    *out_num_samples  = num_samples;
    *out_num_channels = num_channels;
    *out_sample_rate  = (int)sample_rate;

    return samples;
}

#endif


// ============================================================
// Tone estimator
//
// Ricava il tone FT8 dominante per ogni simbolo.
// Non viene usato per misurare il delay: serve solamente
// per fornire alla refine_ft8_delay() i 79 tones del WAV.
//
// Frequenze:
//     freq + tone * 6.25 Hz
//
// Viene utilizzata una correlazione complessa di 160 ms.
// ============================================================

static uint8_t estimate_tone(
    const float* samples,
    int num_samples,
    int start,
    float freq)
{
    float best_power = -1.0f;
    int best_tone = 0;

    for (int tone = 0; tone < 8; ++tone)
    {
        const float f =
            freq +
            tone * TONE_SPACING;

        const float w =
            2.0f * PI * f /
            SAMPLE_RATE;

        const float c =
            cosf(w);

        const float s =
            sinf(w);

        float ci = 0.0f;
        float cq = 0.0f;

        float phase_c = 1.0f;
        float phase_s = 0.0f;

        const int count =
            std::min(
                SPS,
                num_samples - start);

        if (count <= 0)
            continue;

        for (int n = 0; n < count; ++n)
        {
            const float x =
                samples[start + n];

            ci += x * phase_c;
            cq -= x * phase_s;

            const float nc =
                phase_c * c -
                phase_s * s;

            const float ns =
                phase_s * c +
                phase_c * s;

            phase_c = nc;
            phase_s = ns;
        }

        const float power =
            ci * ci + cq * cq;

        if (power > best_power)
        {
            best_power = power;
            best_tone = tone;
        }
    }

    return (uint8_t)best_tone;
}


// ============================================================
// Extract all 79 tones
//
// ATTENZIONE:
//
// delay0 è la posizione nota del messaggio.
// I tone vengono cercati partendo da delay0.
//
// Questo è deliberatamente separato dal test del delay.
// ============================================================

static bool extract_ft8_tones(
    const float* samples,
    int num_samples,
    float delay0,
    float freq,
    uint8_t tones[NTONES])
{
    const int start =
        (int)lroundf(
            delay0 * SAMPLE_RATE);

        
    if (start < 0 || start >= num_samples)
    {
        printf(
            "ERROR: invalid FT8 start position");

        return false;
    }

    for (int k = 0; k < NTONES; ++k)
    {
        const int symbol_start =
            start + k * SPS;

        tones[k] =
            estimate_tone(
                samples,
                num_samples,
                symbol_start,
                freq);
    }

    return true;
}


// ============================================================
// Print tones
// ============================================================

static void print_tones(
    const uint8_t tones[NTONES])
{
    printf(
        "Detected FT8 tones:");

    for (int k = 0; k < NTONES; ++k)
    {
        printf(
            "%2d:%d ",
            k,
            tones[k]);

        if ((k % 13) == 12)
            printf("\n");
    }

    printf("\n");
}


// ============================================================
// TEST ONE WAV
// ============================================================

static void run_one_test(
    const char* filename,
    float true_delay)
{
    printf("\n");
    printf("============================================================\n");
    printf("FILE: %s\n", filename);
    printf("TRUE DELAY: %.6f s\n", true_delay);
    printf("============================================================\n");

    int num_samples = 0;
    int num_channels = 0;
    int sample_rate = 0;

#if defined(ARDUINO)
    const uint32_t heap_before =
        ESP.getFreeHeap();
#else
    const uint64_t heap_before = 0;
#endif

    float* samples =
        load_wav(
            filename,
            &num_samples,
            &num_channels,
            &sample_rate);

    if (!samples)
    {
        printf("LOAD FAILED\n");
        return;
    }

#if defined(ARDUINO)
    const uint32_t heap_after_load =
        ESP.getFreeHeap();
#else
    const uint64_t heap_after_load = 0;
#endif

    printf("Samples          : %d\n", num_samples);
    printf("Duration         : %.3f s\n",
           (float)num_samples / sample_rate);
    printf("Channels         : %d\n", num_channels);
    printf("Sample rate      : %d Hz\n", sample_rate);

#if defined(ARDUINO)
    printf("Heap before load : %u\n",
           (unsigned)heap_before);
    printf("Heap after load  : %u\n",
           (unsigned)heap_after_load);
#else
    printf("Heap before load : N/A (Linux)\n");
    printf("Heap after load  : N/A (Linux)\n");
#endif

    if (sample_rate != SAMPLE_RATE ||
        num_channels != 1)
    {
        printf("ERROR: WAV format mismatch\n");
        free(samples);
        return;
    }

    // --------------------------------------------------------
    // Ricava i 79 tones dal WAV
    // --------------------------------------------------------

    uint8_t tones[NTONES];

    if (!extract_ft8_tones(
            samples,
            num_samples,
            true_delay,
            TEST_FREQ,
            tones))
    {
        free(samples);
        return;
    }

    print_tones(tones);

    // --------------------------------------------------------
    // TEST DELLA refine_ft8_delay_fir()
    // --------------------------------------------------------

    printf("\n");
    printf("Running ORIGINAL refine_ft8_delay_fir()...\n");

    const uint64_t t0 =
        portable_micros();

    const float estimated_delay =
        refine_ft8_delay_fir(
            samples,
            num_samples,
            tones,
            true_delay,
            TEST_FREQ,
            0);

    const uint64_t elapsed_us =
        portable_micros() - t0;

    const float error_s =
        estimated_delay - true_delay;

    const float error_ms =
        error_s * 1000.0f;

    const float error_samples =
        error_s * SAMPLE_RATE;

    printf("\n");
    printf("------------------------------------------------------------\n");

    printf("True delay       : %10.6f s\n",
           true_delay);

    printf("Estimated delay  : %10.6f s\n",
           estimated_delay);

    printf("Error            : %+10.3f ms\n",
           error_ms);

    printf("Error            : %+10.3f samples\n",
           error_samples);

    printf("Execution time   : %10.3f ms\n",
           elapsed_us / 1000.0);

#if defined(ARDUINO)
    printf("Free heap        : %u bytes\n",
           (unsigned)ESP.getFreeHeap());
#endif

    printf("------------------------------------------------------------\n");

    // --------------------------------------------------------
    // Liberiamo il WAV PRIMA di caricare il successivo.
    // --------------------------------------------------------

    free(samples);

#if defined(ARDUINO)
    printf("After free() heap: %u bytes\n",
           (unsigned)ESP.getFreeHeap());
    printf("WAV released.\n");
#else
    printf("WAV released.\n");
#endif
}


// ============================================================
// TEST SU ENTRAMBE LE PIATTAFORME
// ============================================================

static void run_all_tests()
{
    printf("\n");
    printf("############################################################\n");
    printf(" ORIGINAL refine_ft8_delay() - REAL WAV BENCHMARK\n");
    printf("############################################################\n");

#if defined(ARDUINO)
    printf("Platform         : ESP32-S3 / Arduino\n");
    printf("CPU frequency    : %u MHz\n",
           (unsigned)ESP.getCpuFreqMHz());
#else
    printf("Platform         : Linux / Ubuntu\n");
#endif

    printf("Sample rate      : %d Hz\n",
           SAMPLE_RATE);

    printf("Symbol duration  : %d samples = %.3f ms\n",
           SPS,
           1000.0f * SPS / SAMPLE_RATE);

    printf("Test frequency   : %.3f Hz\n",
           TEST_FREQ);

#if defined(ARDUINO)
    printf("Free heap        : %u bytes\n",
           (unsigned)ESP.getFreeHeap());

    printf("Free PSRAM       : %u bytes\n",
           (unsigned)ESP.getFreePsram());

    printf("\n");
    printf("Mounting LittleFS...\n");

    if (!LittleFS.begin(true))
    {
        printf("ERROR: LittleFS.begin() failed\n");
        return;
    }

    printf("LittleFS mounted.\n");
#endif

    // --------------------------------------------------------
    // TEST 1
    // --------------------------------------------------------

    run_one_test(
        "/test_real_ft8_D1500.wav",
        1.500f);

    // --------------------------------------------------------
    // TEST 2
    // --------------------------------------------------------

    run_one_test(
        "/test_real_ft8_D2000.wav",
        2.000f);

    // --------------------------------------------------------
    // TEST 3
    // --------------------------------------------------------

    run_one_test(
        "/test_real_ft8_D2500.wav",
        2.500f);

    printf("\n");
    printf("############################################################\n");
    printf("BENCHMARK COMPLETE\n");

#if defined(ARDUINO)
    printf("Final free heap : %u bytes\n",
           (unsigned)ESP.getFreeHeap());

    printf("Final free PSRAM: %u bytes\n",
           (unsigned)ESP.getFreePsram());
#endif

    printf("############################################################\n");
}


#if defined(ARDUINO)

// ============================================================
// ARDUINO ENTRY POINT
// ============================================================

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

// ============================================================
// LINUX ENTRY POINT
// ============================================================

int main()
{
    run_all_tests();
    return 0;
}

#endif

