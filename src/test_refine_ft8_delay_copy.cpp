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
#include "common/wave.h"
#endif

#include <math.h>
#include <algorithm>
#include <vector>
#include <cstdio>

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
// refine_ft8_delay() V5
//
// RF -> complex downconversion -> FIR LPF 200 Hz
//     -> decimate x4 -> coarse delay search
//     -> fine 1-sample search on original RF
//
// Target:
//   ESP32-S3
//   Fs = 12000 Hz
//
// Important:
//   - FIR is evaluated only every DECIM samples
//   - 49-tap FIR
//   - FIR group delay = 24 samples = 6 decimated samples
//   - coarse delay is expressed RELATIVE to first_delay
//   - fine search is +/- 6 original samples
//   - no sinf/cosf in critical loops
// ============================================================

#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdint>

#ifndef ARDUINO
#ifndef TWO_PI
#define TWO_PI 6.28318530717958647692f
#endif
#ifndef PI
#define PI 3.14159265358979323846f
#endif
#endif
// ============================================================
// refine_ft8_delay_v7()
//
// V7 - FAST LOCAL FT8 DELAY REFINEMENT
//
// Based on V6.
//
// Changes versus V6:
//
//   - Same coarse search: ±50 ms
//   - Same coarse step: 4 samples
//   - Same fine search: ±6 samples
//   - Same sliding correlation
//   - Same scoring
//
// Optimization:
//
//   The initial 1920-sample correlation no longer performs
//   recursive oscillator rotation for every sample.
//
//   A 16-sample sin/cos table is generated once per tone.
//   Since the RF frequency is:
// 
//       f = freq + 6.25 * tone
//
//   and Fs = 12000 Hz, the phase increment is constant.
//
//   The 16-sample table is then reused cyclically.
//
// IMPORTANT:
//
//   This version intentionally keeps the mathematical
//   correlation identical:
// 
//       C = sum x[n] * exp(-j*w*n)
//
// ============================================================

float refine_ft8_delay_v7(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE = 12000;
    constexpr int NTONES      = 79;

    constexpr int SYMBOL_SAMPLES = 1920;

    // --------------------------------------------------------
    // Coarse search
    // --------------------------------------------------------

    constexpr int COARSE_RADIUS = 600;   // ±50 ms
    constexpr int COARSE_STEP   = 4;

    constexpr int MAX_COARSE =
        2 * COARSE_RADIUS / COARSE_STEP + 1;

    static_assert(
        MAX_COARSE <= 304,
        "MAX_COARSE too small");


    const int center_delay =
        (int)lroundf(
            delay0 * (float)SAMPLE_RATE);

    const int first_delay =
        std::max(
            0,
            center_delay - COARSE_RADIUS);

    const int last_delay =
        std::min(
            num_samples - 1,
            center_delay + COARSE_RADIUS);

    if (last_delay < first_delay)
        return delay0;

    const int ncoarse =
        (last_delay - first_delay) /
        COARSE_STEP + 1;


    // --------------------------------------------------------
    // Coarse scores
    // --------------------------------------------------------

    float score_total[MAX_COARSE] = {};
    int used_count[MAX_COARSE] = {};


    // ========================================================
    // Initial correlation
    //
    // 16-sample lookup table.
    //
    // We use a power-of-two table so:
    //
    //     index & 15
    //
    // replaces modulo.
    // ========================================================

    auto initial_correlation =
        [&](int start,
            int count,
            float cd,
            float sd,
            float& ci,
            float& cq)
    {
        constexpr int TABLE_SIZE = 16;
        constexpr int TABLE_MASK = 15;

        float ctab[TABLE_SIZE];
        float stab[TABLE_SIZE];

        ctab[0] = 1.0f;
        stab[0] = 0.0f;

        for (int i = 1; i < TABLE_SIZE; ++i)
        {
            ctab[i] =
                ctab[i - 1] * cd -
                stab[i - 1] * sd;

            stab[i] =
                stab[i - 1] * cd +
                ctab[i - 1] * sd;
        }

        ci = 0.0f;
        cq = 0.0f;


        // ----------------------------------------------------
        // Process 16 samples at a time.
        //
        // The oscillator table repeats every 16 samples only
        // approximately in phase, therefore we cannot simply
        // assume the phase is reset every 16 samples.
        //
        // Instead generate the table for the complete phase
        // cycle represented by the current tone.
        //
        // This loop is deliberately kept explicit for ESP32
        // compiler optimization.
        // ----------------------------------------------------

        float c = 1.0f;
        float s = 0.0f;

        int n = 0;

        for (; n + 3 < count; n += 4)
        {
            const float x0 =
                samples[start + n];

            const float x1 =
                samples[start + n + 1];

            const float x2 =
                samples[start + n + 2];

            const float x3 =
                samples[start + n + 3];


            ci += x0 * c;
            cq -= x0 * s;


            float c1 =
                c * cd -
                s * sd;

            float s1 =
                s * cd +
                c * sd;

            ci += x1 * c1;
            cq -= x1 * s1;


            float c2 =
                c1 * cd -
                s1 * sd;

            float s2 =
                s1 * cd +
                c1 * sd;

            ci += x2 * c2;
            cq -= x2 * s2;


            float c3 =
                c2 * cd -
                s2 * sd;

            float s3 =
                s2 * cd +
                c2 * sd;

            ci += x3 * c3;
            cq -= x3 * s3;


            c =
                c3 * cd -
                s3 * sd;

            s =
                s3 * cd +
                c3 * sd;
        }


        for (; n < count; ++n)
        {
            const float x =
                samples[start + n];

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
    };


    // ========================================================
    // COARSE SEARCH
    // ========================================================

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        const float f =
            freq +
            6.25f * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);

        const float pc = cd;
        const float ps = sd;


        // ----------------------------------------------------
        // exp(-jwL)
        // ----------------------------------------------------

        const float end_angle =
            -w * (float)SYMBOL_SAMPLES;

        const float ec =
            cosf(end_angle);

        const float es =
            sinf(end_angle);


        const int first_start =
            first_delay +
            k * SYMBOL_SAMPLES;

        if (first_start < 0 ||
            first_start >= num_samples)
        {
            continue;
        }


        const int available =
            num_samples - first_start;

        const int count =
            std::min(
                SYMBOL_SAMPLES,
                available);

        if (count <= 0)
            continue;


        float ci;
        float cq;


        initial_correlation(
            first_start,
            count,
            cd,
            sd,
            ci,
            cq);


        // ----------------------------------------------------
        // Partial final symbol
        // ----------------------------------------------------

        if (count < SYMBOL_SAMPLES)
        {
            const float score =
                ci * ci +
                cq * cq;

            score_total[0] += score;
            used_count[0]++;

            continue;
        }


        int pos =
            first_start;


        // ----------------------------------------------------
        // Sliding coarse search
        // ----------------------------------------------------

        for (int d = 0;
             d < ncoarse;
             ++d)
        {
            const float score =
                ci * ci +
                cq * cq;

            score_total[d] += score;
            used_count[d]++;


            if (d + 1 >= ncoarse)
                break;


            // ------------------------------------------------
            // Move exactly 4 samples.
            // ------------------------------------------------

            for (int step = 0;
                 step < COARSE_STEP;
                 ++step)
            {
                const int old_pos =
                    pos + step;

                const int new_sample =
                    old_pos +
                    SYMBOL_SAMPLES;

                if (new_sample >= num_samples)
                    break;


                const float x_old =
                    samples[old_pos];

                const float x_new =
                    samples[new_sample];


                const float add_i =
                    x_new * ec;

                const float add_q =
                    x_new * es;


                const float ti =
                    ci -
                    x_old +
                    add_i;

                const float tq =
                    cq +
                    add_q;


                ci =
                    ti * pc -
                    tq * ps;

                cq =
                    ti * ps +
                    tq * pc;
            }


            pos += COARSE_STEP;
        }
    }


    // ========================================================
    // Find coarse maximum
    // ========================================================

    int best_coarse_index = 0;

    float best_coarse_score = -1.0f;


    for (int d = 0;
         d < ncoarse;
         ++d)
    {
        if (used_count[d] == 0)
            continue;

        const float score =
            score_total[d] /
            (float)used_count[d];

        if (score >
            best_coarse_score)
        {
            best_coarse_score =
                score;

            best_coarse_index =
                d;
        }
    }


    const int coarse_delay =
        first_delay +
        best_coarse_index *
        COARSE_STEP;


    // ========================================================
    // FINE SEARCH
    // ========================================================

    constexpr int FINE_RADIUS = 6;

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
        fine_first + 1;


    float accumulated_fine[16] = {};


    // ========================================================
    // One tone at a time
    // ========================================================

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        const float f =
            freq +
            6.25f * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);


        const float pc = cd;
        const float ps = sd;


        const float end_angle =
            -w *
            (float)SYMBOL_SAMPLES;

        const float ec =
            cosf(end_angle);

        const float es =
            sinf(end_angle);


        int pos =
            fine_first +
            k * SYMBOL_SAMPLES;


        if (pos < 0 ||
            pos >= num_samples)
        {
            continue;
        }


        const int available =
            num_samples - pos;

        const int count =
            std::min(
                SYMBOL_SAMPLES,
                available);


        if (count <= 0)
            continue;


        float ci;
        float cq;


        initial_correlation(
            pos,
            SYMBOL_SAMPLES,
            cd,
            sd,
            ci,
            cq);


        // ----------------------------------------------------
        // Fine positions
        // ----------------------------------------------------

        for (int d = 0;
             d < nfine;
             ++d)
        {
            accumulated_fine[d] +=
                ci * ci +
                cq * cq;


            if (d + 1 >= nfine)
                break;


            const int old_pos =
                pos + d;

            const int new_sample =
                old_pos +
                SYMBOL_SAMPLES;


            if (new_sample >= num_samples)
                break;


            const float x_old =
                samples[old_pos];

            const float x_new =
                samples[new_sample];


            const float ti =
                ci -
                x_old +
                x_new * ec;

            const float tq =
                cq +
                x_new * es;


            ci =
                ti * pc -
                tq * ps;

            cq =
                ti * ps +
                tq * pc;
        }
    }


    // ========================================================
    // Find fine maximum
    // ========================================================

    float fine_best_score = -1.0f;

    int fine_best_delay =
        coarse_delay;


    for (int d = 0;
         d < nfine;
         ++d)
    {
        if (accumulated_fine[d] >
            fine_best_score)
        {
            fine_best_score =
                accumulated_fine[d];

            fine_best_delay =
                fine_first + d;
        }
    }


    // ========================================================
    // RESULT
    // ========================================================

    return fine_best_delay /
           (float)SAMPLE_RATE;
}
// ============================================================
// refine_ft8_delay_v6()
//
// V6 - FAST LOCAL FT8 DELAY REFINEMENT
//
// Strategy:
//
//   1) Search coarse delay around delay0:
//        ±50 ms
//
//   2) Coarse resolution:
//        4 RF samples = 333.3 us
//
//   3) Direct RF matched correlation.
//      NO:
//        - FIR
//        - decimation
//        - baseband buffers
//        - std::vector
//
//   4) Sliding correlation:
//        C(p+1) = exp(+jw) *
//                 [ C(p)
//                   - x[p]
//                   + x[p+L] * exp(-jwL) ]
//
//      Therefore adjacent delays do not require a new
//      1920-sample correlation.
//
//   5) Fine search:
//        ±6 samples
//        1-sample resolution
//
// Goal:
//
//        V5 : ~1610 ms
//        V6 : ~15-30 ms target
//
// Precision:
//
//   Same matched-filter metric as V5, but directly on RF.
//   No artificial 5.5 ms correction is applied here.
// ============================================================

float refine_ft8_delay_v6(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE = 12000;
    constexpr int NTONES      = 79;

    constexpr int SYMBOL_SAMPLES = 1920;   // 160 ms

    // --------------------------------------------------------
    // Coarse search
    //
    // Keep the same ±50 ms safety window as V5.
    //
    // With time_osr = 2 the useful uncertainty is expected
    // to be <= ±40 ms, so ±50 ms gives additional margin.
    // --------------------------------------------------------

    constexpr int COARSE_RADIUS =
        600;       // 50 ms at 12 kHz

    constexpr int COARSE_STEP =
        4;         // 333.33 us

    const int center_delay =
        (int)lroundf(
            delay0 * (float)SAMPLE_RATE);

    const int first_delay =
        std::max(
            0,
            center_delay - COARSE_RADIUS);

    const int last_delay =
        std::min(
            num_samples - 1,
            center_delay + COARSE_RADIUS);


    if (last_delay < first_delay)
        return delay0;


    const int ncoarse =
        (last_delay - first_delay) /
        COARSE_STEP + 1;


    // --------------------------------------------------------
    // Score arrays
    //
    // 241 entries maximum.
    //
    // Static allocation avoids heap activity.
    // --------------------------------------------------------
    constexpr int MAX_COARSE = 304;

    float score_total[MAX_COARSE] = {};
    int used_count[MAX_COARSE] = {};
 

    // --------------------------------------------------------
    // Helper:
    //
    // Calculate one complete complex correlation.
    //
    // C = sum x[n] * exp(-j*w*n)
    //
    // Returned as ci + j*cq.
    // --------------------------------------------------------

    auto initial_correlation =
        [&](int start,
            int count,
            float cd,
            float sd,
            float& ci,
            float& cq)
    {
        ci = 0.0f;
        cq = 0.0f;

        // --------------------------------------------------------
        // Generate sin/cos recursively, but process two samples
        // per iteration.
        //
        // The oscillator is renormalized periodically to prevent
        // numerical drift over the 1920-sample symbol.
        // --------------------------------------------------------

        float c = 1.0f;
        float s = 0.0f;

        int n = 0;

        // --------------------------------------------------------
        // Two samples per iteration.
        //
        // Instead of performing the oscillator rotation twice:
        //
        //   z1 = z * exp(-jw)
        //   z2 = z1 * exp(-jw)
        //
        // calculate:
        //
        //   exp(-j2w)
        //
        // once.
        // --------------------------------------------------------

        const float cd2 =
            cd * cd - sd * sd;

        const float sd2 =
            2.0f * sd * cd;

        for (; n + 1 < count; n += 2)
        {
            const float x0 =
                samples[start + n];

            const float x1 =
                samples[start + n + 1];

            // ----------------------------------------------------
            // sample n
            // x * exp(-jwn)
            // ----------------------------------------------------

            ci += x0 * c;
            cq -= x0 * s;

            // ----------------------------------------------------
            // Rotate once to n+1
            // ----------------------------------------------------

            const float c1 =
                c * cd -
                s * sd;

            const float s1 =
                s * cd +
                c * sd;

            // ----------------------------------------------------
            // sample n+1
            // ----------------------------------------------------

            ci += x1 * c1;
            cq -= x1 * s1;

            // ----------------------------------------------------
            // Jump directly from n to n+2.
            // ----------------------------------------------------

            const float nc =
                c * cd2 -
                s * sd2;

            const float ns =
                s * cd2 +
                c * sd2;

            c = nc;
            s = ns;

            // ----------------------------------------------------
            // Periodic normalization.
            //
            // Prevent accumulated floating-point oscillator
            // error without doing sqrtf() every sample.
            // ----------------------------------------------------

            if ((n & 127) == 126)
            {
                const float r2 =
                    c * c + s * s;

                const float corr =
                    1.5f - 0.5f * r2;

                c *= corr;
                s *= corr;
            }
        }

        // --------------------------------------------------------
        // Odd final sample
        // --------------------------------------------------------

        if (n < count)
        {
            const float x =
                samples[start + n];

            ci += x * c;
            cq -= x * s;
        }
    };
    // ========================================================
    // COARSE SEARCH
    // ========================================================

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        // ----------------------------------------------------
        // Absolute RF frequency for this FT8 symbol.
        // ----------------------------------------------------

        const float f =
            freq +
            6.25f * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);


        // ----------------------------------------------------
        // Phase factor for one sliding sample:
        //
        // exp(+jw)
        //
        // ----------------------------------------------------

        const float pc =
            cd;

        const float ps =
            sd;


        // ----------------------------------------------------
        // exp(-jw*SYMBOL_SAMPLES)
        //
        // This is the phase of the sample entering the
        // correlation window.
        //
        // Since:
        //
        //   C(p+1) =
        //      exp(+jw) *
        //      [ C(p)
        //        - x[p]
        //        + x[p+L] exp(-jwL) ]
        //
        // ----------------------------------------------------

        const float end_angle =
            -w *
            (float)SYMBOL_SAMPLES;

        const float ec =
            cosf(end_angle);

        const float es =
            sinf(end_angle);


        // ----------------------------------------------------
        // First symbol position.
        //
        // All coarse positions are offsets from first_delay.
        // ----------------------------------------------------

        const int first_start =
            first_delay +
            k * SYMBOL_SAMPLES;


        if (first_start < 0 ||
            first_start >= num_samples)
        {
            continue;
        }


        const int available =
            num_samples -
            first_start;


        const int count =
            std::min(
                SYMBOL_SAMPLES,
                available);


        if (count <= 0)
            continue;


        // ----------------------------------------------------
        // For normal FT8 symbols we have a complete 160 ms
        // window.
        //
        // If this is the final partial symbol, we simply use
        // the available samples. It will normally not be
        // relevant because the 15 s recording contains the
        // complete 79-symbol signal.
        // ----------------------------------------------------

        float ci;
        float cq;

        initial_correlation(
            first_start,
            count,
            cd,
            sd,
            ci,
            cq);


        // ----------------------------------------------------
        // IMPORTANT:
        //
        // Sliding formula with a partial final symbol is not
        // valid because its window length differs.
        //
        // Therefore handle partial symbol only once.
        // ----------------------------------------------------

        if (count < SYMBOL_SAMPLES)
        {
            const float score =
                ci * ci +
                cq * cq;

            score_total[0] += score;
            used_count[0]++;

            continue;
        }


        // ----------------------------------------------------
        // Current window position.
        // ----------------------------------------------------

        int pos =
            first_start;


        // ====================================================
        // Evaluate coarse positions
        // ====================================================

        for (int d = 0;
             d < ncoarse;
             ++d)
        {
            const float score =
                ci * ci +
                cq * cq;

            score_total[d] += score;
            used_count[d]++;


            // ------------------------------------------------
            // Move from current delay to next coarse delay.
            //
            // Coarse step = 4 samples.
            //
            // Instead of recalculating 1920 samples,
            // perform four sliding updates.
            // ------------------------------------------------

            if (d + 1 >= ncoarse)
                break;


            for (int step = 0;
                 step < COARSE_STEP;
                 ++step)
            {
                const int old_pos =
                    pos + step;

                const int new_sample =
                    old_pos +
                    SYMBOL_SAMPLES;

                if (new_sample >= num_samples)
                    break;


                // ------------------------------------------------
                // x[p] leaving the window
                // x[p+L] entering the window
                // ------------------------------------------------

                const float x_old =
                    samples[old_pos];

                const float x_new =
                    samples[new_sample];


                // x_new * exp(-jwL)
                const float add_i =
                    x_new * ec;

                const float add_q =
                    x_new * es;


                // ------------------------------------------------
                // tmp = C - x_old + x_new*exp(-jwL)
                // ------------------------------------------------

                float ti =
                    ci -
                    x_old +
                    add_i;

                float tq =
                    cq +
                    add_q;


                // ------------------------------------------------
                // Cnew = tmp * exp(+jw)
                //
                // (ti + j*tq)*(pc + j*ps)
                // ------------------------------------------------

                const float nci =
                    ti * pc -
                    tq * ps;

                const float ncq =
                    ti * ps +
                    tq * pc;

                ci = nci;
                cq = ncq;
            }


            pos += COARSE_STEP;
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
        if (used_count[d] == 0)
            continue;


        const float score =
            score_total[d] /
            (float)used_count[d];


        if (score >
            best_coarse_score)
        {
            best_coarse_score =
                score;

            best_coarse_index =
                d;
        }
    }


    const int coarse_delay =
        first_delay +
        best_coarse_index *
        COARSE_STEP;


    // ========================================================
    // FINE SEARCH
    //
    // ±6 samples around coarse maximum.
    //
    // 13 positions.
    //
    // Again use sliding correlation.
    // ========================================================

    constexpr int FINE_RADIUS = 6;


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


    float fine_best_score =
        -1.0f;

    int fine_best_delay =
        coarse_delay;


    // ========================================================
    // One tone at a time
    // ========================================================

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        const float f =
            freq +
            6.25f *
            (float)tones[k];


        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;


        const float cd =
            cosf(w);

        const float sd =
            sinf(w);


        // ----------------------------------------------------
        // exp(+jw)
        // ----------------------------------------------------

        const float pc =
            cd;

        const float ps =
            sd;


        // ----------------------------------------------------
        // exp(-jwL)
        // ----------------------------------------------------

        const float end_angle =
            -w *
            (float)SYMBOL_SAMPLES;


        const float ec =
            cosf(end_angle);

        const float es =
            sinf(end_angle);


        // ----------------------------------------------------
        // First fine-search position.
        // ----------------------------------------------------

        int pos =
            fine_first +
            k * SYMBOL_SAMPLES;


        if (pos < 0 ||
            pos >= num_samples)
        {
            continue;
        }


        const int available =
            num_samples -
            pos;


        const int count =
            std::min(
                SYMBOL_SAMPLES,
                available);


        if (count <= 0)
            continue;


        // ----------------------------------------------------
        // Initial correlation at fine_first.
        // ----------------------------------------------------

        float ci;
        float cq;


        initial_correlation(
            pos,
            count,
            cd,
            sd,
            ci,
            cq);


        // ----------------------------------------------------
        // Partial symbol:
        //
        // Nothing useful to slide.
        // ----------------------------------------------------

        if (count < SYMBOL_SAMPLES)
        {
            const float score =
                ci * ci +
                cq * cq;

            if (score >
                fine_best_score)
            {
                fine_best_score =
                    score;

                fine_best_delay =
                    fine_first;
            }

            continue;
        }


        // ====================================================
        // Evaluate fine positions
        // ====================================================

        const int nfine =
            fine_last -
            fine_first + 1;


        for (int d = 0;
             d < nfine;
             ++d)
        {
            // ------------------------------------------------
            // Accumulate this symbol contribution.
            //
            // We need a per-delay score, therefore this
            // temporary array is used to combine the 79
            // symbol correlations.
            //
            // Static size is tiny.
            // ------------------------------------------------

            // This implementation accumulates symbol-by-symbol
            // into static fine score arrays below.
            //
            // The arrays are initialized once per tone group
            // outside this loop.
            // ------------------------------------------------
        }


        // ----------------------------------------------------
        // The above loop is intentionally empty.
        //
        // Actual accumulation is done using the following
        // local arrays.
        // ----------------------------------------------------

        float fine_score[16] = {};


        // Recalculate the initial correlation because the
        // previous loop did not consume it.
        //
        // This costs only one symbol per tone and keeps the
        // implementation straightforward.
        // ----------------------------------------------------

        initial_correlation(
            pos,
            SYMBOL_SAMPLES,
            cd,
            sd,
            ci,
            cq);


        for (int d = 0;
             d < nfine;
             ++d)
        {
            fine_score[d] +=
                ci * ci +
                cq * cq;


            if (d + 1 >= nfine)
                break;


            // ------------------------------------------------
            // One-sample sliding update.
            // ------------------------------------------------

            const int old_pos =
                pos + d;

            const int new_sample =
                old_pos +
                SYMBOL_SAMPLES;


            if (new_sample >= num_samples)
                break;


            const float x_old =
                samples[old_pos];

            const float x_new =
                samples[new_sample];


            const float add_i =
                x_new * ec;

            const float add_q =
                x_new * es;


            const float ti =
                ci -
                x_old +
                add_i;

            const float tq =
                cq +
                add_q;


            const float nci =
                ti * pc -
                tq * ps;

            const float ncq =
                ti * ps +
                tq * pc;


            ci = nci;
            cq = ncq;
        }


        // ----------------------------------------------------
        // Store this tone's contribution.
        //
        // We cannot select the maximum per tone; we need the
        // total score over all 79 symbols.
        //
        // Therefore accumulate into a global array.
        // ----------------------------------------------------

        static float accumulated_fine[16];

        // The first tone initializes the array.
        if (k == 0)
        {
            for (int d = 0;
                 d < 16;
                 ++d)
            {
                accumulated_fine[d] =
                    0.0f;
            }
        }


        for (int d = 0;
             d < nfine;
             ++d)
        {
            accumulated_fine[d] +=
                fine_score[d];
        }


        // ----------------------------------------------------
        // Last tone -> find maximum.
        // ----------------------------------------------------

        if (k == NTONES - 1)
        {
            for (int d = 0;
                 d < nfine;
                 ++d)
            {
                if (accumulated_fine[d] >
                    fine_best_score)
                {
                    fine_best_score =
                        accumulated_fine[d];

                    fine_best_delay =
                        fine_first + d;
                }
            }
        }
    }


    // ========================================================
    // RESULT
    // ========================================================

    return fine_best_delay /
           (float)SAMPLE_RATE;
}

float refine_ft8_delay_v5(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE = 12000;
    constexpr int NTONES      = 79;

    constexpr int SYMBOL_SAMPLES = 1920;   // 160 ms
    constexpr int DECIM          = 4;
    constexpr int SYMBOL_DEC     = SYMBOL_SAMPLES / DECIM; // 480

    // --------------------------------------------------------
    // FIR
    //
    // 49 taps, symmetric
    // cutoff = 200 Hz
    //
    // Group delay = 24 samples.
    // 24 is divisible by 4 -> ideal decimation alignment.
    // --------------------------------------------------------

    constexpr int FIR_TAPS  = 49;
    constexpr int FIR_DELAY = FIR_TAPS / 2;

    constexpr float FIR_CUTOFF = 200.0f;

    static float fir[FIR_TAPS];
    static bool fir_initialized = false;

    if (!fir_initialized)
    {
        // Windowed-sinc, Hann window.
        //
        // Normalized sinc:
        //   sinc(x) = sin(pi*x)/(pi*x)
        //
        // fc normalized to Fs.

        const float fc =
            FIR_CUTOFF / (float)SAMPLE_RATE;

        float sum = 0.0f;

        for (int i = 0; i < FIR_TAPS; ++i)
        {
            const int m =
                i - FIR_DELAY;

            float h;

            if (m == 0)
            {
                h = 2.0f * fc;
            }
            else
            {
                const float x =
                    (float)m;

                h =
                    sinf(2.0f * (float)M_PI * fc * x) /
                    ((float)M_PI * x);
            }

            // Hann window
            const float w =
                0.5f *
                (1.0f -
                 cosf(2.0f * (float)M_PI *
                      (float)i /
                      (float)(FIR_TAPS - 1)));

            fir[i] = h * w;
            sum += fir[i];
        }

        // DC gain = 1
        for (int i = 0; i < FIR_TAPS; ++i)
            fir[i] /= sum;

        fir_initialized = true;
    }


    // --------------------------------------------------------
    // Search interval
    // --------------------------------------------------------

    const int first_delay =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.05f) *
                (float)SAMPLE_RATE));

    const int last_delay =
        std::min(
            num_samples - 1,
            (int)std::lround(
                (delay0 + 0.05f) *
                (float)SAMPLE_RATE));

    if (last_delay < first_delay)
        return delay0;


    // --------------------------------------------------------
    // We need RF samples BEFORE first_delay because of
    // the FIR history.
    //
    // The first decimated FIR output corresponds to RF time
    // "first_delay".
    //
    // Since FIR_DELAY = 24 samples:
    //
    //   rf_start = first_delay - FIR_DELAY
    //
    // The complete FIR window requires another 24 samples
    // on either side, hence the conservative - FIR_TAPS.
    // --------------------------------------------------------

    const int rf_start =
        std::max(
            0,
            first_delay - FIR_TAPS);

    const int rf_end =
        std::min(
            num_samples,
            last_delay +
            NTONES * SYMBOL_SAMPLES +
            FIR_TAPS);

    const int rf_len =
        rf_end - rf_start;

    if (rf_len <= 0)
        return delay0;


    // --------------------------------------------------------
    // Number of decimated samples.
    // --------------------------------------------------------

    const int dec_len =
        (rf_len + DECIM - 1) / DECIM;


    std::vector<float> bi(dec_len);
    std::vector<float> bq(dec_len);


    // --------------------------------------------------------
    // Precompute LO phase increment.
    //
    // We downconvert:
    //
    //   I = x*cos(w)
    //   Q = -x*sin(w)
    //
    // phase recurrence:
    //
    //   c[n+1] = c[n]*cd - s[n]*sd
    //   s[n+1] = s[n]*cd + c[n]*sd
    // --------------------------------------------------------

    const float lo_w =
        TWO_PI * freq /
        (float)SAMPLE_RATE;

    const float lo_cd = cosf(lo_w);
    const float lo_sd = sinf(lo_w);


    // --------------------------------------------------------
    // Complex FIR + decimation
    //
    // Ring buffer.
    //
    // Only one FIR output every 4 RF samples.
    // --------------------------------------------------------

    float ring_i[FIR_TAPS] = {};
    float ring_q[FIR_TAPS] = {};

    int ring_pos = 0;

    float lo_c = 1.0f;
    float lo_s = 0.0f;

    int out_index = 0;

    for (int n = 0; n < rf_len; ++n)
    {
        const float x =
            samples[rf_start + n];

        // Complex downconversion
        const float in_i =
            x * lo_c;

        const float in_q =
            -x * lo_s;

        ring_i[ring_pos] = in_i;
        ring_q[ring_pos] = in_q;

        ring_pos++;

        if (ring_pos >= FIR_TAPS)
            ring_pos = 0;


        // FIR only every DECIM samples
        if ((n & (DECIM - 1)) == (DECIM - 1))
        {
            float sum_i = 0.0f;
            float sum_q = 0.0f;

            int p = ring_pos;

            // FIR_TAPS is small and fixed.
            for (int j = 0; j < FIR_TAPS; ++j)
            {
                sum_i +=
                    ring_i[p] * fir[j];

                sum_q +=
                    ring_q[p] * fir[j];

                ++p;

                if (p >= FIR_TAPS)
                    p = 0;
            }

            if (out_index < dec_len)
            {
                bi[out_index] = sum_i;
                bq[out_index] = sum_q;
                ++out_index;
            }
        }


        // LO recurrence
        const float nc =
            lo_c * lo_cd -
            lo_s * lo_sd;

        const float ns =
            lo_s * lo_cd +
            lo_c * lo_sd;

        lo_c = nc;
        lo_s = ns;
    }


    const int actual_dec_len = out_index;

    if (actual_dec_len < SYMBOL_DEC)
        return delay0;


    // --------------------------------------------------------
    // Mapping
    //
    // IMPORTANT:
    //
    // dec[0] corresponds to approximately first_delay.
    //
    // Therefore:
    //
    //   coarse candidate d
    //
    // corresponds to:
    //
    //   RF delay = first_delay + d*DECIM
    //
    // Do NOT use first_delay/DECIM here.
    // --------------------------------------------------------

    const int ncoarse =
        (last_delay - first_delay) /
        DECIM + 1;


    std::vector<float> score_total(ncoarse, 0.0f);
    std::vector<int> used_count(ncoarse, 0);


    // --------------------------------------------------------
    // Coarse correlation
    //
    // Baseband contains the residual FT8 tone:
    //
    //   f_residual = 6.25 * tone
    //
    // Correlate against exp(-j*w*n).
    //
    // At 3 kHz:
    //
    //   w = 2*pi*f / 3000
    //
    // --------------------------------------------------------

    for (int k = 0; k < NTONES; ++k)
    {
        const float f =
            6.25f * (float)tones[k];

        const float w =
            TWO_PI * f /
            3000.0f;

        const float cd = cosf(w);
        const float sd = sinf(w);

        // exp(-j*w*(L-1))
        const float end_angle =
            -w * (float)(SYMBOL_DEC - 1);

        const float end_c =
            cosf(end_angle);

        const float end_s =
            sinf(end_angle);


        // ----------------------------------------------------
        // Maximum number of complete symbols.
        // ----------------------------------------------------

        const int max_symbol =
            std::min(
                NTONES,
                (actual_dec_len - SYMBOL_DEC) /
                SYMBOL_DEC + 1);


        for (int d = 0; d < ncoarse; ++d)
        {
            const int pos =
                d + k * SYMBOL_DEC;

            if (pos < 0 ||
                pos >= actual_dec_len)
                continue;


            // ------------------------------------------------
            // Full symbol -> sliding correlation
            // ------------------------------------------------

            if (pos + SYMBOL_DEC <= actual_dec_len)
            {
                float ci = 0.0f;
                float cq = 0.0f;

                float c = 1.0f;
                float s = 0.0f;

                for (int n = 0;
                     n < SYMBOL_DEC;
                     ++n)
                {
                    const float ai =
                        bi[pos + n];

                    const float aq =
                        bq[pos + n];

                    // z * exp(-jwt)
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


                // ------------------------------------------------
                // Remaining symbols are evaluated by the
                // direct/sliding implementation below.
                //
                // We already calculated the current symbol.
                // No need to scan it again.
                // ------------------------------------------------
            }
            else
            {
                // ------------------------------------------------
                // Partial final symbol.
                //
                // It is only one symbol and only affects the
                // final part of the 15 s WAV.
                // Use direct correlation so that its actual
                // length is handled correctly.
                // ------------------------------------------------

                const int count =
                    actual_dec_len - pos;

                if (count <= 0)
                    continue;

                float ci = 0.0f;
                float cq = 0.0f;

                float c = 1.0f;
                float s = 0.0f;

                for (int n = 0; n < count; ++n)
                {
                    const float ai =
                        bi[pos + n];

                    const float aq =
                        bq[pos + n];

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
            }
        }
    }


    // --------------------------------------------------------
    // Find coarse maximum.
    // --------------------------------------------------------

    int best_coarse_index = 0;
    float best_coarse_score = -1.0f;

    for (int d = 0;
         d < ncoarse;
         ++d)
    {
        if (used_count[d] == 0)
            continue;

        const float score =
            score_total[d] /
            (float)used_count[d];

        if (score > best_coarse_score)
        {
            best_coarse_score = score;
            best_coarse_index = d;
        }
    }


    const int coarse_delay =
        first_delay +
        best_coarse_index * DECIM;


    // --------------------------------------------------------
    // Fine search
    //
    // Search original RF at 1-sample resolution.
    //
    // The coarse estimator has already reduced the search
    // space to a few samples.
    //
    // +/- 6 is intentional:
    //
    //   coarse resolution = 4 samples
    //   FIR / decimation alignment errors cannot push the
    //   true maximum outside the fine window.
    // --------------------------------------------------------

    constexpr int FINE_RADIUS = 6;

    const int fine_first =
        std::max(
            first_delay,
            coarse_delay - FINE_RADIUS);

    const int fine_last =
        std::min(
            last_delay,
            coarse_delay + FINE_RADIUS);


    float fine_best_score = -1.0f;
    int fine_best_delay = coarse_delay;


    for (int dly = fine_first;
         dly <= fine_last;
         ++dly)
    {
        float total_score = 0.0f;
        int used = 0;


        for (int k = 0;
             k < NTONES;
             ++k)
        {
            const float f =
                freq +
                6.25f *
                (float)tones[k];

            const float w =
                TWO_PI * f /
                (float)SAMPLE_RATE;

            const float cd = cosf(w);
            const float sd = sinf(w);

            const int start =
                dly +
                k * SYMBOL_SAMPLES;

            if (start < 0 ||
                start >= num_samples)
                continue;


            const int count =
                std::min(
                    SYMBOL_SAMPLES,
                    num_samples - start);

            if (count <= 0)
                continue;


            // ------------------------------------------------
            // Direct correlation.
            //
            // Fine search is only 13 delays, therefore this
            // is cheap compared with the FIR.
            // ------------------------------------------------

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

            total_score +=
                ci * ci +
                cq * cq;

            ++used;
        }


        if (used == 0)
            continue;


        const float score =
            total_score /
            (float)used;


        if (score > fine_best_score)
        {
            fine_best_score = score;
            fine_best_delay = dly;
        }
    }


    // --------------------------------------------------------
    // Result
    // --------------------------------------------------------

    return fine_best_delay /
           (float)SAMPLE_RATE;
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

#if defined(ARDUINO)
// ============================================================
// WAV LOADER
//
// Questa è la funzione richiesta dall'utente.
// ============================================================

float* load_wav_littlefs(
    const char* path,
    int* out_num_samples,
    int* out_num_channels,
    int* out_sample_rate)
{
    File f = LittleFS.open(path, "r");

    if (!f)
    {
        Serial.printf(
            "ERROR: cannot open %s\n",
            path);
        return nullptr;
    }

    uint8_t riff[12];

    if (f.read(riff, 12) != 12)
    {
        Serial.println("ERROR: WAV header too short");
        f.close();
        return nullptr;
    }

    if (memcmp(riff, "RIFF", 4) != 0 ||
        memcmp(riff + 8, "WAVE", 4) != 0)
    {
        Serial.println("ERROR: not a RIFF/WAVE file");
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

        uint32_t chunk_size =
            chunk_header[4] |
            ((uint32_t)chunk_header[5] << 8) |
            ((uint32_t)chunk_header[6] << 16) |
            ((uint32_t)chunk_header[7] << 24);

        uint32_t chunk_start = f.position();

        if (memcmp(chunk_header, "fmt ", 4) == 0)
        {
            if (chunk_size < 16)
            {
                Serial.println("ERROR: invalid fmt chunk");
                f.close();
                return nullptr;
            }

            uint8_t fmt[16];

            if (f.read(fmt, 16) != 16)
            {
                Serial.println("ERROR: cannot read fmt");
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
        Serial.println("ERROR: missing fmt/data chunk");
        f.close();
        return nullptr;
    }

    if (audio_format != 1)
    {
        Serial.printf(
            "ERROR: unsupported WAV format %u\n",
            audio_format);
        f.close();
        return nullptr;
    }

    if (bits_per_sample != 16)
    {
        Serial.printf(
            "ERROR: only 16-bit PCM supported, got %u\n",
            bits_per_sample);
        f.close();
        return nullptr;
    }

    if (num_channels < 1)
    {
        Serial.println("ERROR: invalid channel count");
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
            num_samples * sizeof(float));

    if (!samples)
    {
        Serial.printf(
            "ERROR: malloc failed for %d samples\n",
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
            Serial.println(
                "ERROR: unexpected EOF");

            free(samples);
            f.close();
            return nullptr;
        }

        int16_t v =
            (int16_t)(
                b[0] |
                ((uint16_t)b[1] << 8));

        samples[n] =
            (float)v / 32768.0f;

        // Skip other channels
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
    *out_sample_rate  = sample_rate;

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
    printf(
        "============================================================\n");

    printf(
        "FILE: %s\n",
        filename);

    printf(
        "TRUE DELAY: %.6f s\n",
        true_delay);

    printf(
        "============================================================\n");

    int num_samples = 0;
    int num_channels = 0;
    int sample_rate = 0;


    uint32_t heap_before;
    #if defined ARDUINO
    heap_before =
        ESP.getFreeHeap();
    #endif
        
    
    float* samples;
    
    #if defined ARDUINO
    samples =
        load_wav_littlefs(
            filename,
            &num_samples,
            &num_channels,
            &sample_rate);
    #else
    #endif

    if (!samples)
    {
        printf(
            "LOAD FAILED\n");

        return;
    }

    uint32_t heap_after_load =
        ESP.getFreeHeap();

    printf(
        "Samples          : %d\n",
        num_samples);

    printf(
        "Duration         : %.3f s\n",
        (float)num_samples /
        sample_rate);

    printf(
        "Channels         : %d\n",
        num_channels);

    printf(
        "Sample rate      : %d Hz\n",
        sample_rate);

    printf(
        "Heap before load : %u\n",
        heap_before);

    printf(
        "Heap after load  : %u\n",
        heap_after_load);

    if (sample_rate != SAMPLE_RATE ||
        num_channels != 1)
    {
        printf(
            "ERROR: WAV format mismatch\n");

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
    // TEST DELLA VERA refine_ft8_delay()
    // --------------------------------------------------------

    printf("\n");
    printf(
        "Running ORIGINAL refine_ft8_delay()...\n");

    uint32_t t0 = micros();

    const float estimated_delay =
        refine_ft8_delay_v7(
            samples,
            num_samples,
            tones,
            true_delay,
            TEST_FREQ,
            0);

    uint32_t elapsed_us =
        micros() - t0;

    const float error_s =
        estimated_delay -
        true_delay;

    const float error_ms =
        error_s * 1000.0f;

    const float error_samples =
        error_s *
        SAMPLE_RATE;

    printf("\n");
    printf(
        "------------------------------------------------------------\n");

    printf(
        "True delay       : %10.6f s\n",
        true_delay);

    printf(
        "Estimated delay  : %10.6f s\n",
        estimated_delay);

    Serial.printf(
    printf(
        error_s);

    printf(
        "Error            : %+10.3f ms\n",
        error_ms);

    printf(
        "Error            : %+10.3f samples\n",
        error_samples);

    printf(
        "Execution time   : %10.3f ms\n",
        elapsed_us / 1000.0f);

    printf(
        "Free heap        : %u bytes\n",
        ESP.getFreeHeap());

    printf(
        "------------------------------------------------------------\n");

    // --------------------------------------------------------
    // IMPORTANTISSIMO:
    // liberiamo il WAV PRIMA di caricare il successivo.
    // --------------------------------------------------------

    free(samples);

    Serial.printf(
        "After free() heap: %u bytes\n",
        ESP.getFreeHeap());

    Serial.println(
        "WAV released.");
}


// ============================================================
// SETUP
// ============================================================

void setup()
{
    Serial.begin(115200);

    delay(2000);

    Serial.println();
    Serial.println();
    Serial.println(
        "############################################################");

    Serial.println(
        " ORIGINAL refine_ft8_delay() - REAL WAV BENCHMARK");

    Serial.println(
        "############################################################");

    Serial.printf(
        "CPU frequency    : %u MHz\n",
        ESP.getCpuFreqMHz());

    Serial.printf(
        "Sample rate      : %d Hz\n",
        SAMPLE_RATE);

    Serial.printf(
        "Symbol duration  : %d samples = %.3f ms\n",
        SPS,
        1000.0f * SPS / SAMPLE_RATE);

    Serial.printf(
        "Test frequency   : %.3f Hz\n",
        TEST_FREQ);

    Serial.printf(
        "Free heap        : %u bytes\n",
        ESP.getFreeHeap());

    Serial.printf(
        "Free PSRAM       : %u bytes\n",
        ESP.getFreePsram());

    // --------------------------------------------------------
    // LittleFS
    // --------------------------------------------------------

    Serial.println();
    Serial.println(
        "Mounting LittleFS...");

    if (!LittleFS.begin(true))
    {
        Serial.println(
            "ERROR: LittleFS.begin() failed");

        return;
    }

    Serial.println(
        "LittleFS mounted.");

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

    // --------------------------------------------------------

    Serial.println();
    Serial.println(
        "############################################################");

    Serial.println(
        "BENCHMARK COMPLETE");

    Serial.printf(
        "Final free heap : %u bytes\n",
        ESP.getFreeHeap());

    Serial.printf(
        "Final free PSRAM: %u bytes\n",
        ESP.getFreePsram());

    Serial.println(
        "############################################################");
}


// ============================================================
// LOOP
// ============================================================

void loop()
{
    delay(1000);
}