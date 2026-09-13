// ============================================================
// test_refine_ft8_delay_speed.cpp
//
// V3G FT8 DELAY ESTIMATOR BENCHMARK
//
// V3G ISOLATION TEST:
//
//   Four independent correlation tests:
//
//       RF      + reference Q = +sin
//       RF      + reference Q = -sin
//
//       IdealBB + reference Q = +sin
//       IdealBB + reference Q = -sin
//
// Purpose:
//   Determine whether the systematic delay bias is caused by
//   the RF mixer / Q sign convention, or already exists in the
//   common BB / FIR / reference / correlation geometry.
//
// Main diagnostic:
//
//       RF   +sin
//       RF   -sin
//       BB   +sin
//       BB   -sin
//
// If the famous approximately -12.466 sample bias disappears
// when Q=-sin, this is strong evidence of a Q-sign / phase
// convention problem.
//
// Compatible with:
//   - Linux
//   - ESP32-S3 / ESP-IDF / Arduino
//
// ============================================================

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <assert.h>

#ifdef ARDUINO
#include <Arduino.h>
#endif


// ============================================================
// CONFIGURATION
// ============================================================

static constexpr int SAMPLE_RATE = 12000;

static constexpr int NSYM = 79;
static constexpr int SPS  = 1920;       // 160 ms

static constexpr float TEST_FREQ  = 1001.0f;
static constexpr float TEST_DELAY = 2.000f;
static constexpr float AMPLITUDE  = 1000.0f;

static constexpr int NWARMUP = 10;
static constexpr int NITER   = 100;

static constexpr float TONE_SPACING = 6.25f;


// ============================================================
// DELAY SWEEP
// ============================================================

static constexpr int NDELAYS = 8;

static const float TEST_DELAYS[NDELAYS] =
{
    0.5f,
    1.0f,
    1.5f,
    2.0f,
    2.5f,
    3.0f,
    3.5f,
    4.0f
};


// ============================================================
// FT8 COSTAS
// ============================================================

static const int COSTAS[7] =
{
    3, 1, 4, 0, 6, 5, 2
};


// ============================================================
// FIR / CORRELATION GEOMETRY
// ============================================================

static constexpr int FIR_TAPS    = 129;
static constexpr int FIR_DELAY   = 64;
static constexpr float LPF_CUTOFF = 150.0f;

static constexpr int WIN    = 128;
static constexpr int SEARCH = 15;

static constexpr int NTRANS = 22;


// Search geometry:
//
// q = local output sample coordinate
//
// Correlation window:
//     q ... q + WIN - 1
//
// Search:
//     [-SEARCH ... +SEARCH]
//

static constexpr int QMIN =
    -SEARCH - WIN / 2;

static constexpr int QMAX =
     SEARCH + WIN / 2 - 1;

static constexpr int INPUT_MIN =
    QMIN - FIR_DELAY;

static constexpr int INPUT_MAX =
    QMAX + FIR_DELAY;

static constexpr int NINPUT =
    INPUT_MAX - INPUT_MIN + 1;

static constexpr int NOUTPUT =
    NINPUT - FIR_TAPS + 1;


// ============================================================
// STATIC WORK BUFFERS
// ============================================================

static float local_mix_i[NINPUT];
static float local_mix_q[NINPUT];

static float local_bb_i[NOUTPUT];
static float local_bb_q[NOUTPUT];


// ============================================================
// FIR
// ============================================================

static float fir[FIR_TAPS];


// ============================================================
// REFERENCE CACHE
// ============================================================

static float ref_i[NTRANS * NOUTPUT];
static float ref_q[NTRANS * NOUTPUT];

static float ref_power_cache[NTRANS];

static int cached_before[NTRANS];
static int cached_after[NTRANS];

static bool ref_initialized = false;


// ============================================================
// TRANSITION INFORMATION
// ============================================================

struct TransitionInfo
{
    int sample;
    int before_tone;
    int after_tone;
};

static TransitionInfo transitions[NTRANS];


// ============================================================
// DEBUG INFORMATION
// ============================================================

struct RefineDebug
{
    int best_offset;
    float best_score;
    float refined_offset;
};

static RefineDebug g_debug;


// ============================================================
// FT8 TEST TONES
// ============================================================

static void make_test_tones(uint8_t *tones)
{
    uint32_t state = 0x12345678u;

    for (int i = 0; i < NSYM; ++i)
    {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;

        tones[i] = (uint8_t)(state % 8);
    }

    // First Costas
    for (int i = 0; i < 7; ++i)
        tones[i] = COSTAS[i];

    // Middle Costas
    for (int i = 0; i < 7; ++i)
        tones[36 + i] = COSTAS[i];

    // Last Costas
    for (int i = 0; i < 7; ++i)
        tones[72 + i] = COSTAS[i];
}


// ============================================================
// SYNTHETIC FT8 RF GENERATOR
// ============================================================

static void generate_ft8_continuous(
    std::vector<float> &signal,
    float delay_sec,
    float carrier_freq,
    float amplitude)
{
    uint8_t tones[NSYM];

    make_test_tones(tones);

    const int delay_samples =
        (int)lroundf(delay_sec * SAMPLE_RATE);

    const int total_samples =
        delay_samples + NSYM * SPS;

    signal.assign(total_samples, 0.0f);

    double phase = 0.0;

    const double two_pi =
        2.0 * M_PI;

    for (int n = 0; n < total_samples; ++n)
    {
        if (n < delay_samples)
        {
            signal[n] = 0.0f;
            continue;
        }

        const int k =
            n - delay_samples;

        const int sym =
            k / SPS;

        if (sym >= NSYM)
        {
            signal[n] = 0.0f;
            continue;
        }

        const float tone =
            (float)tones[sym];

        const float freq =
            carrier_freq +
            tone * TONE_SPACING;

        phase +=
            two_pi * (double)freq /
            (double)SAMPLE_RATE;

        if (phase > two_pi)
            phase -= two_pi;

        signal[n] =
            amplitude *
            cosf((float)phase);
    }
}


// ============================================================
// BUILD COSTAS TRANSITIONS
// ============================================================

static int build_transitions(
    TransitionInfo *out,
    int max_transitions)
{
    uint8_t tones[NSYM];

    make_test_tones(tones);

    int count = 0;

    const int starts[3] =
    {
        0,
        36,
        72
    };

    for (int b = 0; b < 3; ++b)
    {
        const int start = starts[b];

        for (int i = 0; i < 6; ++i)
        {
            if (count >= max_transitions)
                return count;

            const int before =
                tones[start + i];

            const int after =
                tones[start + i + 1];

            out[count].sample =
                (start + i + 1) * SPS;

            out[count].before_tone = before;
            out[count].after_tone  = after;

            ++count;
        }
    }

    //
    // Four additional transitions around the Costas blocks.
    //

    const int extra[4] =
    {
        7,
        35,
        43,
        71
    };

    for (int i = 0; i < 4 && count < max_transitions; ++i)
    {
        const int s = extra[i];

        if (s > 0 && s < NSYM)
        {
            out[count].sample =
                s * SPS;

            out[count].before_tone =
                tones[s - 1];

            out[count].after_tone =
                tones[s];

            ++count;
        }
    }

    return count;
}


// ============================================================
// INITIALIZE FIR
// ============================================================

static void init_fir()
{
    const int M = FIR_TAPS - 1;

    const float fc =
        LPF_CUTOFF /
        (float)SAMPLE_RATE;

    double sum = 0.0;

    for (int n = 0; n < FIR_TAPS; ++n)
    {
        const int k =
            n - M / 2;

        double h;

        if (k == 0)
        {
            h = 2.0 * fc;
        }
        else
        {
            const double x =
                M_PI *
                2.0 *
                fc *
                (double)k;

            h =
                sin(x) /
                (M_PI * (double)k);
        }

        //
        // Hann window.
        //

        const double w =
            0.5 -
            0.5 *
            cos(
                2.0 * M_PI *
                (double)n /
                (double)M);

        h *= w;

        fir[n] =
            (float)h;

        sum += h;
    }

    //
    // Normalize DC gain to 1.
    //

    for (int n = 0; n < FIR_TAPS; ++n)
        fir[n] /= (float)sum;
}


// ============================================================
// BUILD REFERENCE FOR ONE TRANSITION
// ============================================================

static void build_reference_for_transition(
    int index,
    int before_tone,
    int after_tone)
{
    const float step_before =
        2.0f * M_PI *
        (
            TEST_FREQ +
            before_tone * TONE_SPACING
        ) /
        (float)SAMPLE_RATE;

    const float step_after =
        2.0f * M_PI *
        (
            TEST_FREQ +
            after_tone * TONE_SPACING
        ) /
        (float)SAMPLE_RATE;

    static float input_i[NINPUT];
    static float input_q[NINPUT];

    for (int n = 0; n < NINPUT; ++n)
    {
        const int relative =
            INPUT_MIN + n;

        float phase = 0.0f;

        if (relative < 0)
        {
            phase =
                (float)relative *
                step_before;
        }
        else
        {
            phase =
                (float)relative *
                step_after;
        }

        input_i[n] =
            cosf(phase);

        input_q[n] =
            sinf(phase);
    }

    float power = 0.0f;

    for (int n = 0; n < NOUTPUT; ++n)
    {
        float yi = 0.0f;
        float yq = 0.0f;

        for (int k = 0; k < FIR_TAPS; ++k)
        {
            yi +=
                input_i[n + k] *
                fir[k];

            yq +=
                input_q[n + k] *
                fir[k];
        }

        ref_i[index * NOUTPUT + n] =
            yi;

        ref_q[index * NOUTPUT + n] =
            yq;

        power +=
            yi * yi +
            yq * yq;
    }

    ref_power_cache[index] =
        power;

    cached_before[index] =
        before_tone;

    cached_after[index] =
        after_tone;
}


// ============================================================
// ENSURE REFERENCE CACHE
// ============================================================

static void ensure_reference_cache()
{
    if (ref_initialized)
        return;

    for (int i = 0; i < NTRANS; ++i)
    {
        build_reference_for_transition(
            i,
            transitions[i].before_tone,
            transitions[i].after_tone);
    }

    ref_initialized = true;
}


// ============================================================
// PREPARE LOCAL RF -> COMPLEX BASEBAND
// ============================================================

static void prepare_local_baseband(
    const std::vector<float> &rf,
    int transition_sample,
    int before_tone,
    int after_tone)
{
    const float step_before =
        2.0f * M_PI *
        (
            TEST_FREQ +
            before_tone * TONE_SPACING
        ) /
        (float)SAMPLE_RATE;

    const float step_after =
        2.0f * M_PI *
        (
            TEST_FREQ +
            after_tone * TONE_SPACING
        ) /
        (float)SAMPLE_RATE;

    for (int n = 0; n < NINPUT; ++n)
    {
        const int relative =
            INPUT_MIN + n;

        const int sample =
            transition_sample +
            relative;

        float x = 0.0f;

        if (sample >= 0 &&
            sample < (int)rf.size())
        {
            x = rf[sample];
        }

        float phase = 0.0f;

        if (relative < 0)
        {
            phase =
                (float)relative *
                step_before;
        }
        else
        {
            phase =
                (float)relative *
                step_after;
        }

        const float c =
            cosf(phase);

        const float s =
            sinf(phase);

        //
        // Real RF mixer:
        //
        //     I = RF * cos(fc)
        //     Q = RF * (-sin(fc))
        //
        // Desired BB:
        //
        //     I = +0.5 cos(fb)
        //     Q = -0.5 sin(fb)
        //

        local_mix_i[n] =
            x * c;

        local_mix_q[n] =
            -x * s;
    }

    //
    // FIR.
    //

    for (int n = 0; n < NOUTPUT; ++n)
    {
        float yi = 0.0f;
        float yq = 0.0f;

        for (int k = 0; k < FIR_TAPS; ++k)
        {
            yi +=
                local_mix_i[n + k] *
                fir[k];

            yq +=
                local_mix_q[n + k] *
                fir[k];
        }

        local_bb_i[n] = yi;
        local_bb_q[n] = yq;
    }
}


// ============================================================
// PREPARE IDEAL COMPLEX BASEBAND
// ============================================================
//
// This is deliberately the complex signal corresponding to
// the output of the REAL RF mixer:
//
//     I = +0.5 cos(phi)
//     Q = -0.5 sin(phi)
//
// No RF signal is generated and no RF mixing is performed.
//
// ============================================================

static void prepare_ideal_baseband(
    int before_tone,
    int after_tone)
{
    const float step_before =
        2.0f * M_PI *
        (
            TEST_FREQ +
            before_tone * TONE_SPACING
        ) /
        (float)SAMPLE_RATE;

    const float step_after =
        2.0f * M_PI *
        (
            TEST_FREQ +
            after_tone * TONE_SPACING
        ) /
        (float)SAMPLE_RATE;

    static float ideal_i[NINPUT];
    static float ideal_q[NINPUT];

    for (int n = 0; n < NINPUT; ++n)
    {
        const int relative =
            INPUT_MIN + n;

        float phase = 0.0f;

        if (relative < 0)
        {
            phase =
                (float)relative *
                step_before;
        }
        else
        {
            phase =
                (float)relative *
                step_after;
        }

        //
        // Exact complex output convention of the real mixer.
        //

        ideal_i[n] =
            0.5f *
            cosf(phase);

        ideal_q[n] =
            -0.5f *
            sinf(phase);
    }

    //
    // Same FIR as RF path.
    //

    for (int n = 0; n < NOUTPUT; ++n)
    {
        float yi = 0.0f;
        float yq = 0.0f;

        for (int k = 0; k < FIR_TAPS; ++k)
        {
            yi +=
                ideal_i[n + k] *
                fir[k];

            yq +=
                ideal_q[n + k] *
                fir[k];
        }

        local_bb_i[n] =
            yi;

        local_bb_q[n] =
            yq;
    }
}


// ============================================================
// NORMALIZED COMPLEX CORRELATION
// ============================================================
//
// reference_q_sign:
//
//     +1 -> reference Q = +sin
//     -1 -> reference Q = -sin
//
// The reference itself is cached with Q=+sin.
//
// For Q=-sin we simply invert the cached reference Q.
//
// This is intentional: the only thing being changed is the
// reference Q sign.
//
// ============================================================

static float correlation_score(
    int transition_index,
    int offset,
    int reference_q_sign)
{
    const float *ri =
        &ref_i[
            transition_index *
            NOUTPUT
        ];

    const float *rq =
        &ref_q[
            transition_index *
            NOUTPUT
        ];

    const int start =
        offset + SEARCH;

    if (start < 0 ||
        start + WIN > NOUTPUT)
    {
        return 0.0f;
    }

    float corr_i = 0.0f;
    float corr_q = 0.0f;

    float local_power = 0.0f;

    for (int n = 0; n < WIN; ++n)
    {
        const int p =
            start + n;

        const float ai =
            local_bb_i[p];

        const float aq =
            local_bb_q[p];

        const float bi =
            ri[p];

        //
        // Explicit Q sign selection.
        //

        const float bq =
            (float)reference_q_sign *
            rq[p];

        //
        // Complex correlation:
        //
        //     A * conj(B)
        //
        // real:
        //     ai*bi + aq*bq
        //
        // imag:
        //     aq*bi - ai*bq
        //

        corr_i +=
            ai * bi +
            aq * bq;

        corr_q +=
            aq * bi -
            ai * bq;

        local_power +=
            ai * ai +
            aq * aq;
    }

    const float ref_power =
        ref_power_cache[
            transition_index
        ];

    const float denom =
        sqrtf(
            local_power *
            ref_power);

    if (denom <= 1.0e-20f)
        return 0.0f;

    return
        (
            corr_i * corr_i +
            corr_q * corr_q
        ) /
        (
            denom * denom
        );
}


// ============================================================
// DELAY REFINEMENT
// ============================================================

static float refine_ft8_delay(
    const std::vector<float> &rf,
    float initial_delay)
{
    int global_best_offset = 0;

    float global_best_score = -1.0f;

    for (int t = 0; t < NTRANS; ++t)
    {
        const int transition_sample =
            transitions[t].sample +
            (int)lroundf(
                initial_delay *
                SAMPLE_RATE);

        prepare_local_baseband(
            rf,
            transition_sample,
            transitions[t].before_tone,
            transitions[t].after_tone);

        for (int offset = -SEARCH;
             offset <= SEARCH;
             ++offset)
        {
            const float score =
                correlation_score(
                    t,
                    offset,
                    +1);

            if (score > global_best_score)
            {
                global_best_score =
                    score;

                global_best_offset =
                    offset;
            }
        }
    }

    int best_transition = -1;

    float best_transition_score =
        -1.0f;

    for (int t = 0; t < NTRANS; ++t)
    {
        const int transition_sample =
            transitions[t].sample +
            (int)lroundf(
                initial_delay *
                SAMPLE_RATE);

        prepare_local_baseband(
            rf,
            transition_sample,
            transitions[t].before_tone,
            transitions[t].after_tone);

        const float score =
            correlation_score(
                t,
                global_best_offset,
                +1);

        if (score >
            best_transition_score)
        {
            best_transition_score =
                score;

            best_transition =
                t;
        }
    }

    float refined_offset =
        (float)global_best_offset;

    if (best_transition >= 0 &&
        global_best_offset > -SEARCH &&
        global_best_offset < SEARCH)
    {
        const int transition_sample =
            transitions[best_transition].sample +
            (int)lroundf(
                initial_delay *
                SAMPLE_RATE);

        prepare_local_baseband(
            rf,
            transition_sample,
            transitions[best_transition].before_tone,
            transitions[best_transition].after_tone);

        const float ym =
            correlation_score(
                best_transition,
                global_best_offset - 1,
                +1);

        const float y0 =
            correlation_score(
                best_transition,
                global_best_offset,
                +1);

        const float yp =
            correlation_score(
                best_transition,
                global_best_offset + 1,
                +1);

        const float denom =
            ym -
            2.0f * y0 +
            yp;

        if (fabsf(denom) >
            1.0e-12f)
        {
            float delta =
                0.5f *
                (ym - yp) /
                denom;

            if (delta > 1.0f)
                delta = 1.0f;

            if (delta < -1.0f)
                delta = -1.0f;

            refined_offset +=
                delta;
        }
    }

    g_debug.best_offset =
        global_best_offset;

    g_debug.best_score =
        global_best_score;

    g_debug.refined_offset =
        refined_offset;

    return
        initial_delay +
        refined_offset /
        (float)SAMPLE_RATE;
}


// ============================================================
// V3e DIAGNOSTIC
// ============================================================

static void print_v3e_diagnostic()
{
    printf("\n");
    printf("============================================================\n");
    printf("V3e GEOMETRY DIAGNOSTIC\n");
    printf("============================================================\n");

    printf("Sample rate       : %d Hz\n",
           SAMPLE_RATE);

    printf("Symbol duration   : %d samples = %.3f ms\n",
           SPS,
           1000.0f * SPS / SAMPLE_RATE);

    printf("FIR taps          : %d\n",
           FIR_TAPS);

    printf("FIR delay         : %d samples = %.3f ms\n",
           FIR_DELAY,
           1000.0f * FIR_DELAY / SAMPLE_RATE);

    printf("LPF cutoff        : %.1f Hz\n",
           LPF_CUTOFF);

    printf("Correlation WIN   : %d samples = %.3f ms\n",
           WIN,
           1000.0f * WIN / SAMPLE_RATE);

    printf("Search            : +/- %d samples\n",
           SEARCH);

    printf("Transitions       : %d\n",
           NTRANS);

    printf("INPUT_MIN         : %d\n",
           INPUT_MIN);

    printf("INPUT_MAX         : %d\n",
           INPUT_MAX);

    printf("NINPUT            : %d\n",
           NINPUT);

    printf("NOUTPUT           : %d\n",
           NOUTPUT);

    printf("============================================================\n");

    std::vector<float> signal;

    generate_ft8_continuous(
        signal,
        TEST_DELAY,
        TEST_FREQ,
        AMPLITUDE);

    const float estimated =
        refine_ft8_delay(
            signal,
            TEST_DELAY);

    printf("True delay        : %.6f ms\n",
           TEST_DELAY * 1000.0f);

    printf("Best offset       : %d samples\n",
           g_debug.best_offset);

    printf("Best score        : %.9f\n",
           g_debug.best_score);

    printf("Refined offset    : %.6f samples\n",
           g_debug.refined_offset);

    printf("Estimated delay   : %.6f ms\n",
           estimated * 1000.0f);

    printf("Residual error    : %.6f samples\n",
           (estimated - TEST_DELAY) *
           SAMPLE_RATE);

    printf("Residual error    : %.6f ms\n",
           (estimated - TEST_DELAY) *
           1000.0f);

    printf("============================================================\n");
}


// ============================================================
// DELAY SWEEP
// ============================================================

static void run_delay_sweep()
{
    printf("\n");
    printf("============================================================\n");
    printf("V3e DELAY SWEEP\n");
    printf("============================================================\n");

    printf(
        "%10s %16s %16s %16s %16s\n",
        "True ms",
        "Best samples",
        "Refined samples",
        "Estimated ms",
        "Error samples");

    printf(
        "--------------------------------------------------------------------------\n");

    for (int i = 0; i < NDELAYS; ++i)
    {
        const float true_delay =
            TEST_DELAYS[i];

        std::vector<float> signal;

        generate_ft8_continuous(
            signal,
            true_delay,
            TEST_FREQ,
            AMPLITUDE);

        const float estimated =
            refine_ft8_delay(
                signal,
                true_delay);

        const float error_samples =
            (estimated - true_delay) *
            SAMPLE_RATE;

        printf(
            "%10.3f %16d %16.6f %16.6f %16.6f\n",
            true_delay * 1000.0f,
            g_debug.best_offset,
            g_debug.refined_offset,
            estimated * 1000.0f,
            error_samples);
    }

    printf(
        "==========================================================================\n");
}


// ============================================================
// V3G FOUR-WAY ISOLATION TEST
// ============================================================
//
// Four independent experiments:
//
//   1. RF      + Q=+sin
//   2. RF      + Q=-sin
//   3. IdealBB + Q=+sin
//   4. IdealBB + Q=-sin
//
// IMPORTANT:
//
// The ONLY variable in the correlation is the sign of the
// reference Q component.
//
// The RF path contains:
//
//     synthetic RF
//         |
//         v
//     real mixer
//         |
//         v
//        FIR
//         |
//         v
//     correlation
//
// The Ideal BB path contains:
//
//     ideal complex BB
//         |
//         v
//        FIR
//         |
//         v
//     correlation
//
// Therefore:
//
// RF vs BB isolates the RF mixer.
//
// +sin vs -sin isolates the reference-Q convention.
//
// ============================================================

struct V3GResult
{
    int best_offset;
    float best_score;
};


static V3GResult run_one_correlation_test(
    int transition_index,
    int reference_q_sign)
{
    V3GResult result;

    result.best_offset = 0;
    result.best_score  = -1.0f;

    for (int offset = -SEARCH;
         offset <= SEARCH;
         ++offset)
    {
        const float score =
            correlation_score(
                transition_index,
                offset,
                reference_q_sign);

        if (score > result.best_score)
        {
            result.best_score =
                score;

            result.best_offset =
                offset;
        }
    }

    return result;
}


// ============================================================
// PRINT CORRELATION CURVE
// ============================================================

static void print_correlation_curve(
    int transition_index,
    int reference_q_sign)
{
    printf(
        "      offset        score\n");

    printf(
        "      --------------------------\n");

    for (int offset = -SEARCH;
         offset <= SEARCH;
         ++offset)
    {
        const float score =
            correlation_score(
                transition_index,
                offset,
                reference_q_sign);

        printf(
            "      %+8d    %.9f\n",
            offset,
            score);
    }
}


// ============================================================
// V3G ISOLATION TEST
// ============================================================

static void run_v3g_isolation_test()
{
    printf("\n");
    printf("============================================================\n");
    printf("V3G FT8 DELAY ESTIMATOR FOUR-WAY ISOLATION TEST\n");
    printf("============================================================\n");

    printf("Purpose:\n");
    printf("  Test RF / BB independently with reference Q +/- sin\n");

    printf("\n");

    printf("Four tests:\n");
    printf("  A: RF      + Q = +sin\n");
    printf("  B: RF      + Q = -sin\n");
    printf("  C: IdealBB + Q = +sin\n");
    printf("  D: IdealBB + Q = -sin\n");

    printf("============================================================\n");

    //
    // Infrastructure.
    //

    init_fir();

    const int ntrans =
        build_transitions(
            transitions,
            NTRANS);

    if (ntrans != NTRANS)
    {
        printf(
            "ERROR: expected %d transitions, got %d\n",
            NTRANS,
            ntrans);

        return;
    }

    ensure_reference_cache();

    //
    // Fixed transition.
    //

    const int TEST_TRANSITION =
        3;

    const float TEST_DELAY_LOCAL =
        2.000f;

    printf("Sample rate       : %d Hz\n",
           SAMPLE_RATE);

    printf("Test frequency    : %.3f Hz\n",
           TEST_FREQ);

    printf("Test delay        : %.3f ms\n",
           TEST_DELAY_LOCAL * 1000.0f);

    printf("Transition        : %d\n",
           TEST_TRANSITION);

    printf("Transition sample : %d\n",
           transitions[TEST_TRANSITION].sample);

    printf("Before tone       : %d\n",
           transitions[TEST_TRANSITION].before_tone);

    printf("After tone        : %d\n",
           transitions[TEST_TRANSITION].after_tone);

    printf("============================================================\n");


    // ========================================================
    // Generate synthetic RF.
    // ========================================================

    std::vector<float> rf;

    generate_ft8_continuous(
        rf,
        TEST_DELAY_LOCAL,
        TEST_FREQ,
        AMPLITUDE);


    // ========================================================
    // PATH A: RF
    // ========================================================

    printf("\n");
    printf("============================================================\n");
    printf("PATH A - RF\n");
    printf("============================================================\n");

    const int transition_sample =
        transitions[TEST_TRANSITION].sample +
        (int)lroundf(
            TEST_DELAY_LOCAL *
            SAMPLE_RATE);

    prepare_local_baseband(
        rf,
        transition_sample,
        transitions[TEST_TRANSITION].before_tone,
        transitions[TEST_TRANSITION].after_tone);

    //
    // RF +sin
    //

    const V3GResult rf_plus =
        run_one_correlation_test(
            TEST_TRANSITION,
            +1);

    //
    // RF -sin
    //

    const V3GResult rf_minus =
        run_one_correlation_test(
            TEST_TRANSITION,
            -1);

    printf("\n");
    printf("RF + reference Q = +sin\n");
    printf("  best offset : %+d samples\n",
           rf_plus.best_offset);

    printf("  best score  : %.9f\n",
           rf_plus.best_score);

    printf("  time offset : %+.6f ms\n",
           1000.0f *
           rf_plus.best_offset /
           SAMPLE_RATE);

    printf("\n");

    printf("RF + reference Q = -sin\n");
    printf("  best offset : %+d samples\n",
           rf_minus.best_offset);

    printf("  best score  : %.9f\n",
           rf_minus.best_score);

    printf("  time offset : %+.6f ms\n",
           1000.0f *
           rf_minus.best_offset /
           SAMPLE_RATE);


    // ========================================================
    // PATH B: IDEAL BASEBAND
    // ========================================================

    printf("\n");
    printf("============================================================\n");
    printf("PATH B - IDEAL COMPLEX BASEBAND\n");
    printf("============================================================\n");

    prepare_ideal_baseband(
        transitions[TEST_TRANSITION].before_tone,
        transitions[TEST_TRANSITION].after_tone);

    //
    // Ideal BB +sin
    //

    const V3GResult bb_plus =
        run_one_correlation_test(
            TEST_TRANSITION,
            +1);

    //
    // Ideal BB -sin
    //

    const V3GResult bb_minus =
        run_one_correlation_test(
            TEST_TRANSITION,
            -1);

    printf("\n");
    printf("Ideal BB + reference Q = +sin\n");
    printf("  best offset : %+d samples\n",
           bb_plus.best_offset);

    printf("  best score  : %.9f\n",
           bb_plus.best_score);

    printf("  time offset : %+.6f ms\n",
           1000.0f *
           bb_plus.best_offset /
           SAMPLE_RATE);

    printf("\n");

    printf("Ideal BB + reference Q = -sin\n");
    printf("  best offset : %+d samples\n",
           bb_minus.best_offset);

    printf("  best score  : %.9f\n",
           bb_minus.best_score);

    printf("  time offset : %+.6f ms\n",
           1000.0f *
           bb_minus.best_offset /
           SAMPLE_RATE);


    // ========================================================
    // SUMMARY TABLE
    // ========================================================

    printf("\n");
    printf("============================================================\n");
    printf("V3G SUMMARY\n");
    printf("============================================================\n");

    printf(
        "%-18s %16s %16s %16s\n",
        "Path / Ref",
        "Best offset",
        "Time ms",
        "Score");

    printf(
        "--------------------------------------------------------------------------\n");

    printf(
        "%-18s %+16d %16.6f %16.9f\n",
        "RF / +sin",
        rf_plus.best_offset,
        1000.0f *
            rf_plus.best_offset /
            SAMPLE_RATE,
        rf_plus.best_score);

    printf(
        "%-18s %+16d %16.6f %16.9f\n",
        "RF / -sin",
        rf_minus.best_offset,
        1000.0f *
            rf_minus.best_offset /
            SAMPLE_RATE,
        rf_minus.best_score);

    printf(
        "%-18s %+16d %16.6f %16.9f\n",
        "IdealBB / +sin",
        bb_plus.best_offset,
        1000.0f *
            bb_plus.best_offset /
            SAMPLE_RATE,
        bb_plus.best_score);

    printf(
        "%-18s %+16d %16.6f %16.9f\n",
        "IdealBB / -sin",
        bb_minus.best_offset,
        1000.0f *
            bb_minus.best_offset /
            SAMPLE_RATE,
        bb_minus.best_score);

    printf(
        "==========================================================================\n");


    // ========================================================
    // DIFFERENTIAL ANALYSIS
    // ========================================================

    printf("\n");
    printf("============================================================\n");
    printf("V3G DIFFERENTIAL ANALYSIS\n");
    printf("============================================================\n");

    printf(
        "RF Q-sign difference       : %+d samples\n",
        rf_minus.best_offset -
        rf_plus.best_offset);

    printf(
        "IdealBB Q-sign difference  : %+d samples\n",
        bb_minus.best_offset -
        bb_plus.best_offset);

    printf(
        "RF-vs-BB difference (+sin) : %+d samples\n",
        rf_plus.best_offset -
        bb_plus.best_offset);

    printf(
        "RF-vs-BB difference (-sin) : %+d samples\n",
        rf_minus.best_offset -
        bb_minus.best_offset);


    // ========================================================
    // INTERPRETATION
    // ========================================================

    printf("\n");
    printf("============================================================\n");
    printf("V3G INTERPRETATION\n");
    printf("============================================================\n");

    const bool rf_plus_shifted =
        (rf_plus.best_offset != 0);

    const bool rf_minus_zero =
        (rf_minus.best_offset == 0);

    const bool bb_plus_shifted =
        (bb_plus.best_offset != 0);

    const bool bb_minus_zero =
        (bb_minus.best_offset == 0);


    if (rf_plus_shifted &&
        rf_minus_zero &&
        bb_plus_shifted &&
        bb_minus_zero)
    {
        printf(
            "STRONG RESULT:\n");

        printf(
            "  +sin produces the systematic offset.\n");

        printf(
            "  -sin moves both RF and IdealBB to zero.\n");

        printf(
            "  This strongly indicates a Q-sign /\n");

        printf(
            "  phase-convention mismatch in the reference.\n");
    }
    else if (rf_plus_shifted &&
             rf_minus_zero &&
             bb_plus_shifted &&
             !bb_minus_zero)
    {
        printf(
            "RESULT:\n");

        printf(
            "  RF follows the expected Q-sign behavior,\n");

        printf(
            "  but IdealBB does not.\n");

        printf(
            "  Investigate the RF/BB construction separately.\n");
    }
    else if (rf_plus_shifted &&
             !rf_minus_zero &&
             bb_plus_shifted &&
             !bb_minus_zero)
    {
        printf(
            "RESULT:\n");

        printf(
            "  Bias persists for both Q signs.\n");

        printf(
            "  Therefore Q sign alone does not explain it.\n");

        printf(
            "  Investigate FIR/reference/correlation geometry.\n");
    }
    else if (!rf_plus_shifted &&
             !bb_plus_shifted)
    {
        printf(
            "RESULT:\n");

        printf(
            "  No integer-sample systematic bias detected.\n");
    }
    else
    {
        printf(
            "RESULT:\n");

        printf(
            "  Mixed / unexpected behavior.\n");

        printf(
            "  The full correlation curves should be inspected.\n");
    }


    // ========================================================
    // FULL CORRELATION CURVES
    // ========================================================
    //
    // These are intentionally printed LAST so the important
    // four-way result appears first in the serial log.
    //

    printf("\n");
    printf("============================================================\n");
    printf("V3G FULL CORRELATION CURVES\n");
    printf("============================================================\n");

    //
    // RF +sin
    //

    prepare_local_baseband(
        rf,
        transition_sample,
        transitions[TEST_TRANSITION].before_tone,
        transitions[TEST_TRANSITION].after_tone);

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("RF / +sin\n");
    printf("------------------------------------------------------------\n");

    print_correlation_curve(
        TEST_TRANSITION,
        +1);

    //
    // RF -sin
    //

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("RF / -sin\n");
    printf("------------------------------------------------------------\n");

    print_correlation_curve(
        TEST_TRANSITION,
        -1);

    //
    // Ideal BB +sin
    //

    prepare_ideal_baseband(
        transitions[TEST_TRANSITION].before_tone,
        transitions[TEST_TRANSITION].after_tone);

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("IdealBB / +sin\n");
    printf("------------------------------------------------------------\n");

    print_correlation_curve(
        TEST_TRANSITION,
        +1);

    //
    // Ideal BB -sin
    //

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("IdealBB / -sin\n");
    printf("------------------------------------------------------------\n");

    print_correlation_curve(
        TEST_TRANSITION,
        -1);

    printf("\n");
    printf("============================================================\n");
    printf("END V3G ISOLATION TEST\n");
    printf("============================================================\n");
}


// ============================================================
// BENCHMARK
// ============================================================

static void run_benchmark()
{
    printf("\n");
    printf("============================================================\n");
    printf("V3G FT8 DELAY ESTIMATOR BENCHMARK\n");
    printf("============================================================\n");

    printf("Sample rate       : %d Hz\n",
           SAMPLE_RATE);

    printf("Symbol duration   : %d samples = %.3f ms\n",
           SPS,
           1000.0f * SPS / SAMPLE_RATE);

    printf("Test frequency    : %.3f Hz\n",
           TEST_FREQ);

    printf("Test delay        : %.3f ms\n",
           TEST_DELAY * 1000.0f);

    printf("Amplitude         : %.1f\n",
           AMPLITUDE);

    printf("Iterations        : %d\n",
           NITER);

    printf("============================================================\n");


    // ========================================================
    // Initialization
    // ========================================================

    init_fir();

    const int ntrans =
        build_transitions(
            transitions,
            NTRANS);

    if (ntrans != NTRANS)
    {
        printf(
            "ERROR: expected %d transitions, got %d\n",
            NTRANS,
            ntrans);

        return;
    }

    ensure_reference_cache();


    // ========================================================
    // Generate test signal BEFORE benchmark.
    // ========================================================

    std::vector<float> signal;

    generate_ft8_continuous(
        signal,
        TEST_DELAY,
        TEST_FREQ,
        AMPLITUDE);


    // ========================================================
    // Diagnostic
    // ========================================================

    print_v3e_diagnostic();


    // ========================================================
    // Delay sweep
    // ========================================================

    run_delay_sweep();


    // ========================================================
    // V3G isolation test
    // ========================================================

    run_v3g_isolation_test();


    // ========================================================
    // Warm-up
    // ========================================================

    printf("\n");
    printf("============================================================\n");
    printf("WARMUP\n");
    printf("============================================================\n");

    volatile float dummy = 0.0f;

    for (int i = 0; i < NWARMUP; ++i)
    {
        dummy +=
            refine_ft8_delay(
                signal,
                TEST_DELAY);
    }


    // ========================================================
    // Benchmark
    // ========================================================

    printf("\n");
    printf("============================================================\n");
    printf("BENCHMARK\n");
    printf("============================================================\n");

#ifdef ARDUINO

    const uint32_t t0 =
        micros();

    for (int i = 0; i < NITER; ++i)
    {
        dummy +=
            refine_ft8_delay(
                signal,
                TEST_DELAY);
    }

    const uint32_t t1 =
        micros();

    const double total_ms =
        (double)(t1 - t0) /
        1000.0;

#else

    struct timespec ts0;
    struct timespec ts1;

    clock_gettime(
        CLOCK_MONOTONIC,
        &ts0);

    for (int i = 0; i < NITER; ++i)
    {
        dummy +=
            refine_ft8_delay(
                signal,
                TEST_DELAY);
    }

    clock_gettime(
        CLOCK_MONOTONIC,
        &ts1);

    const double total_ms =
        (double)(ts1.tv_sec - ts0.tv_sec) *
            1000.0 +
        (double)(ts1.tv_nsec - ts0.tv_nsec) /
            1000000.0;

#endif

    const double average_ms =
        total_ms /
        (double)NITER;

    const double throughput =
        1000.0 /
        average_ms;

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Results\n");
    printf("------------------------------------------------------------\n");

    printf("Total time         : %12.3f ms\n",
           total_ms);

    printf("Average time       : %12.3f ms/call\n",
           average_ms);

    printf("Throughput         : %12.3f calls/sec\n",
           throughput);

    printf("Dummy              : %12.6f\n",
           (double)dummy);

    printf("------------------------------------------------------------\n");


    // ========================================================
    // Final estimate
    // ========================================================

    const float final_estimate =
        refine_ft8_delay(
            signal,
            TEST_DELAY);

    printf("\n");
    printf("Final estimate\n");
    printf("------------------------------------------------------------\n");

    printf("True delay         : %.6f ms\n",
           TEST_DELAY * 1000.0f);

    printf("Best offset        : %+d samples\n",
           g_debug.best_offset);

    printf("Refined offset     : %.6f samples\n",
           g_debug.refined_offset);

    printf("Estimated delay    : %.6f ms\n",
           final_estimate * 1000.0f);

    printf("Error              : %.6f samples\n",
           (final_estimate - TEST_DELAY) *
           SAMPLE_RATE);

    printf("Error              : %.6f ms\n",
           (final_estimate - TEST_DELAY) *
           1000.0f);

    printf("============================================================\n");
}


// ============================================================
// ARDUINO
// ============================================================

#ifdef ARDUINO

void setup()
{
    delay(1000);

    printf("\n");
    printf("\n");

    run_benchmark();
}

void loop()
{
    delay(1000);
}

#else

// ============================================================
// LINUX
// ============================================================

int main()
{
    run_benchmark();

    return 0;
}

#endif