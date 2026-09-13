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
#include <fstream>
#include <cstdint>
#include <cstring>
#endif

#include <math.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>

#ifndef PI
#define PI 3.14159265358979323846f
#endif
#ifndef TWO_PI
#define TWO_PI 6.28318530717958647692f
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

double elapsed_ms(const struct timespec *t0, const struct timespec *t1)
{
  long sec = t1->tv_sec - t0->tv_sec;
  long nsec = t1->tv_nsec - t0->tv_nsec;
  return (double)sec * 1000.0 + (double)nsec / 1000000.0;
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
// SYMBOL SELECTION
// ============================================================

enum class SymbolMode
{
    FULL,
    COSTAS,
    COSTAS_DATA
};


// ------------------------------------------------------------
// FT8 Costas positions
//
//   0 ... 6
//   36 ... 42
//   72 ... 78
//
// Total = 21 symbols
// ------------------------------------------------------------

static inline bool is_costas_symbol(int k)
{
    return
        (k >= 0  && k <= 6)  ||
        (k >= 36 && k <= 42) ||
        (k >= 72 && k <= 78);
}


// ------------------------------------------------------------
// Selected additional data symbols
//
// Total additional = 12
//
// Together with Costas:
//
//     21 + 12 = 33 symbols
// ------------------------------------------------------------

static inline bool is_selected_data_symbol(int k)
{
    switch (k)
    {
        case 8:
        case 12:
        case 16:
        case 20:
        case 24:
        case 28:
        case 32:
        case 45:
        case 50:
        case 55:
        case 60:
        case 65:
            return true;

        default:
            return false;
    }
}


// ------------------------------------------------------------
// Decide whether symbol k participates in the correlation.
// ------------------------------------------------------------

static inline bool use_symbol(
    int k,
    SymbolMode mode)
{
    switch (mode)
    {
        case SymbolMode::FULL:
            return true;

        case SymbolMode::COSTAS:
            return is_costas_symbol(k);

        case SymbolMode::COSTAS_DATA:
            return
                is_costas_symbol(k) ||
                is_selected_data_symbol(k);
    }

    return false;
}


// ------------------------------------------------------------
// Number of symbols actually used.
// ------------------------------------------------------------

static int count_used_symbols(SymbolMode mode)
{
    int count = 0;

    for (int k = 0; k < NTONES; ++k)
    {
        if (use_symbol(k, mode))
            ++count;
    }

    return count;
}


// ------------------------------------------------------------
// Mode name
// ------------------------------------------------------------

static const char* symbol_mode_name(
    SymbolMode mode)
{
    switch (mode)
    {
        case SymbolMode::FULL:
            return "FULL";

        case SymbolMode::COSTAS:
            return "COSTAS";

        case SymbolMode::COSTAS_DATA:
            return "COSTAS+DATA";
    }

    return "?";
}


// ============================================================
// V10G - EXACT CORRELATION DIAGNOSTIC
//
// Confronta, sulla STESSA identica finestra:
//
//   1) DIRECT   : cosf/sinf per ogni campione
//   2) GROUP4   : initial_correlation() della V8.1
//   3) SLIDING  : update incrementale della V8.1
//
// Scopo:
//   capire ESATTAMENTE da dove nasce il +25 samples.
//
// Costas: primo simbolo (k=0)
// Delay: center-40 ... center+40
//
// ============================================================

float refine_ft8_delay_v10g(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)cand_to_print;

    constexpr int SAMPLE_RATE    = 12000;
    constexpr int SYMBOL_SAMPLES = 1920;

    constexpr int RADIUS = 40;

    // --------------------------------------------------------
    // Costas symbol 0
    // --------------------------------------------------------

    constexpr int k = 0;

    const float f =
        freq +
        TONE_SPACING * (float)tones[k];

    const float w =
        TWO_PI * f /
        (float)SAMPLE_RATE;

    const float pc = cosf(w);
    const float ps = sinf(w);

    const float end_angle =
        -w * (float)SYMBOL_SAMPLES;

    const float ec = cosf(end_angle);
    const float es = sinf(end_angle);

    const int center =
        (int)lroundf(
            delay0 * (float)SAMPLE_RATE);

    const int first =
        center - RADIUS;

    const int last =
        center + RADIUS;

    // ========================================================
    // DIRECT correlation
    //
    // Questa è la reference assoluta.
    //
    // C = sum x[n] exp(-j*w*n)
    // ========================================================

    auto direct_correlation =
        [&](int start,
            float& out_i,
            float& out_q)
    {
        out_i = 0.0f;
        out_q = 0.0f;

        for (int n = 0;
             n < SYMBOL_SAMPLES;
             ++n)
        {
            const float x =
                samples[start + n];

            const float a =
                w * (float)n;

            out_i +=
                x * cosf(a);

            out_q -=
                x * sinf(a);
        }
    };

    // ========================================================
    // EXACT initial_correlation() FROM V8.1
    // ========================================================

    auto group4_correlation =
        [&](int start,
            float& out_i,
            float& out_q)
    {
        out_i = 0.0f;
        out_q = 0.0f;

        const float cw = cosf(w);
        const float sw = sinf(w);

        const float c0 = 1.0f;
        const float s0 = 0.0f;

        const float c1 =
            c0 * cw + s0 * sw;

        const float s1 =
            s0 * cw - c0 * sw;

        const float c2 =
            c1 * cw + s1 * sw;

        const float s2 =
            s1 * cw - c1 * sw;

        const float c3 =
            c2 * cw + s2 * sw;

        const float s3 =
            s2 * cw - c2 * sw;

        const float p4c =
            c3 * cw + s3 * sw;

        const float p4s =
            s3 * cw - c3 * sw;

        float gc = 1.0f;
        float gs = 0.0f;

        int pos = start;
        int remaining = SYMBOL_SAMPLES;

        while (remaining >= 4)
        {
            const float x0 =
                samples[pos + 0];

            const float x1 =
                samples[pos + 1];

            const float x2 =
                samples[pos + 2];

            const float x3 =
                samples[pos + 3];

            const float li =
                x0 +
                x1 * c1 +
                x2 * c2 +
                x3 * c3;

            const float lq =
                x1 * s1 +
                x2 * s2 +
                x3 * s3;

            out_i +=
                li * gc -
                lq * gs;

            out_q +=
                li * gs +
                lq * gc;

            const float ngc =
                gc * p4c -
                gs * p4s;

            const float ngs =
                gc * p4s +
                gs * p4c;

            gc = ngc;
            gs = ngs;

            pos += 4;
            remaining -= 4;
        }

        if (remaining > 0)
        {
            float li = 0.0f;
            float lq = 0.0f;

            if (remaining >= 1)
            {
                const float x0 =
                    samples[pos];

                li += x0;
            }

            if (remaining >= 2)
            {
                const float x1 =
                    samples[pos + 1];

                li += x1 * c1;
                lq += x1 * s1;
            }

            if (remaining >= 3)
            {
                const float x2 =
                    samples[pos + 2];

                li += x2 * c2;
                lq += x2 * s2;
            }

            out_i +=
                li * gc -
                lq * gs;

            out_q +=
                li * gs +
                lq * gc;
        }
    };

    // ========================================================
    // HEADER
    // ========================================================

    Serial.println();
    Serial.println("============================================================");
    Serial.println("V10G - EXACT CORRELATION DIAGNOSTIC");
    Serial.println("============================================================");

    Serial.printf(
        "Center delay      : %d samples\n",
        center);

    Serial.printf(
        "Frequency         : %.6f Hz\n",
        f);

    Serial.printf(
        "Symbol            : %d\n",
        k);

    Serial.printf(
        "Scan              : %d ... %d\n",
        first,
        last);

    Serial.println();
    Serial.println(
        "delay     DIRECT_SCORE     GROUP4_SCORE     SLIDE_SCORE"
        "       D-G        S-G");

    Serial.println(
        "------------------------------------------------------------");

    // ========================================================
    // INITIALIZE SLIDING AT FIRST
    //
    // THIS IS EXACTLY THE V8.1 STARTING POINT.
    // ========================================================

    float slide_i = 0.0f;
    float slide_q = 0.0f;

    group4_correlation(
        first,
        slide_i,
        slide_q);

    int pos = first;

    // ========================================================
    // TRACK MAXIMA
    // ========================================================

    float best_direct = -1.0f;
    float best_group4 = -1.0f;
    float best_slide  = -1.0f;

    int best_direct_delay = first;
    int best_group4_delay = first;
    int best_slide_delay  = first;

    // ========================================================
    // SCAN
    // ========================================================

    for (int d = first;
         d <= last;
         ++d)
    {
        // ----------------------------------------------------
        // DIRECT
        // ----------------------------------------------------

        float di;
        float dq;

        direct_correlation(
            d,
            di,
            dq);

        const float direct_score =
            di * di +
            dq * dq;

        // ----------------------------------------------------
        // GROUP4
        // ----------------------------------------------------

        float gi;
        float gq;

        group4_correlation(
            d,
            gi,
            gq);

        const float group4_score =
            gi * gi +
            gq * gq;

        // ----------------------------------------------------
        // SLIDING
        //
        // At d == first, slide_i/q already correspond
        // exactly to this window.
        // ----------------------------------------------------

        const float slide_score =
            slide_i * slide_i +
            slide_q * slide_q;

        // ----------------------------------------------------
        // Differences
        // ----------------------------------------------------

        const float dg =
            direct_score -
            group4_score;

        const float sg =
            slide_score -
            group4_score;

        // ----------------------------------------------------
        // Print
        //
        // Print every point, but highlight center and +25.
        // ----------------------------------------------------

        if (d == center ||
            d == center + 25 ||
            d == center - 25)
        {
            Serial.printf(
                "%5d *  %14.3f  %14.3f  %14.3f"
                "  %+10.3f  %+10.3f\n",
                d,
                direct_score,
                group4_score,
                slide_score,
                dg,
                sg);
        }

        // ----------------------------------------------------
        // Max DIRECT
        // ----------------------------------------------------

        if (direct_score > best_direct)
        {
            best_direct =
                direct_score;

            best_direct_delay =
                d;
        }

        // ----------------------------------------------------
        // Max GROUP4
        // ----------------------------------------------------

        if (group4_score > best_group4)
        {
            best_group4 =
                group4_score;

            best_group4_delay =
                d;
        }

        // ----------------------------------------------------
        // Max SLIDING
        // ----------------------------------------------------

        if (slide_score > best_slide)
        {
            best_slide =
                slide_score;

            best_slide_delay =
                d;
        }

        // ----------------------------------------------------
        // Advance sliding correlation
        //
        // C(d+1) =
        //
        // exp(+jw) *
        // [ C(d)
        //   - x[d]
        //   + x[d+N] exp(-jwN) ]
        //
        // EXACT V8.1 UPDATE.
        // ----------------------------------------------------

        if (d < last)
        {
            const int old_pos =
                pos;

            const int new_sample =
                old_pos +
                SYMBOL_SAMPLES;

            if (new_sample >= num_samples)
            {
                Serial.println(
                    "ERROR: sliding new_sample out of range");

                break;
            }

            const float x_old =
                samples[old_pos];

            const float x_new =
                samples[new_sample];

            const float add_i =
                x_new * ec;

            const float add_q =
                x_new * es;

            const float ti =
                slide_i -
                x_old +
                add_i;

            const float tq =
                slide_q +
                add_q;

            const float new_i =
                ti * pc -
                tq * ps;

            const float new_q =
                ti * ps +
                tq * pc;

            slide_i = new_i;
            slide_q = new_q;

            pos++;
        }
    }

    // ========================================================
    // RESULTS
    // ========================================================

    Serial.println();
    Serial.println("============================================================");
    Serial.println("MAXIMUMS");
    Serial.println("============================================================");

    Serial.printf(
        "DIRECT  : delay=%d  offset=%+d  score=%.6f\n",
        best_direct_delay,
        best_direct_delay - center,
        best_direct);

    Serial.printf(
        "GROUP4  : delay=%d  offset=%+d  score=%.6f\n",
        best_group4_delay,
        best_group4_delay - center,
        best_group4);

    Serial.printf(
        "SLIDING : delay=%d  offset=%+d  score=%.6f\n",
        best_slide_delay,
        best_slide_delay - center,
        best_slide);

    Serial.println();

    // ========================================================
    // EXTRA CHECK:
    //
    // Compare DIRECT vs GROUP4 exactly at center.
    // ========================================================

    float di0, dq0;
    float gi0, gq0;

    direct_correlation(
        center,
        di0,
        dq0);

    group4_correlation(
        center,
        gi0,
        gq0);

    Serial.println(
        "============================================================");
    Serial.println(
        "CENTER CORRELATION");
    Serial.println(
        "============================================================");

    Serial.printf(
        "DIRECT : I=%+.6f Q=%+.6f SCORE=%.6f\n",
        di0,
        dq0,
        di0 * di0 + dq0 * dq0);

    Serial.printf(
        "GROUP4 : I=%+.6f Q=%+.6f SCORE=%.6f\n",
        gi0,
        gq0,
        gi0 * gi0 + gq0 * gq0);

    Serial.printf(
        "DIFF   : dI=%+.9f dQ=%+.9f\n",
        gi0 - di0,
        gq0 - dq0);

    Serial.println();

    return
        best_slide_delay /
        (float)SAMPLE_RATE;
}

float refine_ft8_delay_v10f(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)tones;
    (void)freq;
    (void)cand_to_print;

    constexpr int SAMPLE_RATE    = 12000;
    constexpr int SYMBOL_SAMPLES = 1920;

    const int center =
        (int)lroundf(delay0 * (float)SAMPLE_RATE);

    Serial.println();
    Serial.println("============================================================");
    Serial.println("V10F - RAW SAMPLE / WINDOW DIAGNOSTIC");
    Serial.println("============================================================");
    Serial.printf("True delay : %d samples\n", center);
    Serial.printf("Samples    : %d\n", num_samples);
    Serial.printf("Symbol     : %d samples\n", SYMBOL_SAMPLES);
    Serial.println();

    // --------------------------------------------------------
    // Analisi di una finestra
    // --------------------------------------------------------

    auto analyze_window =
        [&](const char* name, int start)
    {
        Serial.println();
        Serial.println("------------------------------------------------------------");
        Serial.printf(
            "%s  start=%d  end=%d\n",
            name,
            start,
            start + SYMBOL_SAMPLES - 1);
        Serial.println("------------------------------------------------------------");

        if (start < 0 ||
            start + SYMBOL_SAMPLES > num_samples)
        {
            Serial.println("INVALID WINDOW");
            return;
        }

        float min_val = samples[start];
        float max_val = samples[start];

        double sum = 0.0;
        double sum_abs = 0.0;
        double sum_sq = 0.0;

        // ----------------------------------------------------
        // Primi 16 campioni
        // ----------------------------------------------------

        Serial.println("First 16 samples:");

        for (int n = 0; n < 16; ++n)
        {
            Serial.printf(
                "  [%4d] = %+12.8f\n",
                start + n,
                samples[start + n]);
        }

        // ----------------------------------------------------
        // Statistiche
        // ----------------------------------------------------

        for (int n = 0; n < SYMBOL_SAMPLES; ++n)
        {
            const float x =
                samples[start + n];

            if (x < min_val)
                min_val = x;

            if (x > max_val)
                max_val = x;

            sum += (double)x;
            sum_abs += fabs((double)x);
            sum_sq += (double)x * (double)x;
        }

        const double mean =
            sum / (double)SYMBOL_SAMPLES;

        const double rms =
            sqrt(sum_sq / (double)SYMBOL_SAMPLES);

        Serial.println();
        Serial.printf(
            "MIN       : %+14.8e\n",
            min_val);

        Serial.printf(
            "MAX       : %+14.8e\n",
            max_val);

        Serial.printf(
            "MEAN      : %+14.8e\n",
            mean);

        Serial.printf(
            "RMS       : %14.8e\n",
            rms);

        Serial.printf(
            "SUM       : %+14.8e\n",
            sum);

        Serial.printf(
            "SUM ABS   : %14.8e\n",
            sum_abs);

        // ----------------------------------------------------
        // Hash semplice della finestra
        //
        // Usiamo i bit del float, non il valore convertito.
        // Serve solo per verificare se due finestre sono
        // identiche.
        // ----------------------------------------------------

        uint32_t hash = 2166136261u;

        for (int n = 0; n < SYMBOL_SAMPLES; ++n)
        {
            uint32_t bits;

            memcpy(
                &bits,
                &samples[start + n],
                sizeof(bits));

            hash ^= bits;
            hash *= 16777619u;
        }

        Serial.printf(
            "HASH      : 0x%08lX\n",
            (unsigned long)hash);

        // ----------------------------------------------------
        // Alcuni campioni distribuiti nella finestra
        // ----------------------------------------------------

        Serial.println();
        Serial.println("Selected samples:");

        const int positions[] =
        {
            0,
            1,
            2,
            3,
            10,
            100,
            500,
            960,
            1000,
            1500,
            1919
        };

        constexpr int NPOS =
            sizeof(positions) / sizeof(positions[0]);

        for (int i = 0; i < NPOS; ++i)
        {
            const int n = positions[i];

            Serial.printf(
                "  n=%4d  sample[%d] = %+12.8f\n",
                n,
                start + n,
                samples[start + n]);
        }
    };

    // --------------------------------------------------------
    // Le tre finestre corrispondenti ai tre test:
    //
    // D1500 -> 18000
    // D2000 -> 24000
    // D2500 -> 30000
    //
    // La funzione usa comunque delay0 come centro della
    // finestra principale.
    // --------------------------------------------------------

    const int start0 =
        center;

    const int start_minus_6000 =
        center - 6000;

    const int start_plus_6000 =
        center + 6000;

    Serial.println("REFERENCE WINDOWS");
    Serial.printf(
        "center       = %d\n",
        center);

    Serial.printf(
        "center-6000  = %d\n",
        start_minus_6000);

    Serial.printf(
        "center       = %d\n",
        start0);

    Serial.printf(
        "center+6000  = %d\n",
        start_plus_6000);

    // --------------------------------------------------------
    // Analizza le tre finestre
    // --------------------------------------------------------

    analyze_window(
        "WINDOW -6000",
        start_minus_6000);

    analyze_window(
        "WINDOW CENTER",
        start0);

    analyze_window(
        "WINDOW +6000",
        start_plus_6000);

    // --------------------------------------------------------
    // Confronto diretto tra finestre
    // --------------------------------------------------------

    Serial.println();
    Serial.println("============================================================");
    Serial.println("DIRECT WINDOW COMPARISON");
    Serial.println("============================================================");

    if (start_minus_6000 >= 0 &&
        start_minus_6000 + SYMBOL_SAMPLES <= num_samples &&
        start0 >= 0 &&
        start0 + SYMBOL_SAMPLES <= num_samples &&
        start_plus_6000 >= 0 &&
        start_plus_6000 + SYMBOL_SAMPLES <= num_samples)
    {
        int equal_a_b = 1;
        int equal_b_c = 1;
        int equal_a_c = 1;

        double diff_a_b = 0.0;
        double diff_b_c = 0.0;
        double diff_a_c = 0.0;

        double max_a_b = 0.0;
        double max_b_c = 0.0;
        double max_a_c = 0.0;

        int first_diff_a_b = -1;
        int first_diff_b_c = -1;
        int first_diff_a_c = -1;

        for (int n = 0;
             n < SYMBOL_SAMPLES;
             ++n)
        {
            const float a =
                samples[start_minus_6000 + n];

            const float b =
                samples[start0 + n];

            const float c =
                samples[start_plus_6000 + n];

            const double dab =
                fabs((double)a - (double)b);

            const double dbc =
                fabs((double)b - (double)c);

            const double dac =
                fabs((double)a - (double)c);

            diff_a_b += dab;
            diff_b_c += dbc;
            diff_a_c += dac;

            if (dab > max_a_b)
                max_a_b = dab;

            if (dbc > max_b_c)
                max_b_c = dbc;

            if (dac > max_a_c)
                max_a_c = dac;

            if (a != b)
            {
                equal_a_b = 0;

                if (first_diff_a_b < 0)
                    first_diff_a_b = n;
            }

            if (b != c)
            {
                equal_b_c = 0;

                if (first_diff_b_c < 0)
                    first_diff_b_c = n;
            }

            if (a != c)
            {
                equal_a_c = 0;

                if (first_diff_a_c < 0)
                    first_diff_a_c = n;
            }
        }

        Serial.println();

        Serial.printf(
            "WINDOW A = center-6000 : %d\n",
            start_minus_6000);

        Serial.printf(
            "WINDOW B = center      : %d\n",
            start0);

        Serial.printf(
            "WINDOW C = center+6000 : %d\n",
            start_plus_6000);

        Serial.println();

        Serial.printf(
            "A == B          : %s\n",
            equal_a_b ? "YES" : "NO");

        Serial.printf(
            "B == C          : %s\n",
            equal_b_c ? "YES" : "NO");

        Serial.printf(
            "A == C          : %s\n",
            equal_a_c ? "YES" : "NO");

        Serial.println();

        Serial.printf(
            "A-B sum abs diff: %14.8e\n",
            diff_a_b);

        Serial.printf(
            "B-C sum abs diff: %14.8e\n",
            diff_b_c);

        Serial.printf(
            "A-C sum abs diff: %14.8e\n",
            diff_a_c);

        Serial.println();

        Serial.printf(
            "A-B max diff    : %14.8e\n",
            max_a_b);

        Serial.printf(
            "B-C max diff    : %14.8e\n",
            max_b_c);

        Serial.printf(
            "A-C max diff    : %14.8e\n",
            max_a_c);

        Serial.println();

        Serial.printf(
            "First A-B diff  : %d\n",
            first_diff_a_b);

        Serial.printf(
            "First B-C diff  : %d\n",
            first_diff_b_c);

        Serial.printf(
            "First A-C diff  : %d\n",
            first_diff_a_c);
    }
    else
    {
        Serial.println(
            "One or more windows are outside the buffer.");
    }

    // --------------------------------------------------------
    // IMPORTANTISSIMO:
    //
    // Confrontiamo anche la stessa finestra traslata di
    // pochi campioni. Questo ci dice immediatamente se
    // samples[] contiene realmente la forma d'onda.
    // --------------------------------------------------------

    Serial.println();
    Serial.println("============================================================");
    Serial.println("LOCAL SAMPLE SHIFT CHECK");
    Serial.println("============================================================");

    const int test_start =
        center;

    if (test_start >= 0 &&
        test_start + SYMBOL_SAMPLES + 4 <= num_samples)
    {
        for (int shift = -4;
             shift <= 4;
             ++shift)
        {
            const int s =
                test_start + shift;

            double diff = 0.0;

            for (int n = 0;
                 n < SYMBOL_SAMPLES;
                 ++n)
            {
                diff += fabs(
                    (double)samples[test_start + n] -
                    (double)samples[s + n]);
            }

            Serial.printf(
                "shift=%+d  start=%d  sum_abs_diff=%14.8e\n",
                shift,
                s,
                diff);
        }
    }

    Serial.println();
    Serial.println("============================================================");
    Serial.println("V10F COMPLETE");
    Serial.println("============================================================");

    // V10F è diagnostica: restituisce semplicemente delay0.
    return delay0;
}

float refine_ft8_delay_v10e(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)cand_to_print;

    constexpr int SAMPLE_RATE    = 12000;
    constexpr int SYMBOL_SAMPLES = 1920;
    constexpr int NCOSTAS        = 7;
    constexpr int RADIUS         = 40;
    constexpr float TONE_SPACING = 6.25f;

    // ------------------------------------------------------------
    // IDENTICA initial_correlation() usata in V8.1 / V10C / V10D
    // ------------------------------------------------------------
    auto initial_correlation =
        [&](int start,
            int count,
            float w,
            float& out_i,
            float& out_q)
    {
        out_i = 0.0f;
        out_q = 0.0f;

        if (count <= 0)
            return;

        const float cw = cosf(w);
        const float sw = sinf(w);

        const float c0 = 1.0f;
        const float s0 = 0.0f;

        const float c1 =
            c0 * cw + s0 * sw;
        const float s1 =
            s0 * cw - c0 * sw;

        const float c2 =
            c1 * cw + s1 * sw;
        const float s2 =
            s1 * cw - c1 * sw;

        const float c3 =
            c2 * cw + s2 * sw;
        const float s3 =
            s2 * cw - c2 * sw;

        const float p4c =
            c3 * cw + s3 * sw;

        const float p4s =
            s3 * cw - c3 * sw;

        float gc = 1.0f;
        float gs = 0.0f;

        int pos = start;
        int remaining = count;

        while (remaining >= 4)
        {
            const float x0 = samples[pos + 0];
            const float x1 = samples[pos + 1];
            const float x2 = samples[pos + 2];
            const float x3 = samples[pos + 3];

            const float li =
                x0 +
                x1 * c1 +
                x2 * c2 +
                x3 * c3;

            const float lq =
                x1 * s1 +
                x2 * s2 +
                x3 * s3;

            out_i += li * gc - lq * gs;
            out_q += li * gs + lq * gc;

            const float ngc =
                gc * p4c -
                gs * p4s;

            const float ngs =
                gc * p4s +
                gs * p4c;

            gc = ngc;
            gs = ngs;

            pos += 4;
            remaining -= 4;
        }

        if (remaining > 0)
        {
            float li = 0.0f;
            float lq = 0.0f;

            if (remaining >= 1)
            {
                const float x0 = samples[pos];
                li += x0;
            }

            if (remaining >= 2)
            {
                const float x1 = samples[pos + 1];
                li += x1 * c1;
                lq += x1 * s1;
            }

            if (remaining >= 3)
            {
                const float x2 = samples[pos + 2];
                li += x2 * c2;
                lq += x2 * s2;
            }

            out_i += li * gc - lq * gs;
            out_q += li * gs + lq * gc;
        }
    };

    const int center =
        (int)lroundf(delay0 * (float)SAMPLE_RATE);

    Serial.println();
    Serial.println("============================================================");
    Serial.println("V10E - INDIVIDUAL COSTAS DELAY DIAGNOSTIC");
    Serial.println("============================================================");
    Serial.printf("True delay : %d samples\n", center);
    Serial.printf("Frequency  : %.3f Hz\n", freq);
    Serial.printf("Costas     : %d\n", NCOSTAS);
    Serial.println();

    // ------------------------------------------------------------
    // Per ogni Costas:
    //   cerca il massimo indipendentemente
    // ------------------------------------------------------------

    float costas_best_score[NCOSTAS];
    int   costas_best_delay[NCOSTAS];

    for (int k = 0; k < NCOSTAS; ++k)
    {
        float best_score = -1.0f;
        int best_delay = center;

        const float f =
            freq +
            TONE_SPACING * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        for (int offset = -RADIUS;
             offset <= RADIUS;
             ++offset)
        {
            const int delay =
                center + offset;

            const int start =
                delay + k * SYMBOL_SAMPLES;

            if (start < 0 ||
                start + SYMBOL_SAMPLES > num_samples)
                continue;

            float ci;
            float cq;

            initial_correlation(
                start,
                SYMBOL_SAMPLES,
                w,
                ci,
                cq);

            const float score =
                ci * ci + cq * cq;

            if (score > best_score)
            {
                best_score = score;
                best_delay = delay;
            }
        }

        costas_best_score[k] = best_score;
        costas_best_delay[k] = best_delay;

        Serial.printf(
            "Costas %d  tone=%d  freq=%8.3f Hz  "
            "BEST=%d  offset=%+d  score=%12.6e\n",
            k,
            tones[k],
            f,
            best_delay,
            best_delay - center,
            best_score);
    }

    // ------------------------------------------------------------
    // Ora cerca il massimo della SOMMA dei 7 Costas
    // ------------------------------------------------------------

    float total_best_score = -1.0f;
    int total_best_delay = center;

    Serial.println();
    Serial.println("------------------------------------------------------------");
    Serial.println("AGGREGATED 7-COSTAS METRIC");
    Serial.println("------------------------------------------------------------");

    for (int offset = -RADIUS;
         offset <= RADIUS;
         ++offset)
    {
        const int delay =
            center + offset;

        float total = 0.0f;
        int used = 0;

        for (int k = 0; k < NCOSTAS; ++k)
        {
            const int start =
                delay + k * SYMBOL_SAMPLES;

            if (start < 0 ||
                start + SYMBOL_SAMPLES > num_samples)
                continue;

            const float f =
                freq +
                TONE_SPACING * (float)tones[k];

            const float w =
                TWO_PI * f /
                (float)SAMPLE_RATE;

            float ci;
            float cq;

            initial_correlation(
                start,
                SYMBOL_SAMPLES,
                w,
                ci,
                cq);

            total += ci * ci + cq * cq;
            ++used;
        }

        if (used == 0)
            continue;

        const float score =
            total / (float)used;

        if (score > total_best_score)
        {
            total_best_score = score;
            total_best_delay = delay;
        }

        // Stampiamo solo +/-30 per non generare troppo output
        if (offset >= -30 && offset <= 30)
        {
            Serial.printf(
                "%6d  %+4d  %12.6e\n",
                delay,
                offset,
                score);
        }
    }

    Serial.println("------------------------------------------------------------");

    Serial.printf(
        "AGGREGATED BEST : %d samples\n",
        total_best_delay);

    Serial.printf(
        "AGGREGATED OFFSET: %+d samples\n",
        total_best_delay - center);

    Serial.printf(
        "AGGREGATED SCORE : %12.6e\n",
        total_best_score);

    Serial.println();
    Serial.println("============================================================");
    Serial.println("V10E RESULT");
    Serial.println("============================================================");

    for (int k = 0; k < NCOSTAS; ++k)
    {
        Serial.printf(
            "Costas %d : %+d samples\n",
            k,
            costas_best_delay[k] - center);
    }

    Serial.printf(
        "TOTAL     : %+d samples\n",
        total_best_delay - center);

    Serial.println("============================================================");

    return (float)total_best_delay /
           (float)SAMPLE_RATE;
}

float refine_ft8_delay_v10d(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)cand_to_print;

    constexpr int SAMPLE_RATE    = 12000;
    constexpr int SYMBOL_SAMPLES = 1920;
    constexpr int NCOSTAS        = 7;
    constexpr int RADIUS         = 40;
    constexpr float TONE_SPACING = 6.25f;

    auto initial_correlation =
        [&](int start,
            int count,
            float w,
            float& out_i,
            float& out_q)
    {
        out_i = 0.0f;
        out_q = 0.0f;

        if (count <= 0)
            return;

        const float cw = cosf(w);
        const float sw = sinf(w);

        const float c0 = 1.0f;
        const float s0 = 0.0f;

        const float c1 =
            c0 * cw + s0 * sw;
        const float s1 =
            s0 * cw - c0 * sw;

        const float c2 =
            c1 * cw + s1 * sw;
        const float s2 =
            s1 * cw - c1 * sw;

        const float c3 =
            c2 * cw + s2 * sw;
        const float s3 =
            s2 * cw - c2 * sw;

        const float p4c =
            c3 * cw + s3 * sw;
        const float p4s =
            s3 * cw - c3 * sw;

        float gc = 1.0f;
        float gs = 0.0f;

        int pos = start;
        int remaining = count;

        while (remaining >= 4)
        {
            const float x0 = samples[pos + 0];
            const float x1 = samples[pos + 1];
            const float x2 = samples[pos + 2];
            const float x3 = samples[pos + 3];

            const float li =
                x0 +
                x1 * c1 +
                x2 * c2 +
                x3 * c3;

            const float lq =
                x1 * s1 +
                x2 * s2 +
                x3 * s3;

            out_i += li * gc - lq * gs;
            out_q += li * gs + lq * gc;

            const float ngc =
                gc * p4c -
                gs * p4s;

            const float ngs =
                gc * p4s +
                gs * p4c;

            gc = ngc;
            gs = ngs;

            pos += 4;
            remaining -= 4;
        }

        if (remaining > 0)
        {
            float li = 0.0f;
            float lq = 0.0f;

            if (remaining >= 1)
            {
                const float x0 = samples[pos];
                li += x0;
            }

            if (remaining >= 2)
            {
                const float x1 = samples[pos + 1];
                li += x1 * c1;
                lq += x1 * s1;
            }

            if (remaining >= 3)
            {
                const float x2 = samples[pos + 2];
                li += x2 * c2;
                lq += x2 * s2;
            }

            out_i += li * gc - lq * gs;
            out_q += li * gs + lq * gc;
        }
    };

    const int center =
        (int)lroundf(delay0 * (float)SAMPLE_RATE);

    Serial.println();
    Serial.println("============================================================");
    Serial.println("V10D DIAGNOSTIC - RAW CORRELATION AROUND TRUE DELAY");
    Serial.println("============================================================");
    Serial.printf("Center delay : %d samples\n", center);
    Serial.printf("Frequency     : %.3f Hz\n", freq);
    Serial.printf("Costas        : %d\n", NCOSTAS);
    Serial.println();
    Serial.println(" delay       offset        score                 norm");
    Serial.println("------------------------------------------------------------");

    float best_score = -1.0f;
    int best_delay = center;

    for (int delay = center - RADIUS;
         delay <= center + RADIUS;
         ++delay)
    {
        float total = 0.0f;
        int used = 0;

        for (int k = 0; k < NCOSTAS; ++k)
        {
            const int start =
                delay + k * SYMBOL_SAMPLES;

            if (start < 0 ||
                start + SYMBOL_SAMPLES > num_samples)
                continue;

            const float f =
                freq +
                TONE_SPACING * (float)tones[k];

            const float w =
                TWO_PI * f /
                (float)SAMPLE_RATE;

            float ci;
            float cq;

            initial_correlation(
                start,
                SYMBOL_SAMPLES,
                w,
                ci,
                cq);

            total += ci * ci + cq * cq;
            ++used;
        }

        if (used == 0)
            continue;

        const float score =
            total / (float)used;

        if (score > best_score)
        {
            best_score = score;
            best_delay = delay;
        }

        // Stampiamo soprattutto la zona interessante
        if (delay >= center - 30 &&
            delay <= center + 30)
        {
            Serial.printf(
                "%6d      %+6d      %14.6e\n",
                delay,
                delay - center,
                score);
        }
    }

    Serial.println("------------------------------------------------------------");

    Serial.printf(
        "BEST DELAY : %d samples\n",
        best_delay);

    Serial.printf(
        "OFFSET     : %+d samples\n",
        best_delay - center);

    Serial.printf(
        "DELAY      : %.6f ms\n",
        1000.0f *
        (float)best_delay /
        (float)SAMPLE_RATE);

    Serial.printf(
        "ERROR      : %+d samples\n",
        best_delay - center);

    Serial.println(
        "============================================================");

    return (float)best_delay /
           (float)SAMPLE_RATE;
}

float refine_ft8_delay_v10c(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)cand_to_print;

    constexpr int SAMPLE_RATE    = 12000;
    constexpr int SYMBOL_SAMPLES = 1920;
    constexpr int COARSE_RADIUS  = 600;
    constexpr int COARSE_STEP    = 4;
    constexpr int FINE_RADIUS    = 6;

    constexpr int NCOSTAS = 7;

    constexpr int MAX_COARSE =
        2 * COARSE_RADIUS / COARSE_STEP + 1;

    constexpr int MAX_FINE = 16;

    constexpr float TONE_SPACING = 6.25f;

    // ========================================================
    // INITIAL CORRELATION
    //
    // IDENTICA a V8.1
    //
    // C = sum x[n] * exp(-j*w*n)
    //
    // ========================================================

    auto initial_correlation =
        [&](int start,
            int count,
            float w,
            float& out_i,
            float& out_q)
    {
        out_i = 0.0f;
        out_q = 0.0f;

        if (count <= 0)
            return;

        const float cw = cosf(w);
        const float sw = sinf(w);

        // exp(-j*n*w), n = 0..3

        const float c0 = 1.0f;
        const float s0 = 0.0f;

        const float c1 =
            c0 * cw + s0 * sw;

        const float s1 =
            s0 * cw - c0 * sw;

        const float c2 =
            c1 * cw + s1 * sw;

        const float s2 =
            s1 * cw - c1 * sw;

        const float c3 =
            c2 * cw + s2 * sw;

        const float s3 =
            s2 * cw - c2 * sw;

        // exp(-j*4*w)

        const float p4c =
            c3 * cw + s3 * sw;

        const float p4s =
            s3 * cw - c3 * sw;

        float gc = 1.0f;
        float gs = 0.0f;

        int pos = start;
        int remaining = count;

        // ----------------------------------------------------
        // Groups of 4
        // ----------------------------------------------------

        while (remaining >= 4)
        {
            const float x0 =
                samples[pos + 0];

            const float x1 =
                samples[pos + 1];

            const float x2 =
                samples[pos + 2];

            const float x3 =
                samples[pos + 3];

            const float li =
                x0 +
                x1 * c1 +
                x2 * c2 +
                x3 * c3;

            const float lq =
                x1 * s1 +
                x2 * s2 +
                x3 * s3;

            out_i +=
                li * gc -
                lq * gs;

            out_q +=
                li * gs +
                lq * gc;

            const float ngc =
                gc * p4c -
                gs * p4s;

            const float ngs =
                gc * p4s +
                gs * p4c;

            gc = ngc;
            gs = ngs;

            pos += 4;
            remaining -= 4;
        }

        // ----------------------------------------------------
        // Remaining samples
        // ----------------------------------------------------

        if (remaining > 0)
        {
            float li = 0.0f;
            float lq = 0.0f;

            if (remaining >= 1)
            {
                const float x0 =
                    samples[pos];

                li += x0;
            }

            if (remaining >= 2)
            {
                const float x1 =
                    samples[pos + 1];

                li += x1 * c1;
                lq += x1 * s1;
            }

            if (remaining >= 3)
            {
                const float x2 =
                    samples[pos + 2];

                li += x2 * c2;
                lq += x2 * s2;
            }

            out_i +=
                li * gc -
                lq * gs;

            out_q +=
                li * gs +
                lq * gc;
        }
    };

    // ========================================================
    // CENTER
    // ========================================================

    const int center_delay =
        (int)lroundf(
            delay0 *
            (float)SAMPLE_RATE);

    const int first_delay =
        std::max(
            0,
            center_delay -
            COARSE_RADIUS);

    const int last_delay =
        std::min(
            num_samples - 1,
            center_delay +
            COARSE_RADIUS);

    if (last_delay < first_delay)
        return delay0;

    // ========================================================
    // COARSE SEARCH
    //
    // FULL BRUTE FORCE
    //
    // Ogni candidato viene ricalcolato da zero.
    //
    // ========================================================

    const int ncoarse =
        (last_delay - first_delay) /
        COARSE_STEP + 1;

    float score_total[MAX_COARSE] = {};
    int used_count[MAX_COARSE] = {};

    for (int d = 0;
         d < ncoarse;
         ++d)
    {
        const int delay =
            first_delay +
            d * COARSE_STEP;

        float total = 0.0f;
        int used = 0;

        // ----------------------------------------------------
        // SOLO Costas 0..6
        // ----------------------------------------------------

        for (int k = 0;
             k < NCOSTAS;
             ++k)
        {
            const int start =
                delay +
                k * SYMBOL_SAMPLES;

            if (start < 0 ||
                start + SYMBOL_SAMPLES >
                    num_samples)
            {
                continue;
            }

            const float f =
                freq +
                TONE_SPACING *
                (float)tones[k];

            const float w =
                TWO_PI * f /
                (float)SAMPLE_RATE;

            float ci;
            float cq;

            // IDENTICA correlazione V8.1

            initial_correlation(
                start,
                SYMBOL_SAMPLES,
                w,
                ci,
                cq);

            total +=
                ci * ci +
                cq * cq;

            ++used;
        }

        score_total[d] = total;
        used_count[d] = used;
    }

    // ========================================================
    // FIND COARSE MAX
    // ========================================================

    float best_coarse_score = -1.0f;
    int best_coarse_index = 0;

    for (int d = 0;
         d < ncoarse;
         ++d)
    {
        if (used_count[d] <= 0)
            continue;

        const float metric =
            score_total[d] /
            (float)used_count[d];

        if (metric > best_coarse_score)
        {
            best_coarse_score =
                metric;

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
    // FULL BRUTE FORCE
    //
    // ========================================================

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

    int nfine =
        fine_last -
        fine_first + 1;

    if (nfine > MAX_FINE)
        nfine = MAX_FINE;

    float fine_score[MAX_FINE] = {};

    int fine_used[MAX_FINE] = {};

    for (int d = 0;
         d < nfine;
         ++d)
    {
        const int delay =
            fine_first + d;

        float total = 0.0f;
        int used = 0;

        for (int k = 0;
             k < NCOSTAS;
             ++k)
        {
            const int start =
                delay +
                k * SYMBOL_SAMPLES;

            if (start < 0 ||
                start + SYMBOL_SAMPLES >
                    num_samples)
            {
                continue;
            }

            const float f =
                freq +
                TONE_SPACING *
                (float)tones[k];

            const float w =
                TWO_PI * f /
                (float)SAMPLE_RATE;

            float ci;
            float cq;

            // IDENTICA V8.1

            initial_correlation(
                start,
                SYMBOL_SAMPLES,
                w,
                ci,
                cq);

            total +=
                ci * ci +
                cq * cq;

            ++used;
        }

        fine_score[d] = total;
        fine_used[d] = used;
    }

    // ========================================================
    // FIND FINE MAX
    // ========================================================

    float best_fine_score = -1.0f;

    int best_fine_delay =
        coarse_delay;

    for (int d = 0;
         d < nfine;
         ++d)
    {
        if (fine_used[d] <= 0)
            continue;

        const float metric =
            fine_score[d] /
            (float)fine_used[d];

        if (metric > best_fine_score)
        {
            best_fine_score =
                metric;

            best_fine_delay =
                fine_first + d;
        }
    }

    // ========================================================
    // DEBUG
    // ========================================================

    if (cand_to_print > 0)
    {
        Serial.printf(
            "V10C: coarse=%d  fine=%d  "
            "delay=%.6f ms\n",
            coarse_delay,
            best_fine_delay,
            (float)best_fine_delay *
                1000.0f /
                (float)SAMPLE_RATE);
    }

    // ========================================================
    // RESULT
    // ========================================================

    return
        (float)best_fine_delay /
        (float)SAMPLE_RATE;
}


float refine_ft8_delay_v10b(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE    = 12000;
    constexpr int SYMBOL_SAMPLES = 1920;
    constexpr int COARSE_RADIUS  = 600;
    constexpr int COARSE_STEP    = 4;
    constexpr int FINE_RADIUS    = 6;

    constexpr float TONE_SPACING = 6.25f;
    constexpr int NCOSTAS        = 7;

    // First Costas block: symbols 0..6
    int costas_k[NCOSTAS];

    for (int m = 0; m < NCOSTAS; ++m)
        costas_k[m] = m;

    const int center_delay =
        (int)lroundf(delay0 * SAMPLE_RATE);

    int first_delay =
        center_delay - COARSE_RADIUS;

    int last_delay =
        center_delay + COARSE_RADIUS;

    if (first_delay < 0)
        first_delay = 0;

    if (last_delay >= num_samples)
        last_delay = num_samples - 1;


    // ------------------------------------------------------------
    // Helper: complete correlation of the 7 Costas symbols
    // ------------------------------------------------------------

    auto calculate_metric =
        [&](int delay) -> float
        {
            float total = 0.0f;
            int used = 0;

            for (int m = 0; m < NCOSTAS; ++m)
            {
                const int k =
                    costas_k[m];

                const int start =
                    delay +
                    k * SYMBOL_SAMPLES;

                // All first 7 Costas symbols must be complete.
                if (start < 0 ||
                    start + SYMBOL_SAMPLES > num_samples)
                {
                    continue;
                }

                const float f =
                    freq +
                    TONE_SPACING *
                    (float)tones[k];

                const float w =
                    TWO_PI * f /
                    (float)SAMPLE_RATE;

                const float cw = cosf(w);
                const float sw = sinf(w);

                float c = 1.0f;
                float s = 0.0f;

                float re = 0.0f;
                float im = 0.0f;

                for (int n = 0;
                     n < SYMBOL_SAMPLES;
                     ++n)
                {
                    const float x =
                        samples[start + n];

                    re += x * c;
                    im += x * s;

                    const float nc =
                        c * cw +
                        s * sw;

                    const float ns =
                        s * cw -
                        c * sw;

                    c = nc;
                    s = ns;
                }

                total +=
                    re * re +
                    im * im;

                ++used;
            }

            if (used == 0)
                return 0.0f;

            return total /
                   (float)used;
        };


    // ============================================================
    // COARSE SEARCH
    // ============================================================

    int coarse_delay =
        first_delay;

    float best_coarse_metric =
        -1.0f;

    for (int delay = first_delay;
         delay <= last_delay;
         delay += COARSE_STEP)
    {
        const float metric =
            calculate_metric(delay);

        if (metric > best_coarse_metric)
        {
            best_coarse_metric = metric;
            coarse_delay = delay;
        }
    }


    // ============================================================
    // FINE SEARCH
    // ============================================================

    int fine_first =
        coarse_delay -
        FINE_RADIUS;

    int fine_last =
        coarse_delay +
        FINE_RADIUS;

    if (fine_first < first_delay)
        fine_first = first_delay;

    if (fine_last > last_delay)
        fine_last = last_delay;

    int best_fine_delay =
        coarse_delay;

    float best_fine_metric =
        -1.0f;

    if (cand_to_print > 0)
    {
        Serial.printf(
            "\nV10B FINE SEARCH: "
            "coarse=%d  range=[%d..%d]\n",
            coarse_delay,
            fine_first,
            fine_last);
    }

    for (int delay = fine_first;
         delay <= fine_last;
         ++delay)
    {
        const float metric =
            calculate_metric(delay);

        if (cand_to_print > 0)
        {
            Serial.printf(
                "  delay=%6d  "
                "offset=%+4d  "
                "metric=%12.3f",
                delay,
                delay - coarse_delay,
                metric);
        }

        if (metric > best_fine_metric)
        {
            best_fine_metric = metric;
            best_fine_delay = delay;

            if (cand_to_print > 0)
                Serial.printf("  <--- BEST");
        }

        if (cand_to_print > 0)
            Serial.println();
    }


    // ============================================================
    // RESULT
    // ============================================================

    if (cand_to_print > 0)
    {
        Serial.printf(
            "\nV10B RESULT:\n"
            "  coarse = %d\n"
            "  fine   = %d\n"
            "  delay  = %.6f ms\n"
            "  error  = %+d samples\n",
            coarse_delay,
            best_fine_delay,
            (float)best_fine_delay *
                1000.0f /
                (float)SAMPLE_RATE,
            best_fine_delay -
                center_delay);
    }

    return (float)best_fine_delay /
           (float)SAMPLE_RATE;
}

// ============================================================
// FT8 DELAY REFINEMENT V10
//
// - COSTAS ONLY
// - 7 correlazioni indipendenti
// - ogni simbolo usa la propria frequenza:
//       freq + 6.25 * tones[k]
// - COARSE: shift diretto di 4 campioni
// - FINE: identico approccio V7, 1 sample alla volta
//
// Obiettivo:
//   stessa risposta della V7, ma meno rotazioni nel coarse search
// ============================================================

float refine_ft8_delay_v10(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE    = 12000;
    constexpr int NTONES         = 79;
    constexpr int SYMBOL_SAMPLES = 1920;       // 160 ms
    constexpr int COARSE_RADIUS  = 600;        // +/- 50 ms
    constexpr int COARSE_STEP    = 4;
    constexpr int FINE_RADIUS    = 6;

    constexpr int MAX_COARSE =
        2 * COARSE_RADIUS / COARSE_STEP + 1;   // 301

    constexpr float TONE_SPACING = 6.25f;

    // ========================================================
    // V10
    //
    // SOLO IL PRIMO BLOCCO COSTAS:
    //
    //   0..6
    //
    // In totale: 7 simboli
    // ========================================================

    constexpr int NCOSTAS = 7;

    int costas_k[NCOSTAS];

    for (int m = 0; m < NCOSTAS; ++m)
        costas_k[m] = m;

    // ========================================================
    // Centro della ricerca
    //
    // delay0 è espresso in secondi
    // ========================================================

    const int center_delay =
        (int)lroundf(delay0 * SAMPLE_RATE);

    int first_delay =
        center_delay - COARSE_RADIUS;

    int last_delay =
        center_delay + COARSE_RADIUS;

    if (first_delay < 0)
        first_delay = 0;

    if (last_delay >= num_samples)
        last_delay = num_samples - 1;

    const int ncoarse =
        (last_delay - first_delay) /
        COARSE_STEP + 1;

    // ========================================================
    // Score globale
    // ========================================================

    float score_total[MAX_COARSE] = {};
    uint16_t used_count[MAX_COARSE] = {};

    // ========================================================
    // Correlazioni indipendenti dei 7 Costas
    // ========================================================

    float ci[NCOSTAS];
    float cq[NCOSTAS];

    float pc[NCOSTAS];
    float ps[NCOSTAS];

    // P^1 ... P^4

    float p1c[NCOSTAS], p1s[NCOSTAS];
    float p2c[NCOSTAS], p2s[NCOSTAS];
    float p3c[NCOSTAS], p3s[NCOSTAS];
    float p4c[NCOSTAS], p4s[NCOSTAS];

    // E*P^1 ... E*P^4

    float ep1c[NCOSTAS], ep1s[NCOSTAS];
    float ep2c[NCOSTAS], ep2s[NCOSTAS];
    float ep3c[NCOSTAS], ep3s[NCOSTAS];
    float ep4c[NCOSTAS], ep4s[NCOSTAS];

    // ========================================================
    // INITIAL CORRELATION
    //
    // Ogni simbolo ha la propria frequenza:
    //
    //   freq + 6.25 * tones[k]
    //
    // ========================================================

    for (int m = 0; m < NCOSTAS; ++m)
    {
        const int k = costas_k[m];

        const float f =
            freq +
            TONE_SPACING * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        const float cw = cosf(w);
        const float sw = sinf(w);

        pc[m] = cw;
        ps[m] = sw;

        // ----------------------------------------------------
        // P = exp(+j*w)
        // ----------------------------------------------------

        p1c[m] = cw;
        p1s[m] = sw;

        // ----------------------------------------------------
        // P^2
        // ----------------------------------------------------

        p2c[m] =
            cw * cw -
            sw * sw;

        p2s[m] =
            2.0f * cw * sw;

        // ----------------------------------------------------
        // P^3
        // ----------------------------------------------------

        p3c[m] =
            p2c[m] * cw -
            p2s[m] * sw;

        p3s[m] =
            p2c[m] * sw +
            p2s[m] * cw;

        // ----------------------------------------------------
        // P^4
        // ----------------------------------------------------

        p4c[m] =
            p3c[m] * cw -
            p3s[m] * sw;

        p4s[m] =
            p3c[m] * sw +
            p3s[m] * cw;

        // ----------------------------------------------------
        // E = exp(-j*w*SYMBOL_SAMPLES)
        // ----------------------------------------------------

        const float phase =
            w * (float)SYMBOL_SAMPLES;

        const float ec = cosf(phase);
        const float es = -sinf(phase);

        // ----------------------------------------------------
        // E * P
        // ----------------------------------------------------

        ep1c[m] =
            ec * p1c[m] -
            es * p1s[m];

        ep1s[m] =
            ec * p1s[m] +
            es * p1c[m];

        // ----------------------------------------------------
        // E * P^2
        // ----------------------------------------------------

        ep2c[m] =
            ec * p2c[m] -
            es * p2s[m];

        ep2s[m] =
            ec * p2s[m] +
            es * p2c[m];

        // ----------------------------------------------------
        // E * P^3
        // ----------------------------------------------------

        ep3c[m] =
            ec * p3c[m] -
            es * p3s[m];

        ep3s[m] =
            ec * p3s[m] +
            es * p3c[m];

        // ----------------------------------------------------
        // E * P^4
        // ----------------------------------------------------

        ep4c[m] =
            ec * p4c[m] -
            es * p4s[m];

        ep4s[m] =
            ec * p4s[m] +
            es * p4c[m];

        // ----------------------------------------------------
        // Correlazione iniziale
        //
        // C = sum x[n] * exp(-j*w*n)
        // ----------------------------------------------------

        const int start =
            first_delay +
            k * SYMBOL_SAMPLES;

        const int count =
            SYMBOL_SAMPLES;

        if (start < 0 ||
            start + count > num_samples)
        {
            ci[m] = 0.0f;
            cq[m] = 0.0f;
            continue;
        }

        float cr = 0.0f;
        float ci_local = 0.0f;

        float c = 1.0f;
        float s = 0.0f;

        for (int n = 0; n < count; ++n)
        {
            const float x =
                samples[start + n];

            cr += x * c;
            ci_local += x * s;

            // exp(-j*w*(n+1))

            const float nc =
                c * cw +
                s * sw;

            const float ns =
                s * cw -
                c * sw;

            c = nc;
            s = ns;
        }

        ci[m] = cr;
        cq[m] = ci_local;
    }

    // ========================================================
    // COARSE SEARCH
    //
    // Direct update di 4 campioni.
    //
    // Equivalente a 4 aggiornamenti V7 consecutivi.
    // ========================================================

    int coarse_delay = first_delay;

    float best_score = -1.0f;

    for (int d = 0; d < ncoarse; ++d)
    {
        // ----------------------------------------------------
        // Score corrente
        // ----------------------------------------------------

        float total = 0.0f;
        int used = 0;

        for (int m = 0; m < NCOSTAS; ++m)
        {
            total +=
                ci[m] * ci[m] +
                cq[m] * cq[m];

            ++used;
        }

        score_total[d] = total;
        used_count[d] = used;

        // ----------------------------------------------------
        // Normalizzazione
        // ----------------------------------------------------

        float metric = 0.0f;

        if (used > 0)
            metric =
                total / (float)used;

        if (metric > best_score)
        {
            best_score = metric;

            coarse_delay =
                first_delay +
                d * COARSE_STEP;
        }

        // ----------------------------------------------------
        // Ultimo candidato
        // ----------------------------------------------------

        if (d == ncoarse - 1)
            break;

        // ----------------------------------------------------
        // SHIFT +4
        // ----------------------------------------------------

        for (int m = 0; m < NCOSTAS; ++m)
        {
            const int k =
                costas_k[m];

            const int symbol_start =
                first_delay +
                k * SYMBOL_SAMPLES;

            const int offset =
                d * COARSE_STEP;

            const int old_pos =
                symbol_start +
                offset;

            const int new_pos =
                old_pos +
                SYMBOL_SAMPLES;

            // ------------------------------------------------
            // Campioni uscenti
            // ------------------------------------------------

            const float x0 =
                samples[old_pos + 0];

            const float x1 =
                samples[old_pos + 1];

            const float x2 =
                samples[old_pos + 2];

            const float x3 =
                samples[old_pos + 3];

            // ------------------------------------------------
            // Campioni entranti
            // ------------------------------------------------

            const float y0 =
                samples[new_pos + 0];

            const float y1 =
                samples[new_pos + 1];

            const float y2 =
                samples[new_pos + 2];

            const float y3 =
                samples[new_pos + 3];

            // ------------------------------------------------
            // C0 * P4
            // ------------------------------------------------

            float nr =
                ci[m] * p4c[m] -
                cq[m] * p4s[m];

            float nq =
                ci[m] * p4s[m] +
                cq[m] * p4c[m];

            // ------------------------------------------------
            // -x0 * P4
            // ------------------------------------------------

            nr -= x0 * p4c[m];
            nq -= x0 * p4s[m];

            // ------------------------------------------------
            // -x1 * P3
            // ------------------------------------------------

            nr -= x1 * p3c[m];
            nq -= x1 * p3s[m];

            // ------------------------------------------------
            // -x2 * P2
            // ------------------------------------------------

            nr -= x2 * p2c[m];
            nq -= x2 * p2s[m];

            // ------------------------------------------------
            // -x3 * P
            // ------------------------------------------------

            nr -= x3 * p1c[m];
            nq -= x3 * p1s[m];

            // ------------------------------------------------
            // +y0 * E*P4
            // ------------------------------------------------

            nr += y0 * ep4c[m];
            nq += y0 * ep4s[m];

            // ------------------------------------------------
            // +y1 * E*P3
            // ------------------------------------------------

            nr += y1 * ep3c[m];
            nq += y1 * ep3s[m];

            // ------------------------------------------------
            // +y2 * E*P2
            // ------------------------------------------------

            nr += y2 * ep2c[m];
            nq += y2 * ep2s[m];

            // ------------------------------------------------
            // +y3 * E*P
            // ------------------------------------------------

            nr += y3 * ep1c[m];
            nq += y3 * ep1s[m];

            ci[m] = nr;
            cq[m] = nq;
        }
    }

    // ========================================================
    // FINE SEARCH
    //
    // +/- 6 campioni attorno al risultato coarse.
    // ========================================================

    int fine_first =
        coarse_delay -
        FINE_RADIUS;

    int fine_last =
        coarse_delay +
        FINE_RADIUS;

    if (fine_first < first_delay)
        fine_first = first_delay;

    if (fine_last > last_delay)
        fine_last = last_delay;

    constexpr int MAX_FINE = 16;

    float accumulated_fine[MAX_FINE] = {};
    uint16_t fine_count[MAX_FINE] = {};

    int nfine =
        fine_last -
        fine_first +
        1;

    if (nfine > MAX_FINE)
        nfine = MAX_FINE;

    // ========================================================
    // FINE SEARCH
    //
    // Un Costas alla volta, ma solamente i primi 7.
    // ========================================================

    for (int m = 0; m < NCOSTAS; ++m)
    {
        const int k =
            costas_k[m];

        const float f =
            freq +
            TONE_SPACING * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        const float cw = cosf(w);
        const float sw = sinf(w);

        const int start =
            fine_first +
            k * SYMBOL_SAMPLES;

        if (start < 0 ||
            start + SYMBOL_SAMPLES > num_samples)
        {
            continue;
        }

        float cr = 0.0f;
        float cq_local = 0.0f;

        float c = 1.0f;
        float s = 0.0f;

        // ----------------------------------------------------
        // Correlazione iniziale
        // ----------------------------------------------------

        for (int n = 0;
             n < SYMBOL_SAMPLES;
             ++n)
        {
            const float x =
                samples[start + n];

            cr += x * c;
            cq_local += x * s;

            const float nc =
                c * cw +
                s * sw;

            const float ns =
                s * cw -
                c * sw;

            c = nc;
            s = ns;
        }

        // ----------------------------------------------------
        // Shift di 1 campione
        // ----------------------------------------------------

        for (int d = 0;
             d < nfine;
             ++d)
        {
            accumulated_fine[d] +=
                cr * cr +
                cq_local * cq_local;

            fine_count[d]++;

            if (d == nfine - 1)
                break;

            const int old_pos =
                start + d;

            const int new_pos =
                old_pos +
                SYMBOL_SAMPLES;

            const float x_old =
                samples[old_pos];

            const float x_new =
                samples[new_pos];

            // ------------------------------------------------
            // V7 one-sample update
            // ------------------------------------------------

            const float ti =
                cr -
                x_old +
                x_new * c;

            const float tq =
                cq_local +
                x_new * s;

            const float nr =
                ti * cw -
                tq * sw;

            const float nq =
                ti * sw +
                tq * cw;

            cr = nr;
            cq_local = nq;
        }
    }

    // ========================================================
    // Miglior punto fine
    // ========================================================

    int best_fine_delay =
        coarse_delay;

    float best_fine_metric =
        -1.0f;

    for (int d = 0;
         d < nfine;
         ++d)
    {
        if (fine_count[d] == 0)
            continue;

        const float metric =
            accumulated_fine[d] /
            (float)fine_count[d];

        if (metric > best_fine_metric)
        {
            best_fine_metric = metric;

            best_fine_delay =
                fine_first + d;
        }
    }

    // ========================================================
    // Risultato
    // ========================================================

        if (cand_to_print > 0)
    {
        Serial.printf(
            "\nV10 FINE SEARCH: coarse=%d  "
            "range=[%d .. %d]\n",
            coarse_delay,
            fine_first,
            fine_last);

        for (int d = 0; d < nfine; ++d)
        {
            const int delay =
                fine_first + d;

            const float metric =
                (fine_count[d] > 0)
                    ? accumulated_fine[d] /
                      (float)fine_count[d]
                    : 0.0f;

            Serial.printf(
                "  delay=%6d  "
                "offset=%+4d  "
                "metric=%12.3f%s\n",
                delay,
                delay - coarse_delay,
                metric,
                (delay == best_fine_delay)
                    ? "  <--- BEST"
                    : "");
        }

        Serial.printf(
            "V10 RESULT: coarse=%d  fine=%d  "
            "delay=%.6f ms\n",
            coarse_delay,
            best_fine_delay,
            (float)best_fine_delay *
                1000.0f /
                (float)SAMPLE_RATE);
    }

    return (float)best_fine_delay /
           (float)SAMPLE_RATE;
}

// ============================================================
// FT8 DELAY REFINEMENT V9
//
// - COSTAS ONLY
// - 21 correlazioni indipendenti
// - ogni simbolo usa la propria frequenza:
//       freq + 6.25 * tones[k]
// - COARSE: shift diretto di 4 campioni
// - FINE: identico approccio V7, 1 sample alla volta
//
// Obiettivo:
//   stessa risposta della V7, ma meno rotazioni nel coarse search
// ============================================================

float refine_ft8_delay_v9(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE    = 12000;
    constexpr int NTONES         = 79;
    constexpr int SYMBOL_SAMPLES = 1920;       // 160 ms
    constexpr int COARSE_RADIUS  = 600;        // +/- 50 ms
    constexpr int COARSE_STEP    = 4;
    constexpr int FINE_RADIUS    = 6;

    constexpr int MAX_COARSE =
        2 * COARSE_RADIUS / COARSE_STEP + 1;   // 301

    constexpr float TONE_SPACING = 6.25f;

    // --------------------------------------------------------
    // Costas symbols:
    //
    //   0..6
    //   36..42
    //   72..78
    // --------------------------------------------------------

    auto is_costas = [](int k) -> bool
    {
        return (k < 7) ||
               (k >= 36 && k < 43) ||
               (k >= 72 && k < 79);
    };

    // --------------------------------------------------------
    // Elenco dei 21 simboli Costas
    // --------------------------------------------------------

    int costas_k[21];
    int ncostas = 0;

    for (int k = 0; k < NTONES; ++k)
    {
        if (is_costas(k))
            costas_k[ncostas++] = k;
    }

    // --------------------------------------------------------
    // Centro della ricerca
    // --------------------------------------------------------

    int center_delay = (int)lroundf(delay0 * SAMPLE_RATE);

    int first_delay = center_delay - COARSE_RADIUS;
    int last_delay  = center_delay + COARSE_RADIUS;

    if (first_delay < 0)
        first_delay = 0;

    if (last_delay >= num_samples)
        last_delay = num_samples - 1;

    int ncoarse =
        (last_delay - first_delay) / COARSE_STEP + 1;

    if (ncoarse > MAX_COARSE)
        ncoarse = MAX_COARSE;

    // --------------------------------------------------------
    // Score globale del coarse search
    // --------------------------------------------------------

    float score_total[MAX_COARSE] = {};
    uint16_t used_count[MAX_COARSE] = {};

    // --------------------------------------------------------
    // Correlazioni indipendenti dei 21 simboli Costas
    //
    // ci[k] + j*cq[k]
    // --------------------------------------------------------

    float ci[21];
    float cq[21];

    float pc[21];
    float ps[21];

    // P  = exp(+j*w)
    // E  = exp(-j*w*SYMBOL_SAMPLES)
    //
    // Per lo shift diretto di 4 campioni servono:
    //
    // P
    // P^2
    // P^3
    // P^4
    //
    // e:
    //
    // E*P
    // E*P^2
    // E*P^3
    // E*P^4

    float p1c[21], p1s[21];
    float p2c[21], p2s[21];
    float p3c[21], p3s[21];
    float p4c[21], p4s[21];

    float ep1c[21], ep1s[21];
    float ep2c[21], ep2s[21];
    float ep3c[21], ep3s[21];
    float ep4c[21], ep4s[21];

    // --------------------------------------------------------
    // Initial correlation
    //
    // IMPORTANTE:
    // eseguita separatamente per ogni simbolo Costas.
    //
    // Questo mantiene la stessa semantica della V7.
    // --------------------------------------------------------

    for (int m = 0; m < ncostas; ++m)
    {
        const int k = costas_k[m];

        const float f =
            freq + TONE_SPACING * (float)tones[k];

        const float w =
            TWO_PI * f / (float)SAMPLE_RATE;

        const float cw = cosf(w);
        const float sw = sinf(w);

        pc[m] = cw;
        ps[m] = sw;

        // ----------------------------------------------------
        // P = exp(+j*w)
        // ----------------------------------------------------

        p1c[m] = cw;
        p1s[m] = sw;

        // P2
        p2c[m] = cw * cw - sw * sw;
        p2s[m] = 2.0f * cw * sw;

        // P3 = P2 * P
        p3c[m] =
            p2c[m] * cw -
            p2s[m] * sw;

        p3s[m] =
            p2c[m] * sw +
            p2s[m] * cw;

        // P4 = P3 * P
        p4c[m] =
            p3c[m] * cw -
            p3s[m] * sw;

        p4s[m] =
            p3c[m] * sw +
            p3s[m] * cw;

        // ----------------------------------------------------
        // E = exp(-j*w*SYMBOL_SAMPLES)
        //
        // Per FT8, 1920 campioni = 160 ms.
        // ----------------------------------------------------

        const float phase =
            w * (float)SYMBOL_SAMPLES;

        const float ec = cosf(phase);
        const float es = -sinf(phase);

        // E * P
        ep1c[m] =
            ec * p1c[m] -
            es * p1s[m];

        ep1s[m] =
            ec * p1s[m] +
            es * p1c[m];

        // E * P2
        ep2c[m] =
            ec * p2c[m] -
            es * p2s[m];

        ep2s[m] =
            ec * p2s[m] +
            es * p2c[m];

        // E * P3
        ep3c[m] =
            ec * p3c[m] -
            es * p3s[m];

        ep3s[m] =
            ec * p3s[m] +
            es * p3c[m];

        // E * P4
        ep4c[m] =
            ec * p4c[m] -
            es * p4s[m];

        ep4s[m] =
            ec * p4s[m] +
            es * p4c[m];

        // ----------------------------------------------------
        // Initial correlation
        //
        // C = sum x[n] * exp(-j*w*n)
        // ----------------------------------------------------

        const int start =
            first_delay + k * SYMBOL_SAMPLES;

        const int count =
            SYMBOL_SAMPLES;

        if (start < 0 ||
            start + count > num_samples)
        {
            ci[m] = 0.0f;
            cq[m] = 0.0f;
            continue;
        }

        float cr = 0.0f;
        float ci_local = 0.0f;

        float c = 1.0f;
        float s = 0.0f;

        for (int n = 0; n < count; ++n)
        {
            const float x = samples[start + n];

            cr += x * c;
            ci_local += x * s;

            // exp(-j*w*(n+1))
            const float nc =
                c * cw + s * sw;

            const float ns =
                s * cw - c * sw;

            c = nc;
            s = ns;
        }

        ci[m] = cr;
        cq[m] = ci_local;
    }

    // ========================================================
    // COARSE SEARCH
    //
    // V7:
    //   4 volte:
    //      remove 1 sample
    //      add 1 sample
    //      rotate C
    //
    // V9:
    //   direttamente:
    //
    //   C4 = C0*P4
    //        - x0*P4
    //        - x1*P3
    //        - x2*P2
    //        - x3*P
    //        + y0*E*P4
    //        + y1*E*P3
    //        + y2*E*P2
    //        + y3*E*P
    //
    // dove:
    //
    //   x0..x3 = campioni che escono
    //   y0..y3 = campioni che entrano
    // ========================================================

    int coarse_delay = first_delay;
    float best_score = -1.0f;

    for (int d = 0; d < ncoarse; ++d)
    {
        // ----------------------------------------------------
        // 1. Score corrente
        // ----------------------------------------------------

        float total = 0.0f;
        int used = 0;

        for (int m = 0; m < ncostas; ++m)
        {
            total +=
                ci[m] * ci[m] +
                cq[m] * cq[m];

            ++used;
        }

        score_total[d] = total;
        used_count[d] = used;

        // ----------------------------------------------------
        // 2. Verifica best
        // ----------------------------------------------------

        float metric = 0.0f;

        if (used > 0)
            metric = total / (float)used;

        if (metric > best_score)
        {
            best_score = metric;
            coarse_delay =
                first_delay + d * COARSE_STEP;
        }

        // ----------------------------------------------------
        // 3. Ultimo candidato?
        // ----------------------------------------------------

        if (d == ncoarse - 1)
            break;

        // ----------------------------------------------------
        // 4. Shift di +4 campioni
        //
        // Ogni Costas symbol è completamente indipendente.
        // ----------------------------------------------------

        for (int m = 0; m < ncostas; ++m)
        {
            const int k = costas_k[m];

            const int symbol_start =
                first_delay +
                k * SYMBOL_SAMPLES;

            const int offset =
                d * COARSE_STEP;

            const int old_pos =
                symbol_start + offset;

            const int new_pos =
                old_pos + SYMBOL_SAMPLES;

            // ------------------------------------------------
            // Campioni che escono
            // ------------------------------------------------

            const float x0 = samples[old_pos + 0];
            const float x1 = samples[old_pos + 1];
            const float x2 = samples[old_pos + 2];
            const float x3 = samples[old_pos + 3];

            // ------------------------------------------------
            // Campioni che entrano
            // ------------------------------------------------

            const float y0 = samples[new_pos + 0];
            const float y1 = samples[new_pos + 1];
            const float y2 = samples[new_pos + 2];
            const float y3 = samples[new_pos + 3];

            // ------------------------------------------------
            // C0 * P4
            // ------------------------------------------------

            float nr =
                ci[m] * p4c[m] -
                cq[m] * p4s[m];

            float nq =
                ci[m] * p4s[m] +
                cq[m] * p4c[m];

            // ------------------------------------------------
            // -x0 * P4
            // ------------------------------------------------

            nr -= x0 * p4c[m];
            nq -= x0 * p4s[m];

            // ------------------------------------------------
            // -x1 * P3
            // ------------------------------------------------

            nr -= x1 * p3c[m];
            nq -= x1 * p3s[m];

            // ------------------------------------------------
            // -x2 * P2
            // ------------------------------------------------

            nr -= x2 * p2c[m];
            nq -= x2 * p2s[m];

            // ------------------------------------------------
            // -x3 * P
            // ------------------------------------------------

            nr -= x3 * p1c[m];
            nq -= x3 * p1s[m];

            // ------------------------------------------------
            // +y0 * E*P4
            // ------------------------------------------------

            nr += y0 * ep4c[m];
            nq += y0 * ep4s[m];

            // ------------------------------------------------
            // +y1 * E*P3
            // ------------------------------------------------

            nr += y1 * ep3c[m];
            nq += y1 * ep3s[m];

            // ------------------------------------------------
            // +y2 * E*P2
            // ------------------------------------------------

            nr += y2 * ep2c[m];
            nq += y2 * ep2s[m];

            // ------------------------------------------------
            // +y3 * E*P
            // ------------------------------------------------

            nr += y3 * ep1c[m];
            nq += y3 * ep1s[m];

            ci[m] = nr;
            cq[m] = nq;
        }
    }

    // ========================================================
    // FINE SEARCH
    //
    // Manteniamo volutamente la logica V7:
    //
    //   +/- 6 samples attorno al coarse result
    //
    // partendo da:
    //
    //   coarse_delay - 6
    //
    // e avanzando di 1 sample.
    //
    // In questo modo la V9 può essere confrontata
    // direttamente con la V7.
    // ========================================================

    int fine_first =
        coarse_delay - FINE_RADIUS;

    int fine_last =
        coarse_delay + FINE_RADIUS;

    if (fine_first < first_delay)
        fine_first = first_delay;

    if (fine_last > last_delay)
        fine_last = last_delay;

    constexpr int MAX_FINE = 16;

    float accumulated_fine[MAX_FINE] = {};
    uint16_t fine_count[MAX_FINE] = {};

    int nfine =
        fine_last - fine_first + 1;

    if (nfine > MAX_FINE)
        nfine = MAX_FINE;

    // --------------------------------------------------------
    // Fine search: inizializzazione separata per ogni simbolo
    // --------------------------------------------------------

    for (int m = 0; m < ncostas; ++m)
    {
        const int k = costas_k[m];

        const float f =
            freq + TONE_SPACING * (float)tones[k];

        const float w =
            TWO_PI * f / (float)SAMPLE_RATE;

        const float cw = cosf(w);
        const float sw = sinf(w);

        const int start =
            fine_first + k * SYMBOL_SAMPLES;

        if (start < 0 ||
            start + SYMBOL_SAMPLES > num_samples)
        {
            continue;
        }

        float cr = 0.0f;
        float cq_local = 0.0f;

        float c = 1.0f;
        float s = 0.0f;

        // ----------------------------------------------------
        // Correlazione iniziale
        // ----------------------------------------------------

        for (int n = 0; n < SYMBOL_SAMPLES; ++n)
        {
            const float x =
                samples[start + n];

            cr += x * c;
            cq_local += x * s;

            const float nc =
                c * cw + s * sw;

            const float ns =
                s * cw - c * sw;

            c = nc;
            s = ns;
        }

        // ----------------------------------------------------
        // Scansione fine, 1 sample alla volta
        // ----------------------------------------------------

        for (int d = 0; d < nfine; ++d)
        {
            accumulated_fine[d] +=
                cr * cr +
                cq_local * cq_local;

            fine_count[d]++;

            if (d == nfine - 1)
                break;

            const int old_pos =
                start + d;

            const int new_pos =
                old_pos + SYMBOL_SAMPLES;

            const float x_old =
                samples[old_pos];

            const float x_new =
                samples[new_pos];

            // ------------------------------------------------
            // V7 exact one-sample update
            // ------------------------------------------------

            const float ti =
                cr - x_old + x_new * c;

            const float tq =
                cq_local + x_new * s;

            const float nr =
                ti * cw - tq * sw;

            const float nq =
                ti * sw + tq * cw;

            cr = nr;
            cq_local = nq;
        }
    }

    // ========================================================
    // Cerca il miglior punto fine
    // ========================================================

    int best_fine_delay =
        coarse_delay;

    float best_fine_metric =
        -1.0f;

    for (int d = 0; d < nfine; ++d)
    {
        if (fine_count[d] == 0)
            continue;

        const float metric =
            accumulated_fine[d] /
            (float)fine_count[d];

        if (metric > best_fine_metric)
        {
            best_fine_metric = metric;

            best_fine_delay =
                fine_first + d;
        }
    }

    // ========================================================
    // Risultato
    // ========================================================

    const float result =
        (float)best_fine_delay;

    // --------------------------------------------------------
    // Debug opzionale
    // --------------------------------------------------------

    if (cand_to_print > 0)
    {
        Serial.printf(
            "V9: coarse=%d  fine=%d  "
            "delay=%.6f ms\n",
            coarse_delay,
            best_fine_delay,
            result * 1000.0f /
                (float)SAMPLE_RATE
        );
    }

    return result /
           (float)SAMPLE_RATE;
}

// ============================================================
// FT8 DELAY REFINEMENT V8.2
//
// - Costas only
// - Initial correlation: optimized group-of-4
// - Coarse search: 4-sample direct sliding update
// - Fine search: 1-sample sliding update
//
// IMPORTANT:
// The correlation is:
//
//     C(d) = sum x[d+n] * exp(-j*w*n)
//
// When the window moves forward by one sample:
//
//     C(d+1) = (C(d) - x_old + x_new*E) * exp(+j*w)
//
// Therefore the sliding rotation is +w.
// ============================================================

float refine_ft8_delay_v8_2(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE   = 12000;
    constexpr int NTONES        = 79;
    constexpr int SYMBOL_SAMPLES = 1920;

    constexpr int COARSE_RADIUS = 600;
    constexpr int COARSE_STEP   = 4;
    constexpr int MAX_COARSE    = 301;

    constexpr int FINE_RADIUS   = 6;
    constexpr int MAX_FINE      = 13;

    // --------------------------------------------------------
    // Costas symbols only
    // --------------------------------------------------------

    auto is_costas =
        [](int k) -> bool
        {
            return
                (k <= 6) ||
                (k >= 36 && k <= 42) ||
                (k >= 72 && k <= 78);
        };

    // --------------------------------------------------------
    // Angular frequency
    // --------------------------------------------------------

    const float w =
        2.0f * PI * freq /
        (float)SAMPLE_RATE;

    // --------------------------------------------------------
    // Initial correlation
    //
    // C = sum x[n] exp(-j*w*n)
    //
    // Process four samples at a time.
    // --------------------------------------------------------

    auto initial_correlation =
        [&](int start,
            float& out_i,
            float& out_q)
        {
            float ci = 0.0f;
            float cq = 0.0f;

            const float cw = cosf(w);
            const float sw = sinf(w);

            // exp(-j*w)
            const float c1 = cw;
            const float s1 = -sw;

            // exp(-j*2w)
            const float c2 =
                cw * cw - sw * sw;

            const float s2 =
                -2.0f * sw * cw;

            // exp(-j*3w)
            const float c3 =
                c2 * cw - s2 * sw;

            const float s3 =
                s2 * cw + c2 * (-sw);

            // exp(-j*4w)
            const float c4 =
                c3 * cw - s3 * sw;

            const float s4 =
                s3 * cw + c3 * (-sw);

            // Group phase:
            //
            // exp(-j*4*w*g)

            float gc = 1.0f;
            float gs = 0.0f;

            int n = 0;

            // ------------------------------------------------
            // Full 4-sample groups
            // ------------------------------------------------

            while (n + 3 < NTONES * SYMBOL_SAMPLES)
            {
                const int absolute =
                    start + n;

                if (absolute + 3 >= num_samples)
                    break;

                const int symbol =
                    n / SYMBOL_SAMPLES;

                if (!is_costas(symbol))
                {
                    n += 4;
                    continue;
                }

                const float x0 =
                    samples[absolute];

                const float x1 =
                    samples[absolute + 1];

                const float x2 =
                    samples[absolute + 2];

                const float x3 =
                    samples[absolute + 3];

                // Local correlation:
                //
                // x0
                // + x1*exp(-jw)
                // + x2*exp(-j2w)
                // + x3*exp(-j3w)

                const float li =
                    x0 +
                    x1 * c1 +
                    x2 * c2 +
                    x3 * c3;

                const float lq =
                    x1 * s1 +
                    x2 * s2 +
                    x3 * s3;

                // Multiply by group phase
                //
                // (li+j*lq)*(gc+j*gs)

                ci +=
                    li * gc -
                    lq * gs;

                cq +=
                    li * gs +
                    lq * gc;

                // Advance by four samples:
                //
                // exp(-j4w)

                const float ngc =
                    gc * c4 -
                    gs * s4;

                const float ngs =
                    gc * s4 +
                    gs * c4;

                gc = ngc;
                gs = ngs;

                n += 4;
            }

            // ------------------------------------------------
            // Remaining samples
            // ------------------------------------------------

            while (n < NTONES * SYMBOL_SAMPLES)
            {
                const int absolute =
                    start + n;

                if (absolute >= num_samples)
                    break;

                const int symbol =
                    n / SYMBOL_SAMPLES;

                if (is_costas(symbol))
                {
                    // exp(-j*w*n)
                    //
                    // Current group phase corresponds to
                    // exp(-j*w*n) only if n is at group start.
                    //
                    // For the normally aligned case this loop
                    // handles only the tail.

                    const int r =
                        n & 3;

                    float rc;
                    float rs;

                    if (r == 0)
                    {
                        rc = gc;
                        rs = gs;
                    }
                    else if (r == 1)
                    {
                        rc =
                            gc * c1 -
                            gs * s1;

                        rs =
                            gc * s1 +
                            gs * c1;
                    }
                    else if (r == 2)
                    {
                        rc =
                            gc * c2 -
                            gs * s2;

                        rs =
                            gc * s2 +
                            gs * c2;
                    }
                    else
                    {
                        rc =
                            gc * c3 -
                            gs * s3;

                        rs =
                            gc * s3 +
                            gs * c3;
                    }

                    const float x =
                        samples[absolute];

                    ci += x * rc;
                    cq += x * rs;
                }

                ++n;
            }

            out_i = ci;
            out_q = cq;
        };

    // --------------------------------------------------------
    // Initial delay
    // --------------------------------------------------------

    const int delay0_samples =
        (int)lroundf(
            delay0 * SAMPLE_RATE);

    // --------------------------------------------------------
    // Initial correlation
    // --------------------------------------------------------

    float ci = 0.0f;
    float cq = 0.0f;

    initial_correlation(
        delay0_samples - COARSE_RADIUS,
        ci,
        cq);

    // --------------------------------------------------------
    // Phase factors
    //
    // E = exp(-j*w*SYMBOL_SAMPLES)
    //
    // P = exp(+jw)
    //
    // P^2, P^3, P^4
    // --------------------------------------------------------

    const float end_angle =
        -w * SYMBOL_SAMPLES;

    const float ec =
        cosf(end_angle);

    const float es =
        sinf(end_angle);

    // P = exp(+jw)

    const float pc =
        cosf(w);

    const float ps =
        sinf(w);

    // P2 = exp(+j2w)

    const float p2c =
        pc * pc -
        ps * ps;

    const float p2s =
        2.0f * pc * ps;

    // P3 = exp(+j3w)

    const float p3c =
        p2c * pc -
        p2s * ps;

    const float p3s =
        p2c * ps +
        p2s * pc;

    // P4 = exp(+j4w)

    const float p4c =
        p3c * pc -
        p3s * ps;

    const float p4s =
        p3c * ps +
        p3s * pc;

    // --------------------------------------------------------
    // E * P
    // E * P2
    // E * P3
    // E * P4
    //
    // These are the phase factors of the four new samples.
    // --------------------------------------------------------

    const float ep1c =
        ec * pc -
        es * ps;

    const float ep1s =
        ec * ps +
        es * pc;

    const float ep2c =
        ec * p2c -
        es * p2s;

    const float ep2s =
        ec * p2s +
        es * p2c;

    const float ep3c =
        ec * p3c -
        es * p3s;

    const float ep3s =
        ec * p3s +
        es * p3c;

    const float ep4c =
        ec * p4c -
        es * p4s;

    const float ep4s =
        ec * p4s +
        es * p4c;

    // --------------------------------------------------------
    // Helper: magnitude squared
    // --------------------------------------------------------

    auto metric =
        [](float i, float q) -> float
        {
            return i * i + q * q;
        };

    // --------------------------------------------------------
    // COARSE SEARCH
    //
    // One update moves four samples.
    //
    // Starting:
    //
    //     C0 = sum x[start+n] exp(-jwn)
    //
    // After 4 samples:
    //
    //     C4 =
    //       C0 * P4
    //
    //       - x0 * P4
    //       - x1 * P3
    //       - x2 * P2
    //       - x3 * P
    //
    //       + x0new * E*P4
    //       + x1new * E*P3
    //       + x2new * E*P2
    //       + x3new * E*P
    //
    // where P = exp(+jw).
    // --------------------------------------------------------

    float best_metric =
        metric(ci, cq);

    int best_coarse =
        0;

    int current_start =
        delay0_samples -
        COARSE_RADIUS;

    int best_position =
        0;

    for (int step = 1;
         step <= MAX_COARSE - 1;
         ++step)
    {
        // ----------------------------------------------------
        // Four samples leaving the left side
        // ----------------------------------------------------

        const int old0 =
            current_start;

        const int old1 =
            current_start + 1;

        const int old2 =
            current_start + 2;

        const int old3 =
            current_start + 3;

        // ----------------------------------------------------
        // Four samples entering at the right side
        //
        // Each is one full symbol later.
        // ----------------------------------------------------

        const int new0 =
            old0 + SYMBOL_SAMPLES;

        const int new1 =
            old1 + SYMBOL_SAMPLES;

        const int new2 =
            old2 + SYMBOL_SAMPLES;

        const int new3 =
            old3 + SYMBOL_SAMPLES;

        // ----------------------------------------------------
        // Read samples
        // ----------------------------------------------------

        const float x_old0 =
            samples[old0];

        const float x_old1 =
            samples[old1];

        const float x_old2 =
            samples[old2];

        const float x_old3 =
            samples[old3];

        const float x_new0 =
            samples[new0];

        const float x_new1 =
            samples[new1];

        const float x_new2 =
            samples[new2];

        const float x_new3 =
            samples[new3];

        // ----------------------------------------------------
        // Rotate old correlation by P4
        // ----------------------------------------------------

        float ti =
            ci * p4c -
            cq * p4s;

        float tq =
            ci * p4s +
            cq * p4c;

        // ----------------------------------------------------
        // Remove four old samples
        //
        // old0 has phase P4
        // old1 has phase P3
        // old2 has phase P2
        // old3 has phase P1
        // ----------------------------------------------------

        ti -=
            x_old0 * p4c +
            x_old1 * p3c +
            x_old2 * p2c +
            x_old3 * pc;

        tq -=
            x_old0 * p4s +
            x_old1 * p3s +
            x_old2 * p2s +
            x_old3 * ps;

        // ----------------------------------------------------
        // Add four new samples
        //
        // new0 -> E*P4
        // new1 -> E*P3
        // new2 -> E*P2
        // new3 -> E*P
        // ----------------------------------------------------

        ti +=
            x_new0 * ep4c +
            x_new1 * ep3c +
            x_new2 * ep2c +
            x_new3 * ep1c;

        tq +=
            x_new0 * ep4s +
            x_new1 * ep3s +
            x_new2 * ep2s +
            x_new3 * ep1s;

        ci = ti;
        cq = tq;

        current_start +=
            COARSE_STEP;

        const float m =
            metric(ci, cq);

        if (m > best_metric)
        {
            best_metric = m;
            best_position = step;
        }
    }

    // --------------------------------------------------------
    // Coarse delay
    // --------------------------------------------------------

    const int coarse_delay =
        delay0_samples -
        COARSE_RADIUS +
        best_position * COARSE_STEP;

    // --------------------------------------------------------
    // Recompute correlation exactly at coarse optimum.
    //
    // This makes the fine search independent of accumulated
    // floating-point error in the 4-sample coarse updates.
    // --------------------------------------------------------

    initial_correlation(
        coarse_delay,
        ci,
        cq);

    // --------------------------------------------------------
    // FINE SEARCH
    //
    // Exact V7 one-sample update.
    // --------------------------------------------------------

    float fine_best_metric =
        metric(ci, cq);

    int fine_best_offset =
        0;

    int fine_start =
        coarse_delay;

    for (int d = 1;
         d <= FINE_RADIUS;
         ++d)
    {
        const int old_pos =
            fine_start + d - 1;

        const int new_pos =
            old_pos + SYMBOL_SAMPLES;

        const float x_old =
            samples[old_pos];

        const float x_new =
            samples[new_pos];

        // ----------------------------------------------------
        // C' = (C - x_old + x_new*E) * P
        //
        // P = exp(+jw)
        // ----------------------------------------------------

        float ti =
            ci - x_old;

        float tq =
            cq;

        ti +=
            x_new * ec;

        tq +=
            x_new * es;

        const float ni =
            ti * pc -
            tq * ps;

        const float nq =
            ti * ps +
            tq * pc;

        ci = ni;
        cq = nq;

        const float m =
            metric(ci, cq);

        if (m > fine_best_metric)
        {
            fine_best_metric = m;
            fine_best_offset = d;
        }
    }

    // --------------------------------------------------------
    // Search negative side
    //
    // Recompute at coarse position.
    // --------------------------------------------------------

    initial_correlation(
        coarse_delay,
        ci,
        cq);

    for (int d = 1;
         d <= FINE_RADIUS;
         ++d)
    {
        const int old_pos =
            coarse_delay - d;

        const int new_pos =
            old_pos + SYMBOL_SAMPLES;

        const float x_old =
            samples[old_pos];

        const float x_new =
            samples[new_pos];

        // Move window one sample backwards.
        //
        // Inverse operation:
        //
        // C(d-1) =
        //     (C(d) * exp(-jw))
        //     + x_old*exp(-jw)
        //     - x_new*E*exp(-jw)

        const float inv_pc =
            pc;

        const float inv_ps =
            -ps;

        float ti =
            ci * inv_pc -
            cq * inv_ps;

        float tq =
            ci * inv_ps +
            cq * inv_pc;

        ti +=
            x_old;

        // x_new contribution is removed after rotating
        // the correlation backwards.

        const float rx =
            x_new * ec;

        const float rq =
            x_new * es;

        ti -=
            rx * inv_pc -
            rq * inv_ps;

        tq -=
            rx * inv_ps +
            rq * inv_pc;

        ci = ti;
        cq = tq;

        const float m =
            metric(ci, cq);

        if (m > fine_best_metric)
        {
            fine_best_metric = m;
            fine_best_offset = -d;
        }
    }

    // --------------------------------------------------------
    // Final delay
    // --------------------------------------------------------

    const int final_delay_samples =
        coarse_delay +
        fine_best_offset;

    const float final_delay =
        (float)final_delay_samples /
        SAMPLE_RATE;

    return final_delay;
}

// ============================================================
// refine_ft8_delay_v8_1()
//
// V8.1 - CONSERVATIVE COSTAS-ONLY VERSION
//
// Based on V7 COSTAS.
//
// Optimization:
//   - optimized initial correlation using groups of 4 samples
//
// Deliberately NOT optimized:
//   - coarse sliding update remains exactly the V7
//     one-sample-at-a-time update
//   - fine sliding update remains exactly the V7 update
//
// Purpose:
//   Compare V7 vs V8.1 and isolate the gain from the
//   optimized initial correlation without changing the
//   sliding-update mathematics.
//
// ============================================================

float refine_ft8_delay_v8_1(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)cand_to_print;

    constexpr int SAMPLE_RATE    = 12000;
    constexpr int NTONES         = 79;
    constexpr int SYMBOL_SAMPLES = 1920;

    // --------------------------------------------------------
    // Coarse search
    // --------------------------------------------------------

    constexpr int COARSE_RADIUS = 600;
    constexpr int COARSE_STEP   = 4;

    constexpr int MAX_COARSE =
        2 * COARSE_RADIUS / COARSE_STEP + 1;

    static_assert(
        MAX_COARSE <= 304,
        "MAX_COARSE too small");

    // --------------------------------------------------------
    // Fine search
    // --------------------------------------------------------

    constexpr int FINE_RADIUS = 6;
    constexpr int MAX_FINE    = 16;

    // --------------------------------------------------------
    // COSTAS ONLY
    // --------------------------------------------------------

    auto is_costas = [](int k) -> bool
    {
        return
            (k >= 0  && k <= 6)  ||
            (k >= 36 && k <= 42) ||
            (k >= 72 && k <= 78);
    };

    // ========================================================
    // Search limits
    // ========================================================

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

    // ========================================================
    // Score arrays
    // ========================================================

    float score_total[MAX_COARSE] = {};
    int   used_count[MAX_COARSE]  = {};

    // ========================================================
    // INITIAL CORRELATION V8.1
    //
    // Same mathematical correlation as V7:
    //
    //      C = sum x[n] * exp(-j*w*n)
    //
    // but the oscillator is evaluated using groups of four
    // samples.
    //
    // IMPORTANT:
    //
    // This does NOT change the sliding update.
    // ========================================================

    auto initial_correlation =
        [&](int start,
            int count,
            float w,
            float& out_i,
            float& out_q)
    {
        out_i = 0.0f;
        out_q = 0.0f;

        if (count <= 0)
            return;

        // ----------------------------------------------------
        // exp(-j*w)
        // ----------------------------------------------------

        const float cw = cosf(w);
        const float sw = sinf(w);

        // ----------------------------------------------------
        // Local coefficients:
        //
        // q0 = exp(-j*0*w)
        // q1 = exp(-j*1*w)
        // q2 = exp(-j*2*w)
        // q3 = exp(-j*3*w)
        // ----------------------------------------------------

        const float c0 = 1.0f;
        const float s0 = 0.0f;

        const float c1 =
            c0 * cw + s0 * sw;

        const float s1 =
            s0 * cw - c0 * sw;

        const float c2 =
            c1 * cw + s1 * sw;

        const float s2 =
            s1 * cw - c1 * sw;

        const float c3 =
            c2 * cw + s2 * sw;

        const float s3 =
            s2 * cw - c2 * sw;

        // ----------------------------------------------------
        // exp(-j*4*w)
        // ----------------------------------------------------

        const float p4c =
            c3 * cw + s3 * sw;

        const float p4s =
            s3 * cw - c3 * sw;

        // ----------------------------------------------------
        // Phase of current group.
        //
        // At group 0:
        //
        //      exp(-j*0*w) = 1
        //
        // Each group advances by 4 samples.
        // ----------------------------------------------------

        float gc = 1.0f;
        float gs = 0.0f;

        int pos = start;
        int remaining = count;

        // ----------------------------------------------------
        // Groups of four
        // ----------------------------------------------------

        while (remaining >= 4)
        {
            const float x0 =
                samples[pos + 0];

            const float x1 =
                samples[pos + 1];

            const float x2 =
                samples[pos + 2];

            const float x3 =
                samples[pos + 3];

            // ------------------------------------------------
            // Local correlation:
            //
            // L = x0*q0 + x1*q1 + x2*q2 + x3*q3
            //
            // q0 is real.
            // ------------------------------------------------

            const float li =
                x0 +
                x1 * c1 +
                x2 * c2 +
                x3 * c3;

            const float lq =
                x1 * s1 +
                x2 * s2 +
                x3 * s3;

            // ------------------------------------------------
            // Rotate local result by group phase.
            //
            // (li + j*lq) * (gc + j*gs)
            // ------------------------------------------------

            out_i +=
                li * gc -
                lq * gs;

            out_q +=
                li * gs +
                lq * gc;

            // ------------------------------------------------
            // Next group
            // ------------------------------------------------

            const float ngc =
                gc * p4c -
                gs * p4s;

            const float ngs =
                gc * p4s +
                gs * p4c;

            gc = ngc;
            gs = ngs;

            pos += 4;
            remaining -= 4;
        }

        // ----------------------------------------------------
        // Remaining 1..3 samples
        // ----------------------------------------------------

        if (remaining > 0)
        {
            float li = 0.0f;
            float lq = 0.0f;

            if (remaining >= 1)
            {
                const float x0 =
                    samples[pos + 0];

                li += x0 * c0;
            }

            if (remaining >= 2)
            {
                const float x1 =
                    samples[pos + 1];

                li += x1 * c1;
                lq += x1 * s1;
            }

            if (remaining >= 3)
            {
                const float x2 =
                    samples[pos + 2];

                li += x2 * c2;
                lq += x2 * s2;
            }

            out_i +=
                li * gc -
                lq * gs;

            out_q +=
                li * gs +
                lq * gc;
        }
    };

    // ========================================================
    // COARSE SEARCH
    //
    // IMPORTANT:
    //
    // The update below is deliberately copied from V7.
    //
    // One sample at a time.
    //
    // No V8 four-sample update here.
    // ========================================================

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        if (!is_costas(k))
            continue;

        const float f =
            freq +
            TONE_SPACING * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        // ----------------------------------------------------
        // V7 oscillator coefficient
        //
        // exp(+j*w)
        // ----------------------------------------------------

        const float pc =
            cosf(w);

        const float ps =
            sinf(w);

        // ----------------------------------------------------
        // End rotation:
        //
        // exp(-j*w*SYMBOL_SAMPLES)
        // ----------------------------------------------------

        const float end_angle =
            -w *
            (float)SYMBOL_SAMPLES;

        const float ec =
            cosf(end_angle);

        const float es =
            sinf(end_angle);

        // ----------------------------------------------------
        // Initial position
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
            w,
            ci,
            cq);

        // ----------------------------------------------------
        // Partial final symbol
        //
        // Same handling as V7.
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

        int pos = first_start;

        // ----------------------------------------------------
        // Coarse positions
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
            // EXACT V7 UPDATE
            //
            // Four iterations are NOT collapsed.
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

                // ------------------------------------------------
                // EXACT V7 complex rotation
                // ------------------------------------------------

                const float new_ci =
                    ti * pc -
                    tq * ps;

                const float new_cq =
                    ti * ps +
                    tq * pc;

                ci = new_ci;
                cq = new_cq;
            }

            pos += COARSE_STEP;
        }
    }

    // ========================================================
    // FIND COARSE MAXIMUM
    // ========================================================

    float best_coarse_score = -1.0f;
    int best_coarse_index = 0;

    for (int d = 0;
         d < ncoarse;
         ++d)
    {
        if (used_count[d] <= 0)
            continue;

        const float score =
            score_total[d] /
            (float)used_count[d];

        if (score > best_coarse_score)
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
    // FINE SEARCH LIMITS
    // ========================================================

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

    float accumulated_fine[MAX_FINE] = {};

    // ========================================================
    // FINE SEARCH
    //
    // EXACT V7 ONE-SAMPLE UPDATE
    // ========================================================

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        if (!is_costas(k))
            continue;

        const float f =
            freq +
            TONE_SPACING * (float)tones[k];

        const float w =
            TWO_PI * f /
            (float)SAMPLE_RATE;

        const float pc =
            cosf(w);

        const float ps =
            sinf(w);

        const float end_angle =
            -w *
            (float)SYMBOL_SAMPLES;

        const float ec =
            cosf(end_angle);

        const float es =
            sinf(end_angle);

        // ----------------------------------------------------
        // Initial position
        // ----------------------------------------------------

        const int pos =
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
            count,
            w,
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

            // ------------------------------------------------
            // EXACT V7 UPDATE
            // ------------------------------------------------

            const float ti =
                ci -
                x_old +
                x_new * ec;

            const float tq =
                cq +
                x_new * es;

            const float new_ci =
                ti * pc -
                tq * ps;

            const float new_cq =
                ti * ps +
                tq * pc;

            ci = new_ci;
            cq = new_cq;
        }
    }

    // ========================================================
    // FIND FINE MAXIMUM
    // ========================================================

    float best_fine_score = -1.0f;
    int best_fine_delay = coarse_delay;

    for (int d = 0;
         d < nfine;
         ++d)
    {
        if (accumulated_fine[d] >
            best_fine_score)
        {
            best_fine_score =
                accumulated_fine[d];

            best_fine_delay =
                fine_first + d;
        }
    }

    // ========================================================
    // RESULT
    // ========================================================

    return
        best_fine_delay /
        (float)SAMPLE_RATE;
}

// ============================================================
// V8 FT8 DELAY REFINEMENT
//
// Optimizations:
//
//   - COSTAS symbols only: 21 / 79
//   - Initial correlation processed in groups of 4 samples
//   - Coarse sliding correlation advanced directly by 4 samples
//   - Fine search remains sample-by-sample
//
// The scoring philosophy and delay search range remain unchanged.
//
// COSTAS positions:
//   0..6
//   36..42
//   72..78
//
// ============================================================

float refine_ft8_delay_v8(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    int cand_to_print)
{
    (void)cand_to_print;

    constexpr int SAMPLE_RATE   = 12000;
    constexpr int NTONES        = 79;
    constexpr int SYMBOL_SAMPLES = 1920;

    constexpr int COARSE_RADIUS = 600;   // +/- 50 ms
    constexpr int COARSE_STEP   = 4;
    constexpr int MAX_COARSE    = 301;

    constexpr int FINE_RADIUS   = 6;
    constexpr int MAX_FINE      = 13;

    // --------------------------------------------------------
    // Only the 21 Costas symbols
    // --------------------------------------------------------

    auto is_costas = [](int k) -> bool
    {
        return
            (k >= 0  && k <= 6)  
            // ||
            // (k >= 36 && k <= 42) ||
            // (k >= 72 && k <= 78)
            ;
    };

    // --------------------------------------------------------
    // Search limits
    // --------------------------------------------------------

    const int center_delay =
        (int)lroundf(delay0 * SAMPLE_RATE);

    const int first_delay =
        std::max(0, center_delay - COARSE_RADIUS);

    const int last_delay =
        std::min(num_samples - 1,
                 center_delay + COARSE_RADIUS);

    // --------------------------------------------------------
    // Coarse accumulation
    // --------------------------------------------------------

    float score_total[MAX_COARSE] = {};
    int   used_count[MAX_COARSE]  = {};

    const int coarse_count =
        ((last_delay - first_delay) / COARSE_STEP) + 1;

    // --------------------------------------------------------
    // Helper:
    //
    // Initial correlation of one symbol.
    //
    // Instead of rotating the oscillator for every sample,
    // calculate four local coefficients once and process
    // samples in groups of four.
    //
    // This substantially reduces the number of complex
    // oscillator rotations.
    // --------------------------------------------------------

    auto initial_correlation_v8 =
        [&](int start,
            int count,
            float w,
            float& out_i,
            float& out_q)
    {
        out_i = 0.0f;
        out_q = 0.0f;

        if (count <= 0)
            return;

        // ----------------------------------------------------
        // Four coefficients:
        //
        // q[r] = exp(-j*w*r)
        //
        // r = 0,1,2,3
        // ----------------------------------------------------

        const float cw = cosf(w);
        const float sw = sinf(w);

        float c0 = 1.0f;
        float s0 = 0.0f;

        float c1 = c0 * cw + s0 * sw;
        float s1 = s0 * cw - c0 * sw;

        float c2 = c1 * cw + s1 * sw;
        float s2 = s1 * cw - c1 * sw;

        float c3 = c2 * cw + s2 * sw;
        float s3 = s2 * cw - c2 * sw;

        // ----------------------------------------------------
        // Phase corresponding to the beginning of each
        // 4-sample group.
        //
        // group_phase = exp(-j*w*(4*g))
        // ----------------------------------------------------

        float gc = 1.0f;
        float gs = 0.0f;

        // exp(-j*4*w)
        const float p4c =
            c3 * cw + s3 * sw;

        const float p4s =
            s3 * cw - c3 * sw;

        int pos = start;
        int remaining = count;

        while (remaining >= 4)
        {
            const float x0 = samples[pos + 0];
            const float x1 = samples[pos + 1];
            const float x2 = samples[pos + 2];
            const float x3 = samples[pos + 3];

            // ------------------------------------------------
            // Local correlation:
            //
            // L = x0*q0 + x1*q1 + x2*q2 + x3*q3
            // ------------------------------------------------

            const float li =
                x0 +
                x1 * c1 +
                x2 * c2 +
                x3 * c3;

            const float lq =
                x1 * s1 +
                x2 * s2 +
                x3 * s3;

            // ------------------------------------------------
            // Apply group phase:
            //
            // L * exp(-j*w*(4*g))
            //
            // (li+j*lq)*(gc+j*gs)
            // ------------------------------------------------

            out_i += li * gc - lq * gs;
            out_q += li * gs + lq * gc;

            // ------------------------------------------------
            // Advance group phase by 4 samples
            // ------------------------------------------------

            const float ngc =
                gc * p4c - gs * p4s;

            const float ngs =
                gc * p4s + gs * p4c;

            gc = ngc;
            gs = ngs;

            pos += 4;
            remaining -= 4;
        }

        // ----------------------------------------------------
        // Remaining 1..3 samples
        // ----------------------------------------------------

        if (remaining > 0)
        {
            float lc = 0.0f;
            float ls = 0.0f;

            if (remaining >= 1)
            {
                const float x = samples[pos + 0];
                lc += x * c0;
                ls += x * s0;
            }

            if (remaining >= 2)
            {
                const float x = samples[pos + 1];
                lc += x * c1;
                ls += x * s1;
            }

            if (remaining >= 3)
            {
                const float x = samples[pos + 2];
                lc += x * c2;
                ls += x * s2;
            }

            out_i += lc * gc - ls * gs;
            out_q += lc * gs + ls * gc;
        }
    };

    // --------------------------------------------------------
    // COARSE SEARCH
    // --------------------------------------------------------

    for (int k = 0; k < NTONES; ++k)
    {
        if (!is_costas(k))
            continue;

        const float f =
            freq + TONE_SPACING * tones[k];

        const float w =
            TWO_PI * f / SAMPLE_RATE;

        const float pc = cosf(w);
        const float ps = -sinf(w);

        // Same endpoint rotation convention as V7
        const float end_angle =
            -w * SYMBOL_SAMPLES;

        const float ec = cosf(end_angle);
        const float es = sinf(end_angle);

        // ----------------------------------------------------
        // Initial correlation at first coarse delay
        // ----------------------------------------------------

        int pos = first_delay + k * SYMBOL_SAMPLES;

        if (pos < 0 || pos >= num_samples)
            continue;

        int count =
            std::min(SYMBOL_SAMPLES,
                     num_samples - pos);

        if (count <= 0)
            continue;

        float ci, cq;

        initial_correlation_v8(
            pos,
            count,
            w,
            ci,
            cq);

        // ----------------------------------------------------
        // Score at current delay
        // ----------------------------------------------------

        const float mag2 =
            ci * ci + cq * cq;

        score_total[0] += mag2;
        used_count[0]++;

        // ----------------------------------------------------
        // Coarse sliding search
        //
        // Advance directly by 4 samples.
        //
        // This is algebraically equivalent to four V7
        // one-sample sliding updates, but avoids performing
        // four complete complex rotations.
        // ----------------------------------------------------

        int d = 0;

        while (d + COARSE_STEP < coarse_count)
        {
            // Current absolute position
            const int old_pos =
                first_delay + d * COARSE_STEP
                + k * SYMBOL_SAMPLES;

            // Need four new samples
            if (old_pos + SYMBOL_SAMPLES + 3 >= num_samples)
                break;

            // ------------------------------------------------
            // Four consecutive V7 updates collapsed into one.
            //
            // C1 = (C0 - old0 + new0*E) * P
            // C2 = (C1 - old1 + new1*E) * P
            // C3 = (C2 - old2 + new2*E) * P
            // C4 = (C3 - old3 + new3*E) * P
            //
            // Precompute P^1..P^4.
            // ------------------------------------------------

            const float p1c = pc;
            const float p1s = ps;

            const float p2c =
                p1c * p1c - p1s * p1s;

            const float p2s =
                2.0f * p1c * p1s;

            const float p3c =
                p2c * p1c - p2s * p1s;

            const float p3s =
                p2c * p1s + p2s * p1c;

            const float p4c =
                p3c * p1c - p3s * p1s;

            const float p4s =
                p3c * p1s + p3s * p1c;

            // ------------------------------------------------
            // New-sample complex terms:
            //
            // (new * E) * P^n
            // ------------------------------------------------

            const float q0c =
                ec * p4c - es * p4s;

            const float q0s =
                ec * p4s + es * p4c;

            const float q1c =
                ec * p3c - es * p3s;

            const float q1s =
                ec * p3s + es * p3c;

            const float q2c =
                ec * p2c - es * p2s;

            const float q2s =
                ec * p2s + es * p2c;

            const float q3c =
                ec * p1c - es * p1s;

            const float q3s =
                ec * p1s + es * p1c;

            // ------------------------------------------------
            // Old samples:
            //
            // -old * P^n
            // ------------------------------------------------

            const float o0 = samples[old_pos + 0];
            const float o1 = samples[old_pos + 1];
            const float o2 = samples[old_pos + 2];
            const float o3 = samples[old_pos + 3];

            const float n0 =
                samples[old_pos + SYMBOL_SAMPLES + 0];

            const float n1 =
                samples[old_pos + SYMBOL_SAMPLES + 1];

            const float n2 =
                samples[old_pos + SYMBOL_SAMPLES + 2];

            const float n3 =
                samples[old_pos + SYMBOL_SAMPLES + 3];

            float ti =
                ci * p4c - cq * p4s;

            float tq =
                ci * p4s + cq * p4c;

            ti +=
                -o0 * p4c +
                 n0 * q0c;

            tq +=
                -o0 * p4s +
                 n0 * q0s;

            ti +=
                -o1 * p3c +
                 n1 * q1c;

            tq +=
                -o1 * p3s +
                 n1 * q1s;

            ti +=
                -o2 * p2c +
                 n2 * q2c;

            tq +=
                -o2 * p2s +
                 n2 * q2s;

            ti +=
                -o3 * p1c +
                 n3 * q3c;

            tq +=
                -o3 * p1s +
                 n3 * q3s;

            ci = ti;
            cq = tq;

            d += COARSE_STEP;

            const int idx =
                d / COARSE_STEP;

            const float m2 =
                ci * ci + cq * cq;

            score_total[idx] += m2;
            used_count[idx]++;
        }
    }

    // --------------------------------------------------------
    // Find coarse maximum
    // --------------------------------------------------------

    float best_score = -1.0f;
    int best_index = 0;

    for (int i = 0; i < coarse_count; ++i)
    {
        if (used_count[i] <= 0)
            continue;

        const float score =
            score_total[i] /
            (float)used_count[i];

        if (score > best_score)
        {
            best_score = score;
            best_index = i;
        }
    }

    const int coarse_delay =
        first_delay +
        best_index * COARSE_STEP;

    // --------------------------------------------------------
    // FINE SEARCH
    //
    // Keep the original V7 one-sample sliding correlation.
    // This is the precision-critical part.
    // --------------------------------------------------------

    const int fine_first =
        std::max(first_delay,
                 coarse_delay - FINE_RADIUS);

    const int fine_last =
        std::min(last_delay,
                 coarse_delay + FINE_RADIUS);

    const int fine_count =
        fine_last - fine_first + 1;

    float accumulated_fine[MAX_FINE] = {};
    int   fine_used[MAX_FINE]        = {};

    for (int k = 0; k < NTONES; ++k)
    {
        if (!is_costas(k))
            continue;

        const float f =
            freq + TONE_SPACING * tones[k];

        const float w =
            TWO_PI * f / SAMPLE_RATE;

        const float pc = cosf(w);
        const float ps = -sinf(w);

        const float end_angle =
            -w * SYMBOL_SAMPLES;

        const float ec = cosf(end_angle);
        const float es = sinf(end_angle);

        int pos =
            fine_first +
            k * SYMBOL_SAMPLES;

        if (pos < 0 || pos >= num_samples)
            continue;

        // ----------------------------------------------------
        // IMPORTANT:
        //
        // Use the actual available count, not
        // SYMBOL_SAMPLES unconditionally.
        // ----------------------------------------------------

        const int count =
            std::min(SYMBOL_SAMPLES,
                     num_samples - pos);

        if (count <= 0)
            continue;

        float ci, cq;

        initial_correlation_v8(
            pos,
            count,
            w,
            ci,
            cq);

        for (int d = 0; d < fine_count; ++d)
        {
            if (d > 0)
            {
                const int old_pos =
                    fine_first +
                    k * SYMBOL_SAMPLES +
                    d - 1;

                const int new_sample =
                    old_pos +
                    SYMBOL_SAMPLES;

                if (new_sample >= num_samples)
                    break;

                const float x_old =
                    samples[old_pos];

                const float x_new =
                    samples[new_sample];

                // Same one-sample sliding update as V7
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

            const float mag2 =
                ci * ci + cq * cq;

            accumulated_fine[d] += mag2;
            fine_used[d]++;
        }
    }

    // --------------------------------------------------------
    // Find fine maximum
    // --------------------------------------------------------

    float best_fine_score = -1.0f;
    int best_fine_index = 0;

    for (int d = 0; d < fine_count; ++d)
    {
        if (fine_used[d] <= 0)
            continue;

        const float score =
            accumulated_fine[d] /
            (float)fine_used[d];

        if (score > best_fine_score)
        {
            best_fine_score = score;
            best_fine_index = d;
        }
    }

    const int final_delay_samples =
        fine_first + best_fine_index;

    return
        (float)final_delay_samples /
        (float)SAMPLE_RATE;
}




// ============================================================
// refine_ft8_delay_v7()
//
// V7 - SELECTIVE SYMBOL VERSION
//
// Modes:
//
//   FULL
//       all 79 symbols
//
//   COSTAS
//       21 Costas symbols
//
//   COSTAS_DATA
//       21 Costas + 12 selected data symbols = 33
//
// IMPORTANT:
//
//   The correlation algorithm itself is unchanged.
//
//   Only the symbols participating in the accumulated score
//   are changed.
//
// Coarse:
//   ±50 ms
//   4-sample resolution
//
// Fine:
//   ±6 samples around coarse result
//   1-sample resolution
//
// ============================================================

float refine_ft8_delay_v7(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    float delay0,
    float freq,
    SymbolMode mode,
    int cand_to_print)
{
    constexpr int SAMPLE_RATE = 12000;
    constexpr int NTONES      = 79;
    constexpr int SYMBOL_SAMPLES = 1920;

    // --------------------------------------------------------
    // Coarse search
    // --------------------------------------------------------

    constexpr int COARSE_RADIUS = 600;
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
    // IMPORTANT:
    //
    // This is deliberately the same mathematical operation
    // as before.
    //
    // No LUT optimization is introduced here.
    // ========================================================

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

        float c = 1.0f;
        float s = 0.0f;

        int n = 0;

        // ----------------------------------------------------
        // Four samples per iteration.
        // ----------------------------------------------------

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

        // ----------------------------------------------------
        // Remaining samples.
        // ----------------------------------------------------

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
        // ----------------------------------------------------
        // SELECTIVE SYMBOL PROCESSING
        // ----------------------------------------------------

        if (!use_symbol(k, mode))
            continue;

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

        if (score > best_coarse_score)
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
    // Fine search
    // ========================================================

    for (int k = 0;
         k < NTONES;
         ++k)
    {
        // ----------------------------------------------------
        // SELECTIVE SYMBOL PROCESSING
        // ----------------------------------------------------

        if (!use_symbol(k, mode))
            continue;

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

        // ----------------------------------------------------
        // IMPORTANT:
        //
        // pos is the beginning of the symbol for the FIRST
        // fine delay.
        //
        // Each d below corresponds to:
        //
        //     fine_first + d
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
            num_samples - pos;

        const int count =
            std::min(
                SYMBOL_SAMPLES,
                available);

        if (count <= 0)
            continue;

        float ci;
        float cq;

        // ----------------------------------------------------
        // BUG FIX:
        //
        // Previously this passed SYMBOL_SAMPLES instead of
        // count. That could overread a partial final symbol.
        // ----------------------------------------------------

        initial_correlation(
            pos,
            count,
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

    return
        fine_best_delay /
        (float)SAMPLE_RATE;
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
// WAV LOADER - ESP32 / LittleFS
// ============================================================

static float* load_wav(
    const char* path,
    int* out_num_samples,
    int* out_num_channels,
    int* out_sample_rate)
{
    File f = LittleFS.open(path, "r");

    if (!f)
    {
        Serial.printf("ERROR: cannot open %s\n", path);
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
        Serial.printf("ERROR: unsupported WAV format %u\n",
                      audio_format);
        f.close();
        return nullptr;
    }

    if (bits_per_sample != 16)
    {
        Serial.printf("ERROR: only 16-bit PCM supported, got %u\n",
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

    const int bytes_per_sample = bits_per_sample / 8;
    const int frame_bytes = bytes_per_sample * num_channels;
    const int num_samples = data_size / frame_bytes;

    float* samples =
        (float*)malloc((size_t)num_samples * sizeof(float));

    if (!samples)
    {
        Serial.printf("ERROR: malloc failed for %d samples\n",
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
            Serial.println("ERROR: unexpected EOF");
            free(samples);
            f.close();
            return nullptr;
        }

        const int16_t v =
            (int16_t)(b[0] | ((uint16_t)b[1] << 8));

        samples[n] = (float)v / 32768.0f;

        for (int ch = 1; ch < num_channels; ++ch)
            f.seek(f.position() + bytes_per_sample);
    }

    f.close();

    *out_num_samples  = num_samples;
    *out_num_channels = num_channels;
    *out_sample_rate  = (int)sample_rate;

    return samples;
}

#else

// ============================================================
// WAV LOADER - Linux / Ubuntu
// ============================================================

static uint16_t rd_u16(const uint8_t* p)
{
    return (uint16_t)p[0] |
           ((uint16_t)p[1] << 8);
}

static uint32_t rd_u32(const uint8_t* p)
{
    return (uint32_t)p[0] |
           ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) |
           ((uint32_t)p[3] << 24);
}

static float* load_wav(
    const char* path,
    int* out_num_samples,
    int* out_num_channels,
    int* out_sample_rate)
{
    std::ifstream f(path, std::ios::binary);

    if (!f && path[0] == '/')
    {
        f.clear();
        f.open(path + 1, std::ios::binary);
    }

    if (!f)
    {
        std::printf("ERROR: cannot open %s\n", path);
        return nullptr;
    }

    uint8_t riff[12];

    if (!f.read((char*)riff, sizeof(riff)))
    {
        std::printf("ERROR: WAV header too short\n");
        return nullptr;
    }

    if (std::memcmp(riff, "RIFF", 4) != 0 ||
        std::memcmp(riff + 8, "WAVE", 4) != 0)
    {
        std::printf("ERROR: not a RIFF/WAVE file\n");
        return nullptr;
    }

    uint16_t audio_format = 0;
    uint16_t num_channels = 0;
    uint32_t sample_rate = 0;
    uint16_t bits_per_sample = 0;
    uint16_t block_align = 0;

    std::streamoff data_pos = 0;
    uint32_t data_size = 0;

    bool have_fmt = false;
    bool have_data = false;

    while (f)
    {
        uint8_t chunk_header[8];

        if (!f.read((char*)chunk_header, sizeof(chunk_header)))
            break;

        const uint32_t chunk_size = rd_u32(chunk_header + 4);
        const std::streamoff chunk_start = f.tellg();

        if (std::memcmp(chunk_header, "fmt ", 4) == 0)
        {
            if (chunk_size < 16)
            {
                std::printf("ERROR: invalid fmt chunk\n");
                return nullptr;
            }

            std::vector<uint8_t> fmt(chunk_size);

            if (!f.read((char*)fmt.data(), chunk_size))
            {
                std::printf("ERROR: cannot read fmt chunk\n");
                return nullptr;
            }

            audio_format = rd_u16(fmt.data());
            num_channels = rd_u16(fmt.data() + 2);
            sample_rate = rd_u32(fmt.data() + 4);
            block_align = rd_u16(fmt.data() + 12);
            bits_per_sample = rd_u16(fmt.data() + 14);

            have_fmt = true;
        }
        else if (std::memcmp(chunk_header, "data", 4) == 0)
        {
            data_pos = chunk_start;
            data_size = chunk_size;
            have_data = true;

            if (have_fmt)
                break;
        }
        else
        {
            f.seekg((std::streamoff)chunk_size, std::ios::cur);
        }

        if (chunk_size & 1)
            f.seekg(1, std::ios::cur);

        if (have_fmt && have_data)
            break;
    }

    if (!have_fmt || !have_data)
    {
        std::printf("ERROR: missing fmt/data chunk\n");
        return nullptr;
    }

    if (audio_format != 1)
    {
        std::printf("ERROR: unsupported WAV format %u\n",
                    (unsigned)audio_format);
        return nullptr;
    }

    if (bits_per_sample != 16)
    {
        std::printf("ERROR: only 16-bit PCM supported, got %u\n",
                    (unsigned)bits_per_sample);
        return nullptr;
    }

    if (num_channels < 1)
    {
        std::printf("ERROR: invalid channel count\n");
        return nullptr;
    }

    if (block_align == 0)
        block_align = (uint16_t)(num_channels * 2);

    const uint32_t num_samples_u32 =
        data_size / block_align;

    if (num_samples_u32 == 0)
    {
        std::printf("ERROR: WAV contains no samples\n");
        return nullptr;
    }

    const int num_samples = (int)num_samples_u32;

    float* samples =
        (float*)std::malloc((size_t)num_samples * sizeof(float));

    if (!samples)
    {
        std::printf("ERROR: malloc failed for %d samples\n",
                    num_samples);
        return nullptr;
    }

    f.clear();
    f.seekg(data_pos);

    for (int n = 0; n < num_samples; ++n)
    {
        uint8_t b[2];

        if (!f.read((char*)b, 2))
        {
            std::printf("ERROR: unexpected EOF\n");
            std::free(samples);
            return nullptr;
        }

        const int16_t v =
            (int16_t)(b[0] | ((uint16_t)b[1] << 8));

        samples[n] = (float)v / 32768.0f;

        const std::streamoff skip =
            (std::streamoff)block_align - 2;

        if (skip > 0)
            f.seekg(skip, std::ios::cur);
    }

    *out_num_samples  = num_samples;
    *out_num_channels = (int)num_channels;
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
    const uint32_t heap_before = ESP.getFreeHeap();
#else
    const uint32_t heap_before = 0;
#endif

    float* samples =
        load_wav(
            filename,
            &num_samples,
            &num_channels,
            &sample_rate);

#if defined(ARDUINO)
    const uint32_t heap_after_load = ESP.getFreeHeap();
#else
    const uint32_t heap_after_load = 0;
#endif

    if (!samples)
    {
        printf("LOAD FAILED\n");
        return;
    }

    printf("Samples          : %d\n", num_samples);
    printf("Duration         : %.3f s\n",
           (float)num_samples / sample_rate);
    printf("Channels         : %d\n", num_channels);
    printf("Sample rate      : %d Hz\n", sample_rate);
    printf("Heap before load : %u\n", heap_before);
    printf("Heap after load  : %u\n", heap_after_load);

    if (sample_rate != SAMPLE_RATE ||
        num_channels != 1)
    {
        printf("ERROR: WAV format mismatch\n");
        free(samples);
        return;
    }

    // --------------------------------------------------------
    // Extract 79 tones ONCE.
    // Outside benchmark.
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

    // ========================================================
    // V8.1 BENCHMARK
    // ========================================================

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("V9 : COSTAS ONLY          symbols: 21\n");
    printf("------------------------------------------------------------\n");

#if defined(ARDUINO)

    const uint32_t t0 = micros();

#else

    const auto t0 =
        std::chrono::steady_clock::now();

#endif

    const float estimated_delay =
        refine_ft8_delay_v10g(
            samples,
            num_samples,
            tones,
            true_delay,
            TEST_FREQ,
            1);

#if defined(ARDUINO)

    const uint32_t elapsed_us =
        micros() - t0;

    const double elapsed_ms =
        (double)elapsed_us / 1000.0;

#else

    const auto elapsed_us =
        std::chrono::duration_cast<
            std::chrono::microseconds>(
                std::chrono::steady_clock::now() -
                t0).count();

    const double elapsed_ms =
        (double)elapsed_us / 1000.0;

#endif

    const float error_s =
        estimated_delay - true_delay;

    const float error_ms =
        error_s * 1000.0f;

    const float error_samples =
        error_s * SAMPLE_RATE;

    printf(
        "Estimated delay  : %10.6f s\n",
        estimated_delay);

    printf(
        "Error            : %+10.3f ms\n",
        error_ms);

    printf(
        "Error            : %+10.3f samples\n",
        error_samples);

    printf(
        "Execution time   : %10.3f ms\n",
        elapsed_ms);

#if defined(ESP32)

    printf(
        "Free heap        : %u bytes\n",
        ESP.getFreeHeap());

#endif

    // --------------------------------------------------------
    // Reference
    // --------------------------------------------------------

    printf("\n");
    printf("============================================================\n");
    printf("V9 SUMMARY: %s\n", filename);
    printf("============================================================\n");

    printf(
        "True delay       : %.6f s\n",
        true_delay);

    printf(
        "Estimated delay  : %.6f s\n",
        estimated_delay);

    printf(
        "Error            : %+0.3f ms\n",
        error_ms);

    printf(
        "Error            : %+0.3f samples\n",
        error_samples);

    printf(
        "Execution time   : %.3f ms\n",
        elapsed_ms);

    printf(
        "Symbols used     : 21 / 79\n");

    printf(
        "Symbol ratio     : %.3f\n",
        21.0 / 79.0);

    printf(
        "Ideal speedup    : %.2fx\n",
        79.0 / 21.0);

    printf(
        "============================================================\n");

    free(samples);

#if defined(ESP32)

    printf(
        "After free() heap: %u bytes\n",
        ESP.getFreeHeap());

#endif

    printf("WAV released.\n");
}

// ============================================================
// BENCHMARK DRIVER
// ============================================================

static void run_all_tests()
{
    printf("\r\n\r\n");

    printf(
        "############################################################\n");

    printf(
        " FT8 DELAY REFINEMENT - SYMBOL A/B/C BENCHMARK\n");

    printf(
        "############################################################\n");

#if defined(ESP32)

    printf(
        "Platform         : ESP32-S3 / Arduino\n");

    printf(
        "CPU frequency    : %u MHz\n",
        ESP.getCpuFreqMHz());

#else

    printf(
        "Platform         : Linux / Ubuntu\n");

#endif

    printf(
        "Sample rate      : %d Hz\n",
        SAMPLE_RATE);

    printf(
        "Symbol duration  : %d samples = %.3f ms\n",
        SPS,
        1000.0f * SPS / SAMPLE_RATE);

    printf(
        "Test frequency   : %.3f Hz\n",
        TEST_FREQ);

    printf(
        "FULL             : %d symbols\n",
        count_used_symbols(SymbolMode::FULL));

    printf(
        "COSTAS           : %d symbols\n",
        count_used_symbols(SymbolMode::COSTAS));

    printf(
        "COSTAS+DATA      : %d symbols\n",
        count_used_symbols(SymbolMode::COSTAS_DATA));

#if defined(ESP32)

    printf(
        "Free heap        : %u bytes\n",
        ESP.getFreeHeap());

    printf(
        "Free PSRAM       : %u bytes\n",
        ESP.getFreePsram());

    printf("\r\n");
    printf("Mounting LittleFS...\n");

    if (!LittleFS.begin(true))
    {
        printf(
            "ERROR: LittleFS.begin() failed\n");

        return;
    }

    printf(
        "LittleFS mounted.\n");

#endif

    run_one_test(
        "/test_real_ft8_D1500.wav",
        1.500f);

    run_one_test(
        "/test_real_ft8_D2000.wav",
        2.000f);

    run_one_test(
        "/test_real_ft8_D2500.wav",
        2.500f);




    printf("\r\n\r\n");

    printf(
        "############################################################\n");

    printf(
        "A/B/C BENCHMARK COMPLETE\n");

#if defined(ESP32)

    printf(
        "Final free heap : %u bytes\n",
        ESP.getFreeHeap());

    printf(
        "Final free PSRAM: %u bytes\n",
        ESP.getFreePsram());

#endif

    printf(
        "############################################################\n");
}

// ============================================================
// PLATFORM ENTRY POINT
// ============================================================

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
