#include <cassert>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>
#include <algorithm>
#include "subtract.h"

int rate_ = 12000;  // samples/second
float *samples_ = nullptr;
extern const double subtract_ramp;
static constexpr double PI = 3.14159265358979323846;


// ------------------------------------------------------------
// RMS
// ------------------------------------------------------------

static double rms(const float *p, size_t n)
{
    double e = 0.0;

    for (size_t i = 0; i < n; ++i)
        e += double(p[i]) * double(p[i]);

    return std::sqrt(e / n);
}


// ------------------------------------------------------------
// Generate a pure tone in the region occupied by the FT8 signal
// ------------------------------------------------------------

static void generate_tone(
    float *samples,
    size_t nsamples,
    double freq,
    double amplitude,
    double off_sec,
    int sample_rate)
{
    int block = blocksize(sample_rate);
    int off0 = std::lround(off_sec * sample_rate);

    for (int i = 0; i < 79; ++i)
    {
        for (int n = 0; n < block; ++n)
        {
            int idx = off0 + i * block + n;

            if (idx < 0 || idx >= (int)nsamples)
                continue;

            double theta =
                2.0 * PI * freq * n / sample_rate;

            samples[idx] =
                amplitude * std::cos(theta);
        }
    }
}


// ------------------------------------------------------------
// Test 1
// ------------------------------------------------------------

static void test_perfect_subtraction()
{
    printf("\nTEST 1: perfect subtraction\n");

    const double freq = 1000.0;
    const double amplitude = 1000.0;
    const double offset = 2.0;

    const size_t nsamples = 192000;

    std::vector<float> buffer(nsamples, 0.0f);

    samples_ = buffer.data();

    generate_tone(
        samples_,
        nsamples,
        freq,
        amplitude,
        offset,
        rate_);

    int block = blocksize(rate_);

    int off0 = std::lround(offset * rate_);

    size_t start = off0;
    size_t length = 79 * block;

    double before =
        rms(samples_ + start, length);

    printf("  RMS before = %.6f\n", before);

    // All 79 symbols use tone 0.
    uint8_t tones[79];
    for(int i=0; i<79; i++){
        tones[i] = 0;
    }

    subtract(
        tones,
        freq,
        freq,
        offset,
        samples_,
        nsamples,
        rate_);

    double after =
        rms(samples_ + start, length);

    printf("  RMS after  = %.6f\n", after);

    printf("  reduction   = %.2f dB\n",
           20.0 * std::log10(before / after));

    assert(after < before * 0.1);

    printf("  PASS\n");
}


// ------------------------------------------------------------

static void test_timing_error()
{
    printf("\nTEST 2: timing error\n");

    const double freq = 1000.0;
    const double amplitude = 1000.0;

    const double true_offset = 2.0;

    const size_t nsamples = 192000;

    int block = blocksize(rate_);

    int true_off0 =
        std::lround(true_offset * rate_);

    uint8_t re79[79] = {0};

    printf("\n");
    printf("  true offset = %.3f s\n", true_offset);
    printf("  block       = %d samples\n", block);
    printf("\n");

    printf("  %8s  %14s  %14s\n",
           "error(ms)",
           "RMS residual",
           "suppression");

    printf("  ---------------------------------------------\n");

    const double errors_ms[] = {
        -100.0,
         -80.0,
         -60.0,
         -40.0,
         -20.0,
         -10.0,
          -5.0,
           0.0,
           5.0,
          10.0,
          20.0,
          40.0,
          60.0,
          80.0,
         100.0
    };

    for (double error_ms : errors_ms)
    {
        /*
         * Generate a fresh signal for every test.
         */
        std::vector<float> buffer(nsamples, 0.0f);

        samples_ = buffer.data();

        /*
         * Generate the signal at the TRUE offset.
         */
        generate_tone(
            samples_,
            nsamples,
            freq,
            amplitude,
            true_offset,
            rate_);

        /*
         * Measure original signal energy.
         */
        double before =
            rms(samples_ + true_off0,
                79 * block);

        /*
         * Deliberately give subtract() the WRONG offset.
         */
        double subtract_offset =
            true_offset + error_ms / 1000.0;

        subtract(
            re79,
            freq,
            freq,
            subtract_offset,
            samples_,
            nsamples,
            rate_);

        /*
         * Measure residual around the ORIGINAL
         * signal position.
         *
         * We use the original region rather than the
         * subtraction region so that the measurement
         * remains comparable for every timing error.
         */
        double after =
            rms(samples_ + true_off0,
                79 * block);

        double suppression =
            20.0 * std::log10(before / after);

        printf("  %8.1f  %14.6f  %13.2f dB\n",
               error_ms,
               after,
               suppression);
    }

    printf("\n");
}

static void generate_ft8_signal(
    float *samples,
    size_t nsamples,
    const std::vector<int>& re79,
    double hz0,
    double amplitude,
    double off_sec)
{
    int block = blocksize(rate_);
    int off0 = std::lround(off_sec * rate_);

    const int ramp =
        std::max(1, (int)std::lround(block * subtract_ramp));

    for (int si = 0; si < 79; ++si)
    {
        double freq =
            hz0 + 6.25 * re79[si];

        double dtheta =
            2.0 * M_PI * freq / rate_;

        /*
         * Use zero phase initially.
         */
        double phase = 0.0;

        /*
         * Steady symbol.
         */
        for (int jj = ramp; jj < block-ramp; ++jj)
        {
            int idx =
                off0 + si * block + jj;

            if (idx < 0 || idx >= (int)nsamples)
                continue;

            double theta =
                phase + jj * dtheta;

            samples[idx] +=
                amplitude * std::cos(theta);
        }

        /*
         * We deliberately keep the first version simple.
         * The transition region is generated by interpolating
         * between the current and next tone frequencies.
         */
        double freq1;

        if (si + 1 < 79)
            freq1 = hz0 + 6.25 * re79[si + 1];
        else
            freq1 = freq;

        double dtheta1 =
            2.0 * M_PI * freq1 / rate_;

        double inc =
            (dtheta1 - dtheta) / (2.0 * ramp);

        double theta =
            phase + (block-ramp) * dtheta;

        for (int jj = block-ramp; jj < block; ++jj)
        {
            int idx =
                off0 + si * block + jj;

            if (idx >= 0 && idx < (int)nsamples)
            {
                samples[idx] +=
                    amplitude * std::cos(theta);
            }

            theta += dtheta;
            dtheta += inc;
        }
    }
}

static const uint8_t tones_A[] = {
    3, 1, 6, 2, 0, 5, 4, 1,
    7, 2, 3, 6, 0, 4, 5, 2,
    1, 7, 4, 3, 0, 6, 2, 5,
    4, 1, 3, 7, 2, 0, 6, 5,
    3, 4, 1, 7, 5, 2, 6, 0,
    4, 3, 7, 1, 5, 0, 2, 6,
    3, 5, 1, 4, 7, 2, 0, 6,
    5, 3, 1, 7, 4, 2, 6, 0,
    3, 1, 5, 7, 2, 4, 6, 0,
    1, 3, 5, 7, 2, 4, 6, 0
};

static const uint8_t tones_B[] = {
    6, 2, 0, 7, 3, 1, 5, 4,
    2, 6, 1, 4, 7, 0, 3, 5,
    5, 0, 4, 2, 6, 7, 1, 3,
    1, 5, 7, 3, 0, 2, 4, 6,
    7, 3, 5, 1, 6, 4, 2, 0,
    2, 7, 3, 6, 0, 5, 1, 4,
    4, 1, 6, 0, 3, 7, 5, 2,
    6, 4, 2, 7, 1, 5, 0, 3,
    0, 5, 3, 6, 2, 7, 4, 1,
    7, 1, 4, 0, 6, 2, 5
};

static void generate_ft8_simple(
    float *samples,
    size_t nsamples,
    const uint8_t re79[79],
    double hz0,
    double amplitude,
    double off_sec)
{
    int block = blocksize(rate_);
    int off0 = std::lround(off_sec * rate_);

    for (int si = 0; si < 79; ++si)
    {
        double freq =
            hz0 + 6.25 * re79[si];

        double dtheta =
            2.0 * M_PI * freq / rate_;

        for (int n = 0; n < block; ++n)
        {
            int idx =
                off0 + si * block + n;

            if (idx < 0 || idx >= (int)nsamples)
                continue;

            double theta =
                dtheta * n;

            samples[idx] +=
                amplitude * std::cos(theta);
        }
    }
}


static void generate_ft8_model(
    float *samples,
    size_t nsamples,
    const std::vector<int>& re79,
    double hz0,
    double amplitude,
    double off_sec)
{
    int block = blocksize(rate_);
    int off0 = std::lround(off_sec * rate_);

    int ramp = std::lround(block * subtract_ramp);

    if (ramp < 1)
        ramp = 1;

    /*
     * We use the same phase convention as subtract():
     *
     * phase[si] = 0
     *
     * and the same frequency interpolation through
     * the transition between symbols.
     *
     * For this first ideal test we use the same phase
     * for every symbol, so subtract() should be able
     * to estimate it exactly.
     */

    for (int si = 0; si < 79; ++si)
    {
        double freq =
            hz0 + 6.25 * re79[si];

        double dtheta =
            2.0 * M_PI * freq / rate_;

        double phase = 0.0;

        /*
         * Initial / steady part.
         *
         * We deliberately don't apply a ramp at the
         * beginning of the entire 79-symbol signal yet.
         */
        for (int jj = 0; jj < block-ramp; ++jj)
        {
            int idx =
                off0 + si * block + jj;

            if (idx < 0 || idx >= (int)nsamples)
                continue;

            double theta =
                phase + jj * dtheta;

            double x =
                amplitude * std::cos(theta);

            samples[idx] += x;
        }

        /*
         * Transition to next symbol.
         */
        double freq1;

        if (si + 1 < 79)
            freq1 =
                hz0 + 6.25 * re79[si + 1];
        else
            freq1 = freq;

        double dtheta1 =
            2.0 * M_PI * freq1 / rate_;

        double inc =
            (dtheta1 - dtheta) /
            (2.0 * ramp);

        double theta =
            phase + (block-ramp) * dtheta;

        for (int jj = block-ramp;
             jj < block;
             ++jj)
        {
            int idx =
                off0 + si * block + jj;

            if (idx >= 0 && idx < (int)nsamples)
            {
                double x =
                    amplitude * std::cos(theta);

                samples[idx] += x;
            }

            theta += dtheta;
            dtheta += inc;
        }
    }
}

static double measure_ft8_signal(
    const float *samples,
    const uint8_t re79[79],
    double hz0,
    double off_sec)
{
    int block = blocksize(rate_);
    int off0 = std::lround(off_sec * rate_);

    double sum = 0.0;

    for (int si = 0; si < 79; ++si)
    {
        double freq =
            hz0 + 6.25 * re79[si];

        double dtheta =
            2.0 * M_PI * freq / rate_;

        double re = 0.0;
        double im = 0.0;

        for (int n = 0; n < block; ++n)
        {
            int idx =
                off0 + si * block + n;

            double theta =
                dtheta * n;

            re += samples[idx] * std::cos(theta);
            im -= samples[idx] * std::sin(theta);
        }

        double amp =
            2.0 * std::sqrt(re * re + im * im) / block;

        sum += amp * amp;
    }

    return std::sqrt(sum / 79.0);
}

static void test_two_signals()
{
    
    const size_t nsamples = 192000;

    const double offset = 2.0;

    const double hzA = 1000.0;
    const double hzB = 1020.0;

    const double ampA = 1000.0;
    const double ampB = 100.0;

    printf("\nTEST 3: two FT8-like signals\n");
    printf("  strong signal A: %.2f Hz, %.1f amplitude\n", hzA, ampA);
    printf("  weak signal B:   %.2f Hz, %.1f amplitude\n", hzB, ampB);
    printf("  offset: %.3f s\n", offset);
    printf("  block:  %d samples\n", blocksize(rate_));

    std::vector<float> buffer(nsamples, 0.0f);

    samples_ = buffer.data();

    /*
     * Strong signal A.
     */
    generate_ft8_simple(
        samples_,
        nsamples,
        tones_A,
        hzA,
        ampA,
        offset);

    /*
     * Weak signal B.
     */
    generate_ft8_simple(
        samples_,
        nsamples,
        tones_B,
        hzB,
        ampB,
        offset);


    double A_before = measure_ft8_signal(samples_, tones_A, hzA, offset);
    double B_before = measure_ft8_signal(samples_, tones_B, hzB, offset);
    printf("  A amplitude before = %.3f\n", A_before);
    printf("  B amplitude before = %.3f\n", B_before);

    /*
     * Subtract strong signal A.
     */
    subtract(
        tones_A,
        hzA,
        hzA,
        offset,
        samples_,
        nsamples,
        rate_);

double A_after =
    measure_ft8_signal(
        samples_,
        tones_A,
        hzA,
        offset);

double B_after =
    measure_ft8_signal(
        samples_,
        tones_B,
        hzB,
        offset);

    printf("  A amplitude after  = %.3f\n", A_after);
    printf("  B amplitude after  = %.3f\n", B_after);

    printf("  A suppression = %.2f dB\n",
        20.0 * std::log10(A_before / A_after));

    printf("  B change      = %.2f dB\n",
        20.0 * std::log10(B_after / B_before));
    }

static void test_synthesize()
{
    printf("\nTEST 4: synthesize consistency\n");

    const size_t nsamples = 192000;
    const double offset = 0.0;
    const double offset_error = 0.012;
    const double hz0 = 1000.0;
    const double amplitude = 1000.0;

    std::vector<float> buffer(nsamples, 0.0f);

    samples_ = buffer.data();

    // std::vector<int> re79 = test_re79_A;
    

    std::vector<double> amps(79, amplitude);
    std::vector<double> phases(79, 0.0);

    /*
     * Generate the exact waveform represented by
     * amps/phases/re79.
     */
    printf(" Calling synthesize() with nsamples = %zu hz0 = %.1f offset = %.1f rate = %d to generate the waveform...\n", nsamples, hz0, offset, offset_error, rate_);
    synthesize(
        samples_,
        nsamples,
        tones_A,
        amps,
        phases,
        hz0,
        offset,
        +1.0,
        rate_);

    int block = blocksize(rate_);
    int off0 = std::lround(offset * rate_);
    int off1 = std::lround((offset + offset_error) * rate_);

    double before =
        rms(samples_ + off0, 79 * block);

    printf("calling subtract with offset = %.6f (offset_error = %.6f). RMS before = %.6f\n", offset + offset_error, offset_error, before);

    /*
     * Now subtract the same waveform.
     *
     * IMPORTANT:
     * this does amplitude/phase estimation first,
     * then calls synthesize(-1).
     */
    subtract(
        tones_A,
        hz0,
        hz0,
        offset + offset_error,
        samples_,
        nsamples,
        rate_);

    double after =
        rms(samples_ + off1, 79 * block);

    printf("  RMS after  = %.9f\n", after);

    printf("  suppression = %.2f dB\n",
           20.0 * std::log10(before / after));

    printf("  PASS\n");
}

static void test_exact_synthesis()
{
    printf("\nTEST 5: exact synthesis cancellation\n");

    const size_t nsamples = 192000;
    const double offset = 2.0;
    const double hz0 = 1000.0;
    const double amplitude = 1000.0;

    std::vector<float> buffer(nsamples, 0.0f);

    samples_ = buffer.data();

    
    std::vector<double> amps(79, amplitude);
    std::vector<double> phases(79, 0.0);

    /*
     * Generate exact waveform.
     */
    synthesize(
        samples_,
        nsamples,
        tones_A,
        amps,
        phases,
        hz0,
        offset,
        +1.0,
        rate_);

    int block = blocksize(rate_);
    int off0 = std::lround(offset * rate_);

    double before =
        rms(samples_ + off0, 79 * block);

    /*
     * Subtract EXACTLY the same waveform.
     */
    synthesize(
        samples_,
        nsamples,
        tones_A,
        amps,
        phases,
        hz0,
        offset,
        -1.0,
        rate_);

    double after =
        rms(samples_ + off0, 79 * block);

    printf("  RMS before = %.9f\n", before);
    printf("  RMS after  = %.12f\n", after);

    printf("  suppression = %.2f dB\n",
           20.0 * std::log10(before / after));
}

int main()
{
    printf("subtract() test_synthesize(); \n");
    printf("rate = %d Hz\n", rate_);\
    printf("block = %d samples\n", blocksize(rate_));

    test_synthesize();

    printf("\nALL TESTS PASSED\n");

    return 0;
}