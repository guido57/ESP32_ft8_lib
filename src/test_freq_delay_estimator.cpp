
// test_frequency_estimator.cpp
//
// Linux test for the FT8 fine-frequency estimator.
//
// Build together with subtract.cpp and the rest of your project.
//
// The synthesizer generates a PHASE-CONTINUOUS FT8 waveform.
// The nominal frequency is 1000 Hz and a frequency error is
// added continuously to the whole signal.
//
// The estimator is then asked to recover:
//
//     true_freq = 1000 Hz + freq_error
//
// Expected result:
//     estimated frequency ~= true_freq
//
// ------------------------------------------------------------

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <random>
#include <complex>
#include "subtract.h"

int rate_ = 12000;

static constexpr double PI = 3.14159265358979323846;

// ------------------------------------------------------------
// FT8 Costas pattern
// ------------------------------------------------------------

static const uint8_t costas[7] =
{
    3, 1, 4, 0, 6, 5, 2
};

// ------------------------------------------------------------
// Build a 79-symbol FT8 test sequence.
//
// We only care about the Costas symbols for the frequency
// estimator. Data symbols are set to tone 0.
//
// Costas positions:
//
//     0 ... 6
//    36 ... 42
//    72 ... 78
// ------------------------------------------------------------

static void make_test_tones(uint8_t tones[79])
{
    for (int i = 0; i < 79; ++i)
        tones[i] = 0;

    for (int i = 0; i < 7; ++i)
    {
        tones[i]      = costas[i];
        tones[36 + i] = costas[i];
        tones[72 + i] = costas[i];
    }
}

// ------------------------------------------------------------
// Phase-continuous FT8 synthesizer
//
// IMPORTANT:
//
// The phase is NOT reset at each symbol.
//
// Actual frequency:
//
//     hz0 + freq_error + 6.25 * tone
//
// The frequency changes at symbol boundaries, but the phase
// continues from the previous symbol.
// ------------------------------------------------------------

static void generate_ft8_continuous(
    float *samples,
    size_t nsamples,
    const uint8_t tones[79],
    double hz0,
    double freq_error,
    double amplitude,
    double off_sec)
{
    const int block = blocksize(rate_);
    const int off0 = std::lround(off_sec * rate_);

    // Global phase.
    // This is deliberately maintained across all 79 symbols.
    double phase = 0.0;

    for (int si = 0; si < 79; ++si)
    {
        const double freq =
            hz0 + freq_error + 6.25 * tones[si];

        const double dphase =
            2.0 * PI * freq / rate_;

        for (int n = 0; n < block; ++n)
        {
            const int idx =
                off0 + si * block + n;

            if (idx >= 0 && idx < (int)nsamples)
            {
                samples[idx] +=
                    (float)(amplitude * std::cos(phase));
            }

            phase += dphase;

            // Keep phase numerically bounded.
            if (phase > 2.0 * PI ||
                phase < -2.0 * PI)
            {
                phase = std::fmod(phase, 2.0 * PI);
            }
        }
    }
}

// ------------------------------------------------------------
// Apply a fractional sample delay to an already generated
// waveform.
//
// Output[n] = input[n - fractional_delay]
//
// Linear interpolation is used only to create the TEST
// stimulus. The estimator itself does not use interpolation.
//
// This is deliberately separate from generate_ft8_continuous()
// because the production synthesizer currently rounds off_sec
// to an integer sample.
// ------------------------------------------------------------

static void apply_fractional_delay(
    const float *input,
    float *output,
    size_t nsamples,
    double fractional_delay)
{
    assert(input != nullptr);
    assert(output != nullptr);

    // fractional_delay must be in [0, 1).
    assert(fractional_delay >= 0.0);
    assert(fractional_delay < 1.0);

    for (size_t n = 0; n < nsamples; ++n)
    {
        const double src =
            (double)n - fractional_delay;

        if (src < 0.0 ||
            src >= (double)(nsamples - 1))
        {
            output[n] = 0.0f;
            continue;
        }

        const size_t i =
            (size_t)std::floor(src);

        const double f =
            src - (double)i;

        output[n] =
            (float)(
                (1.0 - f) * input[i] +
                f * input[i + 1]
            );
    }
}

// ------------------------------------------------------------
// Single estimator test
// ------------------------------------------------------------

struct Result
{
    double freq_error;
    double true_freq;
    double estimated_freq;
    double error;
};

static Result test_one_frequency(
    double freq_error)
{
    const double hz0       = 1000.0;
    const double amplitude = 1000.0;
    const double delay     = 2.0;

    const int block = blocksize(rate_);

    // 15 seconds, same size as a normal FT8 recording.
    const size_t nsamples = 15 * rate_;

    std::vector<float> samples(nsamples, 0.0f);

    uint8_t tones[79];
    make_test_tones(tones);

    generate_ft8_continuous(
        samples.data(),
        samples.size(),
        tones,
        hz0,
        freq_error,
        amplitude,
        delay);

    const double true_freq =
        hz0 + freq_error;

    // The estimator receives the nominal/coarse frequency.
    const float estimated_freq =
        refine_ft8_frequency(
            samples.data(),
            (int)nsamples,
            tones,
            (float)delay,
            (float)hz0);

    const double error =
        estimated_freq - true_freq;

    return {
        freq_error,
        true_freq,
        estimated_freq,
        error
    };
}

static void test_frequency_estimator()
{
    printf("\n");
    printf("============================================================\n");
    printf("TEST A: FT8 FINE FREQUENCY ESTIMATOR - NOISELESS\n");
    printf("============================================================\n");

    printf("Sample rate       : %d Hz\n", rate_);
    printf("Symbol block      : %d samples\n", blocksize(rate_));
    printf("Symbol duration   : %.3f ms\n",
           1000.0 * blocksize(rate_) / (double)rate_);
    printf("Tone spacing      : 6.25 Hz\n");
    printf("Nominal frequency : 1000.000 Hz\n");
    printf("Signal delay      : 2.000 s\n");
    printf("Amplitude         : 1000\n");
    printf("Noise             : NONE\n");

    printf("\n");
    printf("%10s %12s %15s %12s\n",
           "Offset",
           "True freq",
           "Estimated",
           "Error");
    printf("------------------------------------------------------------\n");

    std::vector<double> errors;

    // Sweep -2.5 ... +2.5 Hz in 0.25 Hz steps.
    for (int i = -10; i <= 10; ++i)
    {
        const double freq_error =
            0.25 * i;

        Result r =
            test_one_frequency(freq_error);

        errors.push_back(r.error);

        printf("%+10.2f %12.3f %15.6f %+12.6f\n",
               r.freq_error,
               r.true_freq,
               r.estimated_freq,
               r.error);
    }

    // --------------------------------------------------------
    // Statistics
    // --------------------------------------------------------

    double sum = 0.0;
    double sum2 = 0.0;
    double max_abs = 0.0;

    for (double e : errors)
    {
        sum += e;
        sum2 += e * e;
        max_abs = std::max(max_abs, std::abs(e));
    }

    const double mean =
        sum / errors.size();

    const double rms =
        std::sqrt(sum2 / errors.size());

    printf("------------------------------------------------------------\n");
    printf("Mean error        : %+12.6f Hz\n", mean);
    printf("RMS error         : %12.6f Hz\n", rms);
    printf("Maximum error     : %12.6f Hz\n", max_abs);
    printf("------------------------------------------------------------\n");

    // --------------------------------------------------------
    // Acceptance criteria
    //
    // First objective: verify that the estimator itself is
    // unbiased in a perfect noiseless signal.
    //
    // 0.001 Hz is deliberately not used as a hard requirement
    // yet. We first want to see the actual result.
    // --------------------------------------------------------

    constexpr double MAX_ALLOWED_ERROR = 0.01; // 10 mHz

    if (max_abs <= MAX_ALLOWED_ERROR)
    {
        printf("RESULT            : PASS\n");
        printf("Maximum error is below %.3f Hz\n",
               MAX_ALLOWED_ERROR);
    }
    else
    {
        printf("RESULT            : FAIL\n");
        printf("Maximum error exceeds %.3f Hz\n",
               MAX_ALLOWED_ERROR);
    }

    printf("============================================================\n");
}
// ------------------------------------------------------------
// Optional visual/debug test:
//
// Generate one waveform with +1 Hz error.
//
// This can be useful for inspecting the generated WAV later.
// ------------------------------------------------------------

static void test_generate_one()
{
    const double hz0       = 1000.0;
    const double freq_error = +1.0;
    const double amplitude = 1000.0;
    const double delay     = 2.0;

    const size_t nsamples =
        15 * rate_;

    std::vector<float> samples(
        nsamples, 0.0f);

    uint8_t tones[79];
    make_test_tones(tones);

    generate_ft8_continuous(
        samples.data(),
        samples.size(),
        tones,
        hz0,
        freq_error,
        amplitude,
        delay);

    printf("\n");
    printf("Generated test signal:\n");
    printf("  nominal frequency = %.3f Hz\n", hz0);
    printf("  frequency error   = %+.3f Hz\n", freq_error);
    printf("  actual frequency  = %.3f Hz\n",
           hz0 + freq_error);
    printf("  delay              = %.3f s\n", delay);
}


// ------------------------------------------------------------
// TEST B
//
// Frequency estimator accuracy versus AWGN.
//
// The signal is:
//
//     nominal frequency = 1000 Hz
//     true frequency    = 1001 Hz
//     frequency error   = +1 Hz
//     delay             = 2.0 s
//
// Independent AWGN is added for every trial.
//
// SNR is defined using the RMS power of the generated FT8
// signal over its active 79-symbol interval.
//
// ------------------------------------------------------------

struct NoiseTestResult
{
    double snr_db;
    double mean_error;
    double rms_error;
    double max_abs_error;
    int failures;
    int trials;
};


// ------------------------------------------------------------
// Generate Gaussian noise with the requested RMS.
//
// sigma is the standard deviation of the noise.
// ------------------------------------------------------------

static void add_awgn(
    float *samples,
    size_t nsamples,
    double sigma,
    std::mt19937 &rng)
{
    std::normal_distribution<double> normal(0.0, sigma);

    for (size_t i = 0; i < nsamples; ++i)
        samples[i] += (float)normal(rng);
}


// ------------------------------------------------------------
// Calculate signal power over the active FT8 waveform.
//
// The signal starts at delay and occupies 79 symbols.
//
// We deliberately do NOT include the initial silence in the
// SNR calculation.
// ------------------------------------------------------------

static double calculate_signal_power(
    const float *samples,
    size_t nsamples,
    double delay)
{
    const int block = blocksize(rate_);

    const int start =
        std::lround(delay * rate_);

    const size_t active_samples =
        std::min(
            (size_t)(79 * block),
            nsamples > (size_t)start
                ? nsamples - (size_t)start
                : 0);

    if (active_samples == 0)
        return 0.0;

    double sum2 = 0.0;

    for (size_t i = 0; i < active_samples; ++i)
    {
        const double x =
            samples[start + i];

        sum2 += x * x;
    }

    return sum2 / (double)active_samples;
}


// ------------------------------------------------------------
// Run one SNR point.
// ------------------------------------------------------------

static NoiseTestResult test_one_snr(
    double snr_db,
    int trials,
    std::mt19937 &rng)
{
    const double hz0        = 1000.0;
    const double freq_error = +1.0;
    const double amplitude  = 1000.0;
    const double delay      = 2.0;

    const size_t nsamples =
        15 * rate_;

    uint8_t tones[79];
    make_test_tones(tones);

    // --------------------------------------------------------
    // Generate a clean reference signal once.
    // --------------------------------------------------------

    std::vector<float> clean(
        nsamples, 0.0f);

    generate_ft8_continuous(
        clean.data(),
        clean.size(),
        tones,
        hz0,
        freq_error,
        amplitude,
        delay);

    // --------------------------------------------------------
    // Signal power and corresponding noise power.
    //
    // SNR = signal_power / noise_power
    // --------------------------------------------------------

    const double signal_power =
        calculate_signal_power(
            clean.data(),
            clean.size(),
            delay);

    const double noise_power =
        signal_power /
        std::pow(10.0, snr_db / 10.0);

    const double noise_sigma =
        std::sqrt(noise_power);

    // --------------------------------------------------------
    // Run trials.
    // --------------------------------------------------------

    double sum_error  = 0.0;
    double sum_error2 = 0.0;
    double max_abs    = 0.0;

    int failures = 0;

    for (int trial = 0; trial < trials; ++trial)
    {
        // Start from the identical clean signal every time.
        std::vector<float> samples = clean;

        add_awgn(
            samples.data(),
            samples.size(),
            noise_sigma,
            rng);

        const float estimated_freq =
            refine_ft8_frequency(
                samples.data(),
                (int)samples.size(),
                tones,
                (float)delay,
                (float)hz0);

        if (!std::isfinite(estimated_freq))
        {
            ++failures;
            continue;
        }

        const double error =
            estimated_freq -
            (hz0 + freq_error);

        sum_error  += error;
        sum_error2 += error * error;

        max_abs =
            std::max(max_abs, std::abs(error));
    }

    const int valid =
        trials - failures;

    double mean = NAN;
    double rms  = NAN;

    if (valid > 0)
    {
        mean =
            sum_error / valid;

        rms =
            std::sqrt(
                sum_error2 / valid);
    }

    return {
        snr_db,
        mean,
        rms,
        max_abs,
        failures,
        trials
    };
}


// ------------------------------------------------------------
// Complete Test B
// ------------------------------------------------------------

static void test_frequency_estimator_noise()
{
    printf("\n");
    printf("============================================================\n");
    printf("TEST B: FT8 FINE FREQUENCY ESTIMATOR - AWGN\n");
    printf("============================================================\n");

    printf("Sample rate       : %d Hz\n", rate_);
    printf("Symbol block      : %d samples\n",
           blocksize(rate_));
    printf("Symbol duration   : %.3f ms\n",
           1000.0 * blocksize(rate_) / (double)rate_);
    printf("Tone spacing      : 6.25 Hz\n");
    printf("Nominal frequency : 1000.000 Hz\n");
    printf("True frequency    : 1001.000 Hz\n");
    printf("Frequency offset  : +1.000 Hz\n");
    printf("Signal delay      : 2.000 s\n");
    printf("Amplitude         : 1000\n");
    printf("Trials / SNR      : 100\n");
    printf("Noise             : Gaussian AWGN\n");

    printf("\n");
    printf("%8s %14s %14s %14s %10s %10s\n",
           "SNR",
           "Mean error",
           "RMS error",
           "Max error",
           "Failures",
           "Trials");

    printf("------------------------------------------------------------\n");

    // --------------------------------------------------------
    // Fixed seed makes the test completely reproducible.
    // --------------------------------------------------------

    std::mt19937 rng(123456789);

    constexpr int trials = 100;

    const double snr_values[] =
    {
        +10.0,
         +5.0,
          0.0,
         -5.0,
        -10.0,
        -15.0,
        -20.0,
        -25.0,
        -30.0
    };

    for (double snr_db : snr_values)
    {
        NoiseTestResult r =
            test_one_snr(
                snr_db,
                trials,
                rng);

        printf("%+8.1f %14.6f %14.6f %14.6f %10d %10d\n",
               r.snr_db,
               r.mean_error,
               r.rms_error,
               r.max_abs_error,
               r.failures,
               r.trials);
    }

    printf("------------------------------------------------------------\n");
    printf("Test B completed.\n");
    printf("============================================================\n");
}

// ------------------------------------------------------------
// TEST C
//
// FT8 delay estimator accuracy versus AWGN.
//
// The estimator has an exact resolution of:
//
//     1 / 12000 = 83.333 us
//
// Therefore this test uses delays that are EXACTLY aligned
// to integer sample positions.
//
// This isolates the ability of the estimator to select the
// correct sample in the presence of noise.
//
// The key metric is:
//
//     Wrong samples
//
// i.e. the number of trials where:
//
//     estimated_delay_samples != true_delay_samples
//
// ------------------------------------------------------------

struct DelayTestResult
{
    double snr_db;
    double mean_error_ms;
    double rms_error_ms;
    double max_abs_error_ms;
    int wrong_samples;
    int failures;
    int trials;
};


// ------------------------------------------------------------
// Run one SNR point
// ------------------------------------------------------------

static DelayTestResult test_one_delay_snr(
    double snr_db,
    int trials,
    std::mt19937 &rng)
{
    const double hz0        = 1000.0;
    const double freq_error = +1.0;
    const double amplitude  = 1000.0;

    // --------------------------------------------------------
    // IMPORTANT:
    //
    // True delay is EXACTLY an integer number of samples.
    //
    // 24000 samples at 12000 Hz = exactly 2.0 seconds.
    //
    // No fractional-delay ambiguity exists.
    // --------------------------------------------------------

    const int true_delay_samples = 24000;

    const double true_delay =
        true_delay_samples /
        (double)rate_;

    const size_t nsamples =
        15 * rate_;

    uint8_t tones[79];
    make_test_tones(tones);

    // --------------------------------------------------------
    // Generate clean signal.
    // --------------------------------------------------------

    std::vector<float> clean(
        nsamples, 0.0f);

    generate_ft8_continuous(
        clean.data(),
        clean.size(),
        tones,
        hz0,
        freq_error,
        amplitude,
        true_delay);

    // --------------------------------------------------------
    // Signal power.
    // --------------------------------------------------------

    const double signal_power =
        calculate_signal_power(
            clean.data(),
            clean.size(),
            true_delay);

    const double noise_power =
        signal_power /
        std::pow(10.0, snr_db / 10.0);

    const double noise_sigma =
        std::sqrt(noise_power);

    // --------------------------------------------------------
    // Statistics.
    // --------------------------------------------------------

    double sum_error  = 0.0;
    double sum_error2 = 0.0;
    double max_abs    = 0.0;

    int wrong_samples = 0;
    int failures = 0;

    for (int trial = 0;
         trial < trials;
         ++trial)
    {
        std::vector<float> samples = clean;

        add_awgn(
            samples.data(),
            samples.size(),
            noise_sigma,
            rng);

        // ----------------------------------------------------
        // Give the estimator the exact coarse delay.
        //
        // This isolates delay refinement from coarse
        // candidate timing error.
        // ----------------------------------------------------

        const float estimated_delay =
            refine_ft8_delay(
                samples.data(),
                (int)samples.size(),
                tones,
                (float)true_delay,
                (float)hz0,
                -1);

        if (!std::isfinite(estimated_delay))
        {
            ++failures;
            continue;
        }

        // ----------------------------------------------------
        // Convert both delays to sample indices.
        //
        // The estimator itself has integer-sample resolution.
        // ----------------------------------------------------

        const int estimated_delay_samples =
            (int)std::lround(
                estimated_delay * rate_);

        if (estimated_delay_samples !=
            true_delay_samples)
        {
            ++wrong_samples;
        }

        // ----------------------------------------------------
        // Error in milliseconds.
        // ----------------------------------------------------

        const double error_ms =
            1000.0 *
            ((double)estimated_delay -
             true_delay);

        sum_error += error_ms;

        sum_error2 +=
            error_ms * error_ms;

        max_abs =
            std::max(
                max_abs,
                std::abs(error_ms));
    }

    const int valid =
        trials - failures;

    double mean = NAN;
    double rms  = NAN;

    if (valid > 0)
    {
        mean =
            sum_error / valid;

        rms =
            std::sqrt(
                sum_error2 / valid);
    }

    return {
        snr_db,
        mean,
        rms,
        max_abs,
        wrong_samples,
        failures,
        trials
    };
}


// ------------------------------------------------------------
// Complete Test C
// ------------------------------------------------------------

static void test_delay_estimator_noise()
{
    printf("\n");
    printf("============================================================\n");
    printf("TEST C: FT8 DELAY ESTIMATOR - AWGN\n");
    printf("============================================================\n");

    printf("Sample rate       : %d Hz\n", rate_);
    printf("Symbol block      : %d samples\n",
           blocksize(rate_));

    printf("Symbol duration   : %.3f ms\n",
           1000.0 *
           blocksize(rate_) /
           (double)rate_);

    printf("Tone spacing      : 6.25 Hz\n");

    printf("Nominal frequency : 1000.000 Hz\n");
    printf("True frequency    : 1001.000 Hz\n");

    printf("True delay        : %.9f s\n",
           24000.0 / rate_);

    printf("True delay sample  : %d\n",
           24000);

    printf("Delay resolution  : %.6f ms\n",
           1000.0 / rate_);

    printf("Amplitude         : 1000\n");
    printf("Trials / SNR      : 100\n");
    printf("Noise             : Gaussian AWGN\n");

    printf("\n");

    printf("%8s %14s %14s %14s %14s %10s %10s\n",
           "SNR",
           "Mean error",
           "RMS error",
           "Max error",
           "Wrong samples",
           "Failures",
           "Trials");

    printf("----------------------------------------------------------------------------\n");

    // Different seed from Test B.
    std::mt19937 rng(987654321);

    constexpr int trials = 100;

    const double snr_values[] =
    {
        +10.0,
         +5.0,
          0.0,
         -5.0,
        -10.0,
        -15.0,
        -20.0,
        -25.0,
        -30.0
    };

    for (double snr_db : snr_values)
    {
        DelayTestResult r =
            test_one_delay_snr(
                snr_db,
                trials,
                rng);

        printf("%+8.1f %14.6f %14.6f %14.6f %14d %10d %10d\n",
               r.snr_db,
               r.mean_error_ms,
               r.rms_error_ms,
               r.max_abs_error_ms,
               r.wrong_samples,
               r.failures,
               r.trials);
    }

    printf("----------------------------------------------------------------------------\n");

    printf("Delay resolution:\n");
    printf("  1 sample = %.6f ms\n",
           1000.0 / rate_);

    printf("============================================================\n");
    printf("Test C completed.\n");
    printf("============================================================\n");
}
static void test_delay_estimator_noiseless()
{
    printf("\n");
    printf("============================================================\n");
    printf("TEST C0: FT8 DELAY ESTIMATOR - NOISELESS\n");
    printf("============================================================\n");

    const double hz0        = 1000.0;
    const double freq_error = +1.0;
    const double amplitude  = 1000.0;

    const size_t nsamples =
        15 * rate_;

    uint8_t tones[79];
    make_test_tones(tones);

    // --------------------------------------------------------
    // IMPORTANT:
    //
    // All delays are specified as EXACT integer sample
    // positions.
    //
    // At 12000 Hz:
    //
    //     24000 samples = 2.000000 s
    //     24001 samples = 2.000083 s
    //     ...
    //
    // This avoids any ambiguity caused by fractional-sample
    // delays.
    // --------------------------------------------------------

    const int delay_samples[] =
    {
        24000,
        24001,
        24002,
        24005,
        24010,
        24020,
        24030,
        24040,
        24050,
        24060,
        24070,
        24080,
        24090,
        24100
    };

    printf("\n");
    printf("%8s %12s %15s %15s\n",
           "True samp",
           "True delay",
           "Estimated",
           "Error ms");

    printf("------------------------------------------------------------\n");

    for (int true_delay_samples : delay_samples)
    {
        const double true_delay =
            true_delay_samples / (double)rate_;

        // ----------------------------------------------------
        // Generate clean signal.
        // ----------------------------------------------------

        std::vector<float> samples(
            nsamples, 0.0f);

        generate_ft8_continuous(
            samples.data(),
            samples.size(),
            tones,
            hz0,
            freq_error,
            amplitude,
            true_delay);

        // ----------------------------------------------------
        // Estimate delay.
        //
        // Give the estimator the exact true delay as delay0.
        // This isolates the delay refinement itself.
        // ----------------------------------------------------

        const float estimated_delay =
            refine_ft8_delay(
                samples.data(),
                (int)samples.size(),
                tones,
                (float)true_delay,
                (float)hz0,
                -1);

        // ----------------------------------------------------
        // Convert estimated delay back to sample number.
        // ----------------------------------------------------

        const int estimated_delay_samples =
            (int)std::lround(
                estimated_delay * rate_);

        // ----------------------------------------------------
        // Error in milliseconds.
        // ----------------------------------------------------

        const double error_ms =
            1000.0 *
            ((double)estimated_delay -
             true_delay);

        printf("%8d %12.6f %15d %+15.6f\n",
               true_delay_samples,
               true_delay,
               estimated_delay_samples,
               error_ms);
    }

    printf("------------------------------------------------------------\n");

    printf("Delay resolution: %.6f ms/sample\n",
           1000.0 / rate_);

    printf("Test C0 completed.\n");
    printf("============================================================\n");
}
// ------------------------------------------------------------
// ------------------------------------------------------------
// DEBUG: inspect the delay score curve
//
// This reproduces the score calculation used by
// refine_ft8_delay(), but prints the score around the
// true delay.
//
// No modification of the actual estimator.
// ------------------------------------------------------------

static void debug_delay_score_curve(
    const float* samples,
    int num_samples,
    const uint8_t* tones,
    double delay0,
    double freq,
    int radius_samples)
{
    constexpr int sample_rate = 12000;
    constexpr int ntones = 79;

    const int block = blocksize(sample_rate);

    constexpr double TWO_PI =
        6.2831853071795864769;

    const int center =
        (int)std::lround(delay0 * sample_rate);

    const int first_delay =
        std::max(0, center - radius_samples);

    const int last_delay =
        std::min(
            num_samples - 1,
            center + radius_samples);

    const int ndelays =
        last_delay - first_delay + 1;

    std::vector<double> score_total(
        ndelays, 0.0);

    std::vector<int> used_count(
        ndelays, 0);

    // --------------------------------------------------------
    // Calculate score exactly as refine_ft8_delay().
    // --------------------------------------------------------

    for (int k = 0; k < ntones; ++k)
    {
        const double freq_hz =
            freq + 6.25 * tones[k];

        const double w =
            TWO_PI * freq_hz /
            (double)sample_rate;

        const double cd =
            std::cos(w);

        const double sd =
            std::sin(w);

        double end_c = 1.0;
        double end_s = 0.0;

        for (int n = 0; n < block - 1; ++n)
        {
            const double nc =
                end_c * cd + end_s * sd;

            const double ns =
                end_s * cd - end_c * sd;

            end_c = nc;
            end_s = ns;
        }

        double ci = 0.0;
        double cq = 0.0;

        const int start0 =
            first_delay + k * block;

        if (start0 >= 0 &&
            start0 < num_samples)
        {
            const int count =
                std::min(
                    block,
                    num_samples - start0);

            double c = 1.0;
            double s = 0.0;

            for (int n = 0; n < count; ++n)
            {
                const double x =
                    samples[start0 + n];

                ci += x * c;
                cq -= x * s;

                const double nc =
                    c * cd -
                    s * sd;

                const double ns =
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

            const double x_old =
                samples[start - 1];

            const double x_new =
                samples[new_index];

            const double r =
                ci - x_old;

            const double im =
                cq;

            const double rotated_re =
                r * cd - im * sd;

            const double rotated_im =
                r * sd + im * cd;

            ci =
                rotated_re +
                x_new * end_c;

            cq =
                rotated_im +
                x_new * end_s;

            score_total[d] +=
                ci * ci + cq * cq;

            ++used_count[d];
        }
    }

    // --------------------------------------------------------
    // Find maximum.
    // --------------------------------------------------------

    int best_d = 0;
    double best_score = -1.0;

    for (int d = 0; d < ndelays; ++d)
    {
        if (used_count[d] == 0)
            continue;

        const double score =
            score_total[d] /
            (double)used_count[d];

        if (score > best_score)
        {
            best_score = score;
            best_d = d;
        }
    }

    const int best_sample =
        first_delay + best_d;

    // --------------------------------------------------------
    // Print score curve.
    //
    // Normalize to the maximum so the shape is easier to see.
    // --------------------------------------------------------

    printf("\n");
    printf("============================================================\n");
    printf("DELAY SCORE CURVE\n");
    printf("============================================================\n");

    printf("True delay sample : %d\n", center);
    printf("Best delay sample : %d\n", best_sample);
    printf("True delay        : %.6f ms\n",
           1000.0 * center / sample_rate);
    printf("Best delay        : %.6f ms\n",
           1000.0 * best_sample / sample_rate);

    printf("\n");
    printf("%8s %14s %14s\n",
           "Offset",
           "Score",
           "Relative");

    printf("-----------------------------------------------\n");

    for (int d = 0; d < ndelays; ++d)
    {
        if (used_count[d] == 0)
            continue;

        const double score =
            score_total[d] /
            (double)used_count[d];

        const int sample =
            first_delay + d;

        const int offset =
            sample - center;

        const double relative =
            score / best_score;

        printf("%+8d %14.3f %14.6f",
               offset,
               score,
               relative);

        if (sample == center)
            printf("  <-- TRUE");

        if (sample == best_sample)
            printf("  <-- BEST");

        printf("\n");
    }

    printf("============================================================\n");
}

// ------------------------------------------------------------
// DEBUG TEST
//
// Generate one +10 dB noisy signal and inspect the delay
// score around the true delay.
// ------------------------------------------------------------

static void test_delay_score_debug()
{
    printf("\n");
    printf("============================================================\n");
    printf("TEST C DEBUG: DELAY SCORE AT +10 dB\n");
    printf("============================================================\n");

    const double hz0        = 1000.0;
    const double freq_error = +1.0;
    const double amplitude  = 1000.0;

    const int true_delay_samples = 24000;

    const double true_delay =
        true_delay_samples /
        (double)rate_;

    const double snr_db = +10.0;

    const size_t nsamples =
        15 * rate_;

    uint8_t tones[79];
    make_test_tones(tones);

    std::vector<float> samples(
        nsamples, 0.0f);

    generate_ft8_continuous(
        samples.data(),
        samples.size(),
        tones,
        hz0,
        freq_error,
        amplitude,
        true_delay);

    const double signal_power =
        calculate_signal_power(
            samples.data(),
            samples.size(),
            true_delay);

    const double noise_power =
        signal_power /
        std::pow(10.0, snr_db / 10.0);

    const double noise_sigma =
        std::sqrt(noise_power);

    std::mt19937 rng(123456789);

    add_awgn(
        samples.data(),
        samples.size(),
        noise_sigma,
        rng);

    printf("SNR               : %.1f dB\n", snr_db);
    printf("True delay sample : %d\n",
           true_delay_samples);
    printf("Noise sigma       : %.3f\n",
           noise_sigma);

    debug_delay_score_curve(
        samples.data(),
        (int)samples.size(),
        tones,
        true_delay,
        hz0,
        15);

    const float estimated_delay =
        refine_ft8_delay(
            samples.data(),
            (int)samples.size(),
            tones,
            (float)true_delay,
            (float)hz0,
            -1);

    const int estimated_sample =
        (int)std::lround(
            estimated_delay * rate_);

    printf("\nEstimator result:\n");
    printf("  estimated delay  = %.9f s\n",
           estimated_delay);
    printf("  estimated sample = %d\n",
           estimated_sample);
    printf("  error             = %.6f ms\n",
           1000.0 *
           (estimated_delay - true_delay));

    printf("============================================================\n");
}


// ============================================================
// D1 v6
//
// FT8 delay estimation:
//
//     RF
//      |
//      v
// complex downconversion
//      |
//      v
// low-pass filter
//      |
//      v
// Costas transition correlation
//
// The low-pass filter removes the unwanted image produced by
// mixing a REAL RF signal:
//
//     cos(fc) * exp(-j*fc)
//              |
//              +-- desired baseband
//              |
//              +-- image at 2*fc
//
// At fc = 1001 Hz the image is at 2002 Hz, corresponding to
// approximately 6 samples/cycle at Fs = 12000 Hz.
//
// Reference waveform is generated by the EXACT same
// generate_ft8_continuous() function used by the production
// code.
//
// No independent phase reconstruction.
// No FFT.
// No std::complex.
// ============================================================

static void test_D1()
{
    printf("\n");
    printf("============================================================\n");
    printf("D1 v6 - FT8 BASEBAND COSTAS TRANSITION DELAY ESTIMATOR\n");
    printf("============================================================\n");

    constexpr int SAMPLE_RATE = 12000;
    constexpr int NSYM = 79;
    constexpr int SPS = 1920;

    constexpr double FREQ = 1001.0;
    constexpr double TONE_SPACING = 6.25;

    constexpr int TRUE_DELAY = 24000;
    constexpr double AMPLITUDE = 1000.0;

    constexpr int SEARCH = 15;

    // --------------------------------------------------------
    // Transition window
    // --------------------------------------------------------

    constexpr int WIN = 128;
    constexpr int HALF_WIN = WIN / 2;

    // --------------------------------------------------------
    // FIR low-pass filter
    //
    // Passband comfortably contains all FT8 tones:
    //
    //     maximum tone = 7 * 6.25 = 43.75 Hz
    //
    // Cutoff = 150 Hz.
    //
    // The purpose here is primarily to remove the 2*carrier
    // image (~2002 Hz).
    //
    // Odd number of taps -> integer group delay.
    // --------------------------------------------------------

    constexpr int FIR_TAPS = 129;
    constexpr double LPF_CUTOFF = 150.0;

    // --------------------------------------------------------
    // Generate test tones.
    // --------------------------------------------------------

    uint8_t tones[NSYM];
    make_test_tones(tones);

    const size_t total_samples =
        TRUE_DELAY + NSYM * SPS;

    std::vector<float> clean(total_samples);
    std::vector<float> noisy(total_samples);

    // --------------------------------------------------------
    // Generate the exact production waveform.
    // --------------------------------------------------------

    generate_ft8_continuous(
        clean.data(),
        total_samples,
        tones,
        FREQ,
        0.0,
        AMPLITUDE,
        TRUE_DELAY / (double)SAMPLE_RATE
    );

    // --------------------------------------------------------
    // Received signal.
    //
    // Identity test for now.
    //
    // AWGN will be added in the next test.
    // --------------------------------------------------------

    for (size_t i = 0; i < total_samples; ++i)
        noisy[i] = clean[i];

    // ========================================================
    // DESIGN FIR LOW-PASS
    // ========================================================

    float fir[FIR_TAPS];

    const double fc =
        LPF_CUTOFF / SAMPLE_RATE;

    const int mid =
        FIR_TAPS / 2;

    double sum = 0.0;

    for (int k = 0; k < FIR_TAPS; ++k)
    {
        const int n = k - mid;

        double h;

        if (n == 0)
        {
            h = 2.0 * fc;
        }
        else
        {
            h =
                sin(2.0 * M_PI * fc * n) /
                (M_PI * n);
        }

        // Hann window.
        const double window =
            0.5 *
            (1.0 -
             cos(2.0 * M_PI * k /
                 (FIR_TAPS - 1)));

        h *= window;

        fir[k] = (float)h;
        sum += h;
    }

    // Normalize DC gain to 1.
    for (int k = 0; k < FIR_TAPS; ++k)
        fir[k] /= (float)sum;

    printf("\n");
    printf("Sample rate        : %d Hz\n", SAMPLE_RATE);
    printf("Carrier            : %.3f Hz\n", FREQ);
    printf("LPF cutoff         : %.1f Hz\n", LPF_CUTOFF);
    printf("FIR taps           : %d\n", FIR_TAPS);
    printf("FIR group delay    : %d samples\n", mid);

    // ========================================================
    // COMPLEX DOWNCONVERSION
    //
    // x[n] * exp(-j*carrier)
    //
    // I = x*cos
    // Q = -x*sin
    //
    // We do this for the COMPLETE waveform before filtering.
    // ========================================================

    std::vector<float> clean_i(total_samples);
    std::vector<float> clean_q(total_samples);

    std::vector<float> noisy_i(total_samples);
    std::vector<float> noisy_q(total_samples);

    const double w =
        2.0 * M_PI * FREQ / SAMPLE_RATE;

    const float cw = cos(w);
    const float sw = sin(w);

    float c = 1.0f;
    float s = 0.0f;

    for (size_t n = 0; n < total_samples; ++n)
    {
        const float x = clean[n];

        clean_i[n] =
            x * c;

        clean_q[n] =
            -x * s;

        const float nc =
            c * cw -
            s * sw;

        const float ns =
            s * cw +
            c * sw;

        c = nc;
        s = ns;
    }

    // Reset oscillator for received signal.
    c = 1.0f;
    s = 0.0f;

    for (size_t n = 0; n < total_samples; ++n)
    {
        const float x = noisy[n];

        noisy_i[n] =
            x * c;

        noisy_q[n] =
            -x * s;

        const float nc =
            c * cw -
            s * sw;

        const float ns =
            s * cw +
            c * sw;

        c = nc;
        s = ns;
    }

    // ========================================================
    // FIR FILTER
    //
    // Apply exactly the same filter to clean/reference and
    // received signal.
    // ========================================================

    std::vector<float> clean_i_lp(total_samples);
    std::vector<float> clean_q_lp(total_samples);

    std::vector<float> noisy_i_lp(total_samples);
    std::vector<float> noisy_q_lp(total_samples);

    for (size_t n = 0; n < total_samples; ++n)
    {
        double ai = 0.0;
        double aq = 0.0;

        double bi = 0.0;
        double bq = 0.0;

        const int first =
            (n >= FIR_TAPS - 1)
                ? (int)n - FIR_TAPS + 1
                : 0;

        for (int k = first; k <= (int)n; ++k)
        {
            const int h =
                (int)n - k;

            ai +=
                clean_i[k] *
                fir[h];

            aq +=
                clean_q[k] *
                fir[h];

            bi +=
                noisy_i[k] *
                fir[h];

            bq +=
                noisy_q[k] *
                fir[h];
        }

        clean_i_lp[n] = (float)ai;
        clean_q_lp[n] = (float)aq;

        noisy_i_lp[n] = (float)bi;
        noisy_q_lp[n] = (float)bq;
    }

    // ========================================================
    // COSTAS TRANSITIONS
    // ========================================================

    struct Transition
    {
        int sample;
        int tone_before;
        int tone_after;
    };

    Transition transitions[22];

    int ntrans = 0;

    for (int k = 0; k < NSYM - 1; ++k)
    {
        if (tones[k] == tones[k + 1])
            continue;

        const int boundary =
            TRUE_DELAY + (k + 1) * SPS;

        transitions[ntrans++] =
        {
            boundary,
            tones[k],
            tones[k + 1]
        };
    }

    printf("\n");
    printf("Costas transitions : %d\n", ntrans);
    printf("Window             : %d samples\n", WIN);
    printf("Search             : +/- %d samples\n", SEARCH);

    // ========================================================
    // EXTRACT FILTERED BASEBAND REFERENCE
    //
    // This is taken directly from the generated waveform.
    //
    // Therefore there is NO independent phase model.
    // ========================================================

    float ref_i[22][WIN];
    float ref_q[22][WIN];

    for (int t = 0; t < ntrans; ++t)
    {
        const int center =
            transitions[t].sample;

        const int start =
            center - HALF_WIN;

        for (int m = 0; m < WIN; ++m)
        {
            const int idx =
                start + m;

            ref_i[t][m] =
                clean_i_lp[idx];

            ref_q[t][m] =
                clean_q_lp[idx];
        }
    }

    // ========================================================
    // SEARCH DELAY
    //
    // NORMALIZED COMPLEX CORRELATION
    //
    //        |sum(received * conj(reference))|^2
    // rho2 = -----------------------------------
    //        sum(|received|^2) * sum(|reference|^2)
    //
    // This removes the dependence on the energy of the
    // candidate window.
    // ========================================================

    double scores[2 * SEARCH + 1];

    double best_score = -1.0;
    int best_offset = 0;

    for (int d = -SEARCH; d <= SEARCH; ++d)
    {
        double score = 0.0;

        for (int t = 0; t < ntrans; ++t)
        {
            const int center =
                transitions[t].sample + d;

            const int start =
                center - HALF_WIN;

            double ci = 0.0;
            double cq = 0.0;

            double energy_received = 0.0;
            double energy_reference = 0.0;

            for (int m = 0; m < WIN; ++m)
            {
                const int idx =
                    start + m;

                const float bi =
                    noisy_i_lp[idx];

                const float bq =
                    noisy_q_lp[idx];

                const float ri =
                    ref_i[t][m];

                const float rq =
                    ref_q[t][m];

                // ------------------------------------------------
                // received * conj(reference)
                // ------------------------------------------------

                ci +=
                    bi * ri +
                    bq * rq;

                cq +=
                    bq * ri -
                    bi * rq;

                // ------------------------------------------------
                // Energies
                // ------------------------------------------------

                energy_received +=
                    bi * bi +
                    bq * bq;

                energy_reference +=
                    ri * ri +
                    rq * rq;
            }

            const double correlation_power =
                ci * ci +
                cq * cq;

            const double denominator =
                energy_received *
                energy_reference;

            double rho2 = 0.0;

            if (denominator > 0.0)
            {
                rho2 =
                    correlation_power /
                    denominator;
            }

            score += rho2;
        }

        scores[d + SEARCH] =
            score;

        if (score > best_score)
        {
            best_score = score;
            best_offset = d;
        }
    }
    
    // ========================================================
    // RESULT
    // ========================================================

    printf("\n");
    printf("True delay      : %.6f ms\n",
           1000.0 *
           TRUE_DELAY /
           SAMPLE_RATE);

    printf("Best offset     : %+d samples\n",
           best_offset);

    printf("Estimated delay : %.6f ms\n",
           1000.0 *
           (TRUE_DELAY + best_offset) /
           SAMPLE_RATE);

    printf("\n");
    printf("Offset       Relative score\n");
    printf("--------------------------------\n");

    for (int d = -SEARCH; d <= SEARCH; ++d)
    {
        const double normalized =
            scores[d + SEARCH] /
            best_score;

        printf("%+5d       %.9f\n",
               d,
               normalized);
    }

    printf("\n");

    if (best_offset == 0)
    {
        printf("D1 PASS: exact delay found\n");
    }
    else
    {
        printf("D1 FAIL: delay error = %+d samples\n",
               best_offset);
    }

    printf("============================================================\n");
}

// ============================================================
// D2
//
// FT8 DELAY ESTIMATION WITH FRACTIONAL DELAY
//
// The signal is deliberately delayed by:
//
//     integer_delay + fractional_delay
//
// The estimator performs ONLY an integer-sample search.
//
// After finding the integer maximum, a 3-point parabolic
// interpolation estimates the sub-sample position.
//
// The expensive correlation search therefore remains:
//
//     +/- SEARCH samples
//
// with exactly one correlation calculation per integer offset.
//
// ============================================================

static void test_D2()
{
    printf("\n");
    printf("============================================================\n");
    printf("D2 - FT8 FRACTIONAL DELAY + PARABOLIC INTERPOLATION\n");
    printf("============================================================\n");

    constexpr int SAMPLE_RATE = 12000;
    constexpr int NSYM = 79;
    constexpr int SPS = 1920;

    constexpr double FREQ = 1001.0;
    constexpr double AMPLITUDE = 1000.0;

    constexpr int INTEGER_DELAY = 24000;

    constexpr int SEARCH = 15;

    // --------------------------------------------------------
    // Costas transition window
    // --------------------------------------------------------

    constexpr int WIN = 128;
    constexpr int HALF_WIN = WIN / 2;

    // --------------------------------------------------------
    // FIR
    // --------------------------------------------------------

    constexpr int FIR_TAPS = 129;
    constexpr double LPF_CUTOFF = 150.0;

    // --------------------------------------------------------
    // Fractional delays to test.
    //
    // The integer part remains fixed.
    // --------------------------------------------------------

    constexpr int NFRACTIONAL = 11;

    const double fractional_values[NFRACTIONAL] =
    {
        0.00,
        0.10,
        0.20,
        0.30,
        0.40,
        0.50,
        0.60,
        0.70,
        0.80,
        0.90,
        0.99
    };

    // --------------------------------------------------------
    // Generate test tones.
    // --------------------------------------------------------

    uint8_t tones[NSYM];

    make_test_tones(tones);

    // --------------------------------------------------------
    // Need one extra sample for fractional interpolation.
    // --------------------------------------------------------

    const size_t total_samples =
        (size_t)INTEGER_DELAY +
        (size_t)NSYM * SPS +
        2;

    // --------------------------------------------------------
    // Generate integer-aligned clean waveform.
    // --------------------------------------------------------

    std::vector<float> base(
        total_samples,
        0.0f);

    generate_ft8_continuous(
        base.data(),
        base.size(),
        tones,
        FREQ,
        0.0,
        AMPLITUDE,
        INTEGER_DELAY /
            (double)SAMPLE_RATE);

    // ========================================================
    // DESIGN FIR
    // ========================================================

    float fir[FIR_TAPS];

    const double fc =
        LPF_CUTOFF / SAMPLE_RATE;

    const int mid =
        FIR_TAPS / 2;

    double sum = 0.0;

    for (int k = 0; k < FIR_TAPS; ++k)
    {
        const int n =
            k - mid;

        double h;

        if (n == 0)
        {
            h = 2.0 * fc;
        }
        else
        {
            h =
                sin(2.0 * M_PI * fc * n) /
                (M_PI * n);
        }

        const double window =
            0.5 *
            (1.0 -
             cos(2.0 * M_PI * k /
                 (FIR_TAPS - 1)));

        h *= window;

        fir[k] = (float)h;

        sum += h;
    }

    for (int k = 0; k < FIR_TAPS; ++k)
        fir[k] /= (float)sum;

    // ========================================================
    // Build Costas transition list.
    //
    // The transitions are at the INTEGER reference delay.
    // The received waveform is subsequently shifted by the
    // fractional amount.
    // ========================================================

    struct Transition
    {
        int sample;
        int tone_before;
        int tone_after;
    };

    Transition transitions[22];

    int ntrans = 0;

    for (int k = 0; k < NSYM - 1; ++k)
    {
        if (tones[k] == tones[k + 1])
            continue;

        const int boundary =
            INTEGER_DELAY +
            (k + 1) * SPS;

        transitions[ntrans++] =
        {
            boundary,
            tones[k],
            tones[k + 1]
        };
    }

    // ========================================================
    // Downconvert + FIR function.
    //
    // We use a lambda so the exact same processing is applied
    // to the reference and received signal.
    // ========================================================

    auto downconvert_and_filter =
        [&](const std::vector<float>& input,
            std::vector<float>& out_i,
            std::vector<float>& out_q)
    {
        const size_t nsize =
            input.size();

        std::vector<float> i(
            nsize);

        std::vector<float> q(
            nsize);

        const double w =
            2.0 * M_PI * FREQ /
            SAMPLE_RATE;

        const float cw =
            cos(w);

        const float sw =
            sin(w);

        float c = 1.0f;
        float s = 0.0f;

        for (size_t n = 0;
             n < nsize;
             ++n)
        {
            const float x =
                input[n];

            i[n] =
                x * c;

            q[n] =
                -x * s;

            const float nc =
                c * cw -
                s * sw;

            const float ns =
                s * cw +
                c * sw;

            c = nc;
            s = ns;
        }

        out_i.resize(nsize);
        out_q.resize(nsize);

        for (size_t n = 0;
             n < nsize;
             ++n)
        {
            double ai = 0.0;
            double aq = 0.0;

            const int first =
                (n >= FIR_TAPS - 1)
                    ? (int)n -
                      FIR_TAPS + 1
                    : 0;

            for (int k = first;
                 k <= (int)n;
                 ++k)
            {
                const int h =
                    (int)n - k;

                ai +=
                    i[k] * fir[h];

                aq +=
                    q[k] * fir[h];
            }

            out_i[n] =
                (float)ai;

            out_q[n] =
                (float)aq;
        }
    };

    // --------------------------------------------------------
    // Reference baseband.
    // --------------------------------------------------------

    std::vector<float> ref_i;
    std::vector<float> ref_q;

    downconvert_and_filter(
        base,
        ref_i,
        ref_q);

    // --------------------------------------------------------
    // Extract reference transition windows.
    // --------------------------------------------------------

    float reference_i[22][WIN];
    float reference_q[22][WIN];

    for (int t = 0;
         t < ntrans;
         ++t)
    {
        const int center =
            transitions[t].sample;

        const int start =
            center - HALF_WIN;

        for (int m = 0;
             m < WIN;
             ++m)
        {
            const int idx =
                start + m;

            reference_i[t][m] =
                ref_i[idx];

            reference_q[t][m] =
                ref_q[idx];
        }
    }

    // ========================================================
    // Print header.
    // ========================================================

    printf("\n");

    printf("Sample rate        : %d Hz\n",
           SAMPLE_RATE);

    printf("Carrier            : %.3f Hz\n",
           FREQ);

    printf("FIR cutoff         : %.1f Hz\n",
           LPF_CUTOFF);

    printf("FIR taps            : %d\n",
           FIR_TAPS);

    printf("Window             : %d samples\n",
           WIN);

    printf("Search             : +/- %d samples\n",
           SEARCH);

    printf("Costas transitions : %d\n",
           ntrans);

    printf("\n");

    printf("%10s %10s %12s %12s %12s %12s\n",
           "True samp",
           "Frac",
           "Best int",
           "Parabola",
           "Estimated",
           "Error");

    printf("---------------------------------------------------------------------\n");

    // ========================================================
    // FRACTIONAL DELAY SWEEP
    // ========================================================

    for (int fi = 0;
         fi < NFRACTIONAL;
         ++fi)
    {
        const double fraction =
            fractional_values[fi];

        // ----------------------------------------------------
        // Apply fractional delay to the integer-aligned
        // waveform.
        // ----------------------------------------------------

        std::vector<float> received(
            total_samples,
            0.0f);

        apply_fractional_delay(
            base.data(),
            received.data(),
            received.size(),
            fraction);

        // ----------------------------------------------------
        // Downconvert + filter received signal.
        // ----------------------------------------------------

        std::vector<float> received_i;
        std::vector<float> received_q;

        downconvert_and_filter(
            received,
            received_i,
            received_q);

        // ====================================================
        // INTEGER SAMPLE SEARCH
        // ====================================================

        double scores[
            2 * SEARCH + 1
        ];

        double best_score =
            -1.0;

        int best_offset =
            0;

        for (int d = -SEARCH;
             d <= SEARCH;
             ++d)
        {
            double score =
                0.0;

            for (int t = 0;
                 t < ntrans;
                 ++t)
            {
                const int center =
                    transitions[t].sample +
                    d;

                const int start =
                    center -
                    HALF_WIN;

                double ci = 0.0;
                double cq = 0.0;

                double energy_received =
                    0.0;

                double energy_reference =
                    0.0;

                for (int m = 0;
                     m < WIN;
                     ++m)
                {
                    const int idx =
                        start + m;

                    const float bi =
                        received_i[idx];

                    const float bq =
                        received_q[idx];

                    const float ri =
                        reference_i[t][m];

                    const float rq =
                        reference_q[t][m];

                    // received * conj(reference)

                    ci +=
                        bi * ri +
                        bq * rq;

                    cq +=
                        bq * ri -
                        bi * rq;

                    energy_received +=
                        bi * bi +
                        bq * bq;

                    energy_reference +=
                        ri * ri +
                        rq * rq;
                }

                const double correlation_power =
                    ci * ci +
                    cq * cq;

                const double denominator =
                    energy_received *
                    energy_reference;

                if (denominator > 0.0)
                {
                    score +=
                        correlation_power /
                        denominator;
                }
            }

            scores[d + SEARCH] =
                score;

            if (score > best_score)
            {
                best_score =
                    score;

                best_offset =
                    d;
            }
        }

        // ====================================================
        // PARABOLIC INTERPOLATION
        // ====================================================

        double fractional_offset =
            0.0;

        if (best_offset > -SEARCH &&
            best_offset < SEARCH)
        {
            const double ym =
                scores[
                    best_offset -
                    1 +
                    SEARCH
                ];

            const double y0 =
                scores[
                    best_offset +
                    SEARCH
                ];

            const double yp =
                scores[
                    best_offset +
                    1 +
                    SEARCH
                ];

            const double denominator =
                ym -
                2.0 * y0 +
                yp;

            if (std::abs(denominator) >
                1e-15)
            {
                fractional_offset =
                    0.5 *
                    (ym - yp) /
                    denominator;
            }
        }

        // ----------------------------------------------------
        // Estimated delay.
        // ----------------------------------------------------

        const double estimated_samples =
            INTEGER_DELAY +
            best_offset +
            fractional_offset;

        const double true_samples =
            INTEGER_DELAY +
            fraction;

        const double error =
            estimated_samples -
            true_samples;

        printf("%10d %10.2f %12d %12.6f %12.6f %+12.6f\n",
               INTEGER_DELAY,
               fraction,
               best_offset,
               fractional_offset,
               estimated_samples,
               error);
    }

    printf("---------------------------------------------------------------------\n");

    printf("One sample       : %.9f ms\n",
           1000.0 / SAMPLE_RATE);

    printf("One tenth sample : %.9f ms\n",
           100.0 / SAMPLE_RATE);

    printf("\n");
    printf("D2 completed.\n");
    printf("============================================================\n");
}


// ============================================================
// D3 - FT8 DELAY ESTIMATOR: INTEGER vs PARABOLIC REFINEMENT
//
// Monte Carlo test with:
//   - random fractional delay [0, 1) sample
//   - AWGN
//   - integer correlation search +/- SEARCH samples
//   - 3-point parabolic interpolation
//
// Compares:
//   1) integer-grid estimate
//   2) parabolic sub-sample estimate
//
// The estimator itself is not modified.
// ============================================================

static void test_D3()
{
    printf("\n");
    printf("============================================================\n");
    printf("D3 - FT8 DELAY ESTIMATOR: INTEGER vs PARABOLIC\n");
    printf("============================================================\n");

    constexpr int SAMPLE_RATE = 12000;
    constexpr int NSYM        = 79;
    constexpr int SPS         = 1920;

    constexpr double FREQ      = 1001.0;
    constexpr double AMPLITUDE = 1000.0;

    constexpr int INTEGER_DELAY = 24000;
    constexpr int SEARCH        = 15;
    constexpr int WIN           = 128;

    constexpr int FIR_TAPS = 129;
    constexpr double LPF_CUTOFF = 150.0;

    constexpr int NTRIALS = 200;

    const double snr_db_list[] =
    {
        20.0,
        10.0,
        0.0,
        -10.0
    };

    // --------------------------------------------------------
    // FT8 tones
    // --------------------------------------------------------

    uint8_t tones[79];
    make_test_tones(tones);

    // --------------------------------------------------------
    // Generate integer-aligned reference waveform
    // --------------------------------------------------------

    const size_t nsamples =
        INTEGER_DELAY + NSYM * SPS + SAMPLE_RATE;

    std::vector<float> base(nsamples, 0.0f);

    generate_ft8_continuous(
        base.data(),
        base.size(),
        tones,
        FREQ,
        0.0,
        AMPLITUDE,
        (double)INTEGER_DELAY / SAMPLE_RATE);

    // --------------------------------------------------------
    // FIR low-pass filter
    // Same filter as D1 / D2
    // --------------------------------------------------------

    std::vector<float> fir_coeff(FIR_TAPS);

    {
        const int M = FIR_TAPS - 1;
        const double fc =
            LPF_CUTOFF / SAMPLE_RATE;

        for (int n = 0; n < FIR_TAPS; ++n)
        {
            const double x = n - M / 2.0;

            double h;

            if (std::abs(x) < 1e-12)
            {
                h = 2.0 * fc;
            }
            else
            {
                h = std::sin(2.0 * PI * fc * x)
                    / (PI * x);
            }

            // Hann window
            const double w =
                0.5 *
                (1.0 -
                 std::cos(2.0 * PI * n / M));

            fir_coeff[n] =
                (float)(h * w);
        }

        // Normalize DC gain
        double sum = 0.0;

        for (float h : fir_coeff)
            sum += h;

        for (float &h : fir_coeff)
            h = (float)(h / sum);
    }

    // --------------------------------------------------------
    // RF -> complex baseband -> FIR
    // --------------------------------------------------------

    auto downconvert =
        [&](const std::vector<float> &rf,
            std::vector<std::complex<float>> &bb)
    {
        bb.assign(rf.size(),
                  std::complex<float>(0.0f, 0.0f));

        const double w =
            2.0 * PI * FREQ / SAMPLE_RATE;

        std::vector<float> mix_i(rf.size());
        std::vector<float> mix_q(rf.size());

        for (size_t n = 0; n < rf.size(); ++n)
        {
            const double phase = w * n;

            mix_i[n] =
                (float)(rf[n] * std::cos(phase));

            mix_q[n] =
                (float)(-rf[n] * std::sin(phase));
        }

        // FIR
        for (size_t n = 0; n < rf.size(); ++n)
        {
            double si = 0.0;
            double sq = 0.0;

            const size_t kmax =
                std::min<size_t>(
                    FIR_TAPS,
                    n + 1);

            for (size_t k = 0; k < kmax; ++k)
            {
                si +=
                    fir_coeff[k] *
                    mix_i[n - k];

                sq +=
                    fir_coeff[k] *
                    mix_q[n - k];
            }

            bb[n] =
                std::complex<float>(
                    (float)si,
                    (float)sq);
        }
    };

    // --------------------------------------------------------
    // Costas transitions
    //
    // IMPORTANT:
    // exactly the same definition as D1/D2.
    //
    // A transition exists whenever:
    //
    //     tones[k] != tones[k+1]
    //
    // There are 22 transitions.
    // --------------------------------------------------------

    struct Transition
    {
        int sample;
        int tone_before;
        int tone_after;
    };

    Transition transitions[22];

    int ntrans = 0;

    for (int k = 0; k < NSYM - 1; ++k)
    {
        if (tones[k] == tones[k + 1])
            continue;

        const int boundary =
            INTEGER_DELAY + (k + 1) * SPS;

        transitions[ntrans++] =
        {
            boundary,
            tones[k],
            tones[k + 1]
        };
    }

    assert(ntrans == 22);

    printf("Costas transitions : %d\n", ntrans);
    printf("Window             : %d samples\n", WIN);
    printf("Search             : +/- %d samples\n",
           SEARCH);
    printf("Trials             : %d\n", NTRIALS);

    // --------------------------------------------------------
    // Reference transition waveforms
    //
    // The reference is extracted from the same baseband
    // waveform used by the estimator.
    // --------------------------------------------------------

    const int fir_delay =
        FIR_TAPS / 2;

    std::vector<std::complex<float>> ref(
        ntrans * WIN);

    {
        std::vector<std::complex<float>> bb;

        downconvert(base, bb);

        for (int t = 0; t < ntrans; ++t)
        {
            const int center =
                transitions[t].sample +
                fir_delay;

            const int start =
                center - WIN / 2;

            for (int i = 0; i < WIN; ++i)
            {
                const int idx = start + i;

                if (idx >= 0 &&
                    idx < (int)bb.size())
                {
                    ref[t * WIN + i] =
                        bb[idx];
                }
                else
                {
                    ref[t * WIN + i] =
                        std::complex<float>(0.0f, 0.0f);
                }
            }
        }
    }

    // --------------------------------------------------------
    // RNG
    // --------------------------------------------------------

    std::mt19937 rng(1234567);

    std::uniform_real_distribution<double>
        frac_dist(0.0, 1.0);

    // --------------------------------------------------------
    // Results
    // --------------------------------------------------------

    for (double snr_db : snr_db_list)
    {
        double sum_sq_integer = 0.0;
        double sum_sq_parabola = 0.0;

        double max_err_integer = 0.0;
        double max_err_parabola = 0.0;

        int wrong_integer = 0;

        // ----------------------------------------------------
        // Signal power
        //
        // IMPORTANT:
        // calculate only over the active FT8 interval,
        // not over the whole 15-second buffer.
        // ----------------------------------------------------

        const double signal_power =
            calculate_signal_power(
                base.data(),
                base.size(),
                (double)INTEGER_DELAY / SAMPLE_RATE);

        const double noise_power =
            signal_power /
            std::pow(10.0, snr_db / 10.0);

        const double noise_sigma =
            std::sqrt(noise_power);

        std::normal_distribution<double>
            noise_dist(0.0, noise_sigma);

        // ----------------------------------------------------
        // Monte Carlo
        // ----------------------------------------------------

        for (int trial = 0;
             trial < NTRIALS;
             ++trial)
        {
            // Random fractional delay
            const double true_frac =
                frac_dist(rng);

            const double true_delay =
                (double)INTEGER_DELAY +
                true_frac;

            // ------------------------------------------------
            // Create fractionally delayed signal
            //
            // Linear interpolation is used ONLY to create
            // the test stimulus. It is NOT part of the
            // estimator.
            // ------------------------------------------------

            std::vector<float> delayed(nsamples);

            apply_fractional_delay(
                base.data(),
                delayed.data(),
                base.size(),
                true_frac);

            // ------------------------------------------------
            // Add AWGN
            // ------------------------------------------------

            for (float &x : delayed)
                x += (float)noise_dist(rng);

            // ------------------------------------------------
            // Downconvert + LPF
            // ------------------------------------------------

            std::vector<std::complex<float>> bb;

            downconvert(delayed, bb);

            // ------------------------------------------------
            // Integer correlation search
            //
            // Normalized correlation:
            //
            // rho² =
            // |sum rx * conj(ref)|²
            // -----------------------------------------------
            // (sum |rx|²) (sum |ref|²)
            //
            // Score = sum over all transitions.
            // ------------------------------------------------

            std::vector<double> score(
                2 * SEARCH + 1);

            for (int offset = -SEARCH;
                 offset <= SEARCH;
                 ++offset)
            {
                double total_score = 0.0;

                for (int t = 0;
                     t < ntrans;
                     ++t)
                {
                    const int center =
                        transitions[t].sample +
                        fir_delay +
                        offset;

                    const int start =
                        center - WIN / 2;

                    std::complex<double> corr(0.0, 0.0);

                    double rx_power = 0.0;
                    double ref_power = 0.0;

                    for (int i = 0;
                         i < WIN;
                         ++i)
                    {
                        const int idx =
                            start + i;

                        if (idx < 0 ||
                            idx >= (int)bb.size())
                            continue;

                        const auto rx =
                            bb[idx];

                        const auto rr =
                            ref[t * WIN + i];

                        corr +=
                            (std::complex<double>)rx *
                            std::conj(
                                (std::complex<double>)rr);

                        rx_power +=
                            std::norm(
                                (std::complex<double>)rx);

                        ref_power +=
                            std::norm(
                                (std::complex<double>)rr);
                    }

                    if (rx_power > 0.0 &&
                        ref_power > 0.0)
                    {
                        const double c2 =
                            std::norm(corr);

                        total_score +=
                            c2 /
                            (rx_power * ref_power);
                    }
                }

                score[offset + SEARCH] =
                    total_score;
            }

            // ------------------------------------------------
            // Find integer maximum
            // ------------------------------------------------

            int best_offset = 0;

            double best_score =
                score[SEARCH];

            for (int offset = -SEARCH;
                 offset <= SEARCH;
                 ++offset)
            {
                const double s =
                    score[offset + SEARCH];

                if (s > best_score)
                {
                    best_score = s;
                    best_offset = offset;
                }
            }

            const double integer_estimate =
                (double)INTEGER_DELAY +
                best_offset;

            const double integer_error =
                integer_estimate -
                true_delay;

            // ------------------------------------------------
            // Integer-grid error
            // ------------------------------------------------

            sum_sq_integer +=
                integer_error *
                integer_error;

            max_err_integer =
                std::max(
                    max_err_integer,
                    std::abs(integer_error));

            // Count estimates that are not the nearest
            // integer to the true fractional delay.
            const int expected_offset =
                (true_frac >= 0.5) ? 1 : 0;

            if (best_offset != expected_offset)
                ++wrong_integer;

            // ------------------------------------------------
            // Parabolic interpolation
            //
            // delta =
            // 0.5 * (S[-1] - S[+1]) /
            //       (S[-1] - 2*S[0] + S[+1])
            // ------------------------------------------------

            double parabola = 0.0;

            if (best_offset > -SEARCH &&
                best_offset < SEARCH)
            {
                const double sm =
                    score[best_offset - 1 + SEARCH];

                const double s0 =
                    score[best_offset + SEARCH];

                const double sp =
                    score[best_offset + 1 + SEARCH];

                const double denominator =
                    sm - 2.0 * s0 + sp;

                if (std::abs(denominator) > 1e-20)
                {
                    parabola =
                        0.5 *
                        (sm - sp) /
                        denominator;
                }
            }

            const double parabolic_estimate =
                (double)INTEGER_DELAY +
                best_offset +
                parabola;

            const double parabolic_error =
                parabolic_estimate -
                true_delay;

            // ------------------------------------------------
            // Parabolic error
            // ------------------------------------------------

            sum_sq_parabola +=
                parabolic_error *
                parabolic_error;

            max_err_parabola =
                std::max(
                    max_err_parabola,
                    std::abs(parabolic_error));
        }

        // ----------------------------------------------------
        // Statistics
        // ----------------------------------------------------

        const double rms_integer =
            std::sqrt(
                sum_sq_integer /
                NTRIALS);

        const double rms_parabola =
            std::sqrt(
                sum_sq_parabola /
                NTRIALS);

        // ----------------------------------------------------
        // Results
        // ----------------------------------------------------

        printf("\n");
        printf("SNR = %+5.1f dB\n", snr_db);

        printf("  Integer:\n");
        printf("    RMS error       = %10.6f samples\n",
               rms_integer);
        printf("    RMS error       = %10.6f ms\n",
               rms_integer * 1000.0 / SAMPLE_RATE);
        printf("    Max error       = %10.6f samples\n",
               max_err_integer);
        printf("    Wrong integer   = %d / %d\n",
               wrong_integer,
               NTRIALS);

        printf("  Parabolic:\n");
        printf("    RMS error       = %10.6f samples\n",
               rms_parabola);
        printf("    RMS error       = %10.6f ms\n",
               rms_parabola * 1000.0 / SAMPLE_RATE);
        printf("    Max error       = %10.6f samples\n",
               max_err_parabola);

        if (rms_parabola > 0.0)
        {
            printf("  Improvement      = %8.2fx\n",
                   rms_integer / rms_parabola);
        }
    }

    printf("\n");
}

// ============================================================
// D3a - FT8 DELAY ESTIMATOR: NOISELESS FRACTIONAL DELAY
//
// Monte Carlo test with:
//   - random fractional delay [0, 1) sample
//   - NO AWGN
//   - integer correlation search +/- SEARCH samples
//   - 3-point parabolic interpolation
//
// Compares:
//   1) integer-grid estimate
//   2) parabolic sub-sample estimate
//
// Purpose:
//   Measure the intrinsic accuracy of the estimator,
//   independently from noise.
// ============================================================

static void test_D3a()
{
    printf("\n");
    printf("============================================================\n");
    printf("D3a - FT8 DELAY ESTIMATOR: NOISELESS FRACTIONAL DELAY\n");
    printf("============================================================\n");

    constexpr int SAMPLE_RATE = 12000;
    constexpr int NSYM        = 79;
    constexpr int SPS         = 1920;

    constexpr double FREQ      = 1001.0;
    constexpr double AMPLITUDE = 1000.0;

    constexpr int INTEGER_DELAY = 24000;
    constexpr int SEARCH        = 15;
    constexpr int WIN           = 128;

    constexpr int FIR_TAPS = 129;
    constexpr double LPF_CUTOFF = 150.0;

    constexpr int NTRIALS = 200;

    // --------------------------------------------------------
    // FT8 tones
    // --------------------------------------------------------

    uint8_t tones[79];
    make_test_tones(tones);

    // --------------------------------------------------------
    // Generate integer-aligned reference waveform
    // --------------------------------------------------------

    const size_t nsamples =
        INTEGER_DELAY + NSYM * SPS + SAMPLE_RATE;

    std::vector<float> base(nsamples, 0.0f);

    generate_ft8_continuous(
        base.data(),
        base.size(),
        tones,
        FREQ,
        0.0,
        AMPLITUDE,
        (double)INTEGER_DELAY / SAMPLE_RATE);

    // --------------------------------------------------------
    // FIR
    // Same filter used by D1 / D2
    // --------------------------------------------------------

    std::vector<float> fir_coeff(FIR_TAPS);

    {
        const int M = FIR_TAPS - 1;

        const double fc =
            LPF_CUTOFF / SAMPLE_RATE;

        for (int n = 0; n < FIR_TAPS; ++n)
        {
            const double x =
                n - M / 2.0;

            double h;

            if (std::abs(x) < 1e-12)
            {
                h = 2.0 * fc;
            }
            else
            {
                h =
                    std::sin(2.0 * PI * fc * x)
                    / (PI * x);
            }

            const double w =
                0.5 *
                (1.0 -
                 std::cos(2.0 * PI * n / M));

            fir_coeff[n] =
                (float)(h * w);
        }

        // Normalize DC gain
        double sum = 0.0;

        for (float h : fir_coeff)
            sum += h;

        for (float &h : fir_coeff)
            h = (float)(h / sum);
    }

    // --------------------------------------------------------
    // RF -> complex baseband -> FIR
    // --------------------------------------------------------

    auto downconvert =
        [&](const std::vector<float> &rf,
            std::vector<std::complex<float>> &bb)
    {
        bb.assign(
            rf.size(),
            std::complex<float>(0.0f, 0.0f));

        const double w =
            2.0 * PI * FREQ / SAMPLE_RATE;

        std::vector<float> mix_i(rf.size());
        std::vector<float> mix_q(rf.size());

        for (size_t n = 0;
             n < rf.size();
             ++n)
        {
            const double phase =
                w * n;

            mix_i[n] =
                (float)(
                    rf[n] *
                    std::cos(phase));

            mix_q[n] =
                (float)(
                    -rf[n] *
                    std::sin(phase));
        }

        for (size_t n = 0;
             n < rf.size();
             ++n)
        {
            double si = 0.0;
            double sq = 0.0;

            const size_t kmax =
                std::min<size_t>(
                    FIR_TAPS,
                    n + 1);

            for (size_t k = 0;
                 k < kmax;
                 ++k)
            {
                si +=
                    fir_coeff[k] *
                    mix_i[n - k];

                sq +=
                    fir_coeff[k] *
                    mix_q[n - k];
            }

            bb[n] =
                std::complex<float>(
                    (float)si,
                    (float)sq);
        }
    };

    // --------------------------------------------------------
    // Costas transitions
    //
    // Exactly the same definition as D1 / D2.
    // --------------------------------------------------------

    struct Transition
    {
        int sample;
        int tone_before;
        int tone_after;
    };

    Transition transitions[22];

    int ntrans = 0;

    for (int k = 0;
         k < NSYM - 1;
         ++k)
    {
        if (tones[k] == tones[k + 1])
            continue;

        const int boundary =
            INTEGER_DELAY +
            (k + 1) * SPS;

        transitions[ntrans++] =
        {
            boundary,
            tones[k],
            tones[k + 1]
        };
    }

    assert(ntrans == 22);

    printf("Sample rate        : %d Hz\n",
           SAMPLE_RATE);
    printf("Carrier            : %.3f Hz\n",
           FREQ);
    printf("LPF cutoff         : %.1f Hz\n",
           LPF_CUTOFF);
    printf("FIR taps           : %d\n",
           FIR_TAPS);
    printf("FIR group delay    : %d samples\n",
           FIR_TAPS / 2);
    printf("\n");

    printf("Costas transitions : %d\n",
           ntrans);
    printf("Window             : %d samples\n",
           WIN);
    printf("Search             : +/- %d samples\n",
           SEARCH);
    printf("Trials             : %d\n",
           NTRIALS);

    // --------------------------------------------------------
    // Reference baseband
    // --------------------------------------------------------

    std::vector<std::complex<float>> ref(
        ntrans * WIN);

    {
        std::vector<std::complex<float>> bb;

        downconvert(base, bb);

        const int fir_delay =
            FIR_TAPS / 2;

        for (int t = 0;
             t < ntrans;
             ++t)
        {
            const int center =
                transitions[t].sample +
                fir_delay;

            const int start =
                center - WIN / 2;

            for (int i = 0;
                 i < WIN;
                 ++i)
            {
                const int idx =
                    start + i;

                if (idx >= 0 &&
                    idx < (int)bb.size())
                {
                    ref[t * WIN + i] =
                        bb[idx];
                }
                else
                {
                    ref[t * WIN + i] =
                        std::complex<float>(
                            0.0f,
                            0.0f);
                }
            }
        }
    }

    // --------------------------------------------------------
    // RNG
    // --------------------------------------------------------

    std::mt19937 rng(1234567);

    std::uniform_real_distribution<double>
        frac_dist(0.0, 1.0);

    // --------------------------------------------------------
    // Statistics
    // --------------------------------------------------------

    double sum_error_integer = 0.0;
    double sum_error_parabola = 0.0;

    double sum_sq_integer = 0.0;
    double sum_sq_parabola = 0.0;

    double max_abs_integer = 0.0;
    double max_abs_parabola = 0.0;

    int wrong_integer = 0;

    // --------------------------------------------------------
    // Monte Carlo
    // --------------------------------------------------------

    for (int trial = 0;
         trial < NTRIALS;
         ++trial)
    {
        // ----------------------------------------------------
        // Random fractional delay
        // ----------------------------------------------------

        const double true_frac =
            frac_dist(rng);

        const double true_delay =
            INTEGER_DELAY +
            true_frac;

        // ----------------------------------------------------
        // Fractionally delayed signal
        //
        // IMPORTANT:
        // linear interpolation is used ONLY to create
        // the test stimulus.
        // ----------------------------------------------------

        std::vector<float> delayed(nsamples);

        apply_fractional_delay(
            base.data(),
            delayed.data(),
            base.size(),
            true_frac);

        // ----------------------------------------------------
        // Downconvert + LPF
        // ----------------------------------------------------

        std::vector<std::complex<float>> bb;

        downconvert(delayed, bb);

        // ----------------------------------------------------
        // Correlation score for integer offsets
        // ----------------------------------------------------

        std::vector<double> score(
            2 * SEARCH + 1);

        const int fir_delay =
            FIR_TAPS / 2;

        for (int offset = -SEARCH;
             offset <= SEARCH;
             ++offset)
        {
            double total_score = 0.0;

            for (int t = 0;
                 t < ntrans;
                 ++t)
            {
                const int center =
                    transitions[t].sample +
                    fir_delay +
                    offset;

                const int start =
                    center - WIN / 2;

                std::complex<double> corr(
                    0.0,
                    0.0);

                double rx_power = 0.0;
                double ref_power = 0.0;

                for (int i = 0;
                     i < WIN;
                     ++i)
                {
                    const int idx =
                        start + i;

                    if (idx < 0 ||
                        idx >= (int)bb.size())
                        continue;

                    const std::complex<float> rx =
                        bb[idx];

                    const std::complex<float> rr =
                        ref[t * WIN + i];

                    corr +=
                        std::complex<double>(rx.real(),
                                             rx.imag()) *
                        std::conj(
                            std::complex<double>(
                                rr.real(),
                                rr.imag()));

                    rx_power +=
                        (double)rx.real() *
                        rx.real() +
                        (double)rx.imag() *
                        rx.imag();

                    ref_power +=
                        (double)rr.real() *
                        rr.real() +
                        (double)rr.imag() *
                        rr.imag();
                }

                if (rx_power > 0.0 &&
                    ref_power > 0.0)
                {
                    const double c2 =
                        corr.real() * corr.real() +
                        corr.imag() * corr.imag();

                    total_score +=
                        c2 /
                        (rx_power * ref_power);
                }
            }

            score[offset + SEARCH] =
                total_score;
        }

        // ----------------------------------------------------
        // Integer maximum
        // ----------------------------------------------------

        int best_offset = 0;

        double best_score =
            score[SEARCH];

        for (int offset = -SEARCH;
             offset <= SEARCH;
             ++offset)
        {
            const double s =
                score[offset + SEARCH];

            if (s > best_score)
            {
                best_score = s;
                best_offset = offset;
            }
        }

        const double integer_estimate =
            INTEGER_DELAY +
            best_offset;

        const double integer_error =
            integer_estimate -
            true_delay;

        // ----------------------------------------------------
        // Integer statistics
        // ----------------------------------------------------

        sum_error_integer +=
            integer_error;

        sum_sq_integer +=
            integer_error *
            integer_error;

        max_abs_integer =
            std::max(
                max_abs_integer,
                std::abs(integer_error));

        // The ideal integer result is the nearest sample.
        const int expected_offset =
            (true_frac >= 0.5)
                ? 1
                : 0;

        if (best_offset != expected_offset)
            ++wrong_integer;

        // ----------------------------------------------------
        // Parabolic interpolation
        // ----------------------------------------------------

        double parabola = 0.0;

        if (best_offset > -SEARCH &&
            best_offset < SEARCH)
        {
            const double sm =
                score[best_offset - 1 + SEARCH];

            const double s0 =
                score[best_offset + SEARCH];

            const double sp =
                score[best_offset + 1 + SEARCH];

            const double denominator =
                sm - 2.0 * s0 + sp;

            if (std::abs(denominator) > 1e-20)
            {
                parabola =
                    0.5 *
                    (sm - sp) /
                    denominator;
            }
        }

        const double parabolic_estimate =
            INTEGER_DELAY +
            best_offset +
            parabola;

        const double parabolic_error =
            parabolic_estimate -
            true_delay;

        // ----------------------------------------------------
        // Parabolic statistics
        // ----------------------------------------------------

        sum_error_parabola +=
            parabolic_error;

        sum_sq_parabola +=
            parabolic_error *
            parabolic_error;

        max_abs_parabola =
            std::max(
                max_abs_parabola,
                std::abs(parabolic_error));
    }

    // --------------------------------------------------------
    // Final statistics
    // --------------------------------------------------------

    const double mean_integer =
        sum_error_integer /
        NTRIALS;

    const double mean_parabola =
        sum_error_parabola /
        NTRIALS;

    const double rms_integer =
        std::sqrt(
            sum_sq_integer /
            NTRIALS);

    const double rms_parabola =
        std::sqrt(
            sum_sq_parabola /
            NTRIALS);

    // --------------------------------------------------------
    // Results
    // --------------------------------------------------------

    printf("\n");
    printf("------------------------------------------------------------\n");
    printf("Results\n");
    printf("------------------------------------------------------------\n");

    printf("\nInteger-grid estimator:\n");

    printf("  Mean error       = %10.6f samples\n",
           mean_integer);

    printf("  RMS error        = %10.6f samples\n",
           rms_integer);

    printf("  RMS error        = %10.6f ms\n",
           rms_integer * 1000.0 / SAMPLE_RATE);

    printf("  Max abs error    = %10.6f samples\n",
           max_abs_integer);

    printf("  Wrong integer    = %d / %d\n",
           wrong_integer,
           NTRIALS);

    printf("\nParabolic estimator:\n");

    printf("  Mean error       = %10.6f samples\n",
           mean_parabola);

    printf("  RMS error        = %10.6f samples\n",
           rms_parabola);

    printf("  RMS error        = %10.6f ms\n",
           rms_parabola * 1000.0 / SAMPLE_RATE);

    printf("  Max abs error    = %10.6f samples\n",
           max_abs_parabola);

    printf("\n");

    if (rms_parabola > 0.0)
    {
        printf("Improvement        = %10.3fx\n",
               rms_integer / rms_parabola);
    }

    printf("\n");
    printf("One sample         : %.9f ms\n",
           1000.0 / SAMPLE_RATE);

    printf("One tenth sample   : %.9f ms\n",
           100.0 / SAMPLE_RATE);

    printf("\nD3a completed.\n");
    printf("============================================================\n");
}

// ============================================================
// D3b - FT8 DELAY ESTIMATOR: NOISY FRACTIONAL DELAY
//
// Monte Carlo test:
//   - random fractional-sample delay [0,1)
//   - AWGN
//   - integer-grid delay estimator
//   - parabolic interpolation
//
// Same estimator as D3a.
// The fractional delay is used ONLY to create the test signal.
// ============================================================

// ============================================================
// D3b - FT8 DELAY ESTIMATOR: FRACTIONAL DELAY + AWGN
//
// Based EXACTLY on D3a.
//
// Differences from D3a:
//   - random fractional delay [0,1) sample
//   - AWGN added at several SNR levels
//   - integer correlation search +/- SEARCH samples
//   - 3-point parabolic interpolation
//
// Compares:
//   1) integer-grid estimate
//   2) parabolic sub-sample estimate
//
// Purpose:
//   Measure delay-estimator accuracy versus noise,
//   while keeping the estimator itself identical to D3a.
// ============================================================

static void test_D3b()
{
    printf("\n");
    printf("============================================================\n");
    printf("D3b - FT8 DELAY ESTIMATOR: FRACTIONAL DELAY + AWGN\n");
    printf("============================================================\n");

    constexpr int SAMPLE_RATE = 12000;
    constexpr int NSYM        = 79;
    constexpr int SPS         = 1920;

    constexpr double FREQ      = 1001.0;
    constexpr double AMPLITUDE = 1000.0;

    constexpr int INTEGER_DELAY = 24000;
    constexpr int SEARCH        = 15;
    constexpr int WIN           = 128;

    constexpr int FIR_TAPS = 129;
    constexpr double LPF_CUTOFF = 150.0;

    constexpr int NTRIALS = 200;

    // SNR test points
    const double snr_db_list[] =
    {
        20.0,
        15.0,
        10.0,
         5.0,
         0.0,
        -5.0,
        -10.0
    };

    constexpr int NSNR =
        sizeof(snr_db_list) / sizeof(snr_db_list[0]);

    // --------------------------------------------------------
    // FT8 tones
    // --------------------------------------------------------

    uint8_t tones[79];
    make_test_tones(tones);

    // --------------------------------------------------------
    // Generate integer-aligned reference waveform
    // --------------------------------------------------------

    const size_t nsamples =
        INTEGER_DELAY + NSYM * SPS + SAMPLE_RATE;

    std::vector<float> base(nsamples, 0.0f);

    generate_ft8_continuous(
        base.data(),
        base.size(),
        tones,
        FREQ,
        0.0,
        AMPLITUDE,
        (double)INTEGER_DELAY / SAMPLE_RATE);

    // --------------------------------------------------------
    // FIR
    // --------------------------------------------------------

    std::vector<float> fir_coeff(FIR_TAPS);

    {
        const int M = FIR_TAPS - 1;
        const double fc =
            LPF_CUTOFF / SAMPLE_RATE;

        for (int n = 0; n < FIR_TAPS; ++n)
        {
            const double x =
                n - M / 2.0;

            double h;

            if (std::abs(x) < 1e-12)
            {
                h = 2.0 * fc;
            }
            else
            {
                h =
                    std::sin(2.0 * PI * fc * x) /
                    (PI * x);
            }

            const double w =
                0.5 *
                (1.0 -
                 std::cos(2.0 * PI * n / M));

            fir_coeff[n] =
                (float)(h * w);
        }

        // Normalize DC gain
        double sum = 0.0;

        for (float h : fir_coeff)
            sum += h;

        for (float &h : fir_coeff)
            h = (float)(h / sum);
    }

    // --------------------------------------------------------
    // RF -> complex baseband -> FIR
    //
    // EXACTLY as in D3a
    // --------------------------------------------------------

    auto downconvert =
        [&](const std::vector<float> &rf,
            std::vector<std::complex<float>> &bb)
    {
        bb.assign(
            rf.size(),
            std::complex<float>(0.0f, 0.0f));

        const double w =
            2.0 * PI * FREQ / SAMPLE_RATE;

        std::vector<float> mix_i(rf.size());
        std::vector<float> mix_q(rf.size());

        for (size_t n = 0; n < rf.size(); ++n)
        {
            const double phase =
                w * n;

            mix_i[n] =
                (float)(rf[n] *
                        std::cos(phase));

            mix_q[n] =
                (float)(-rf[n] *
                        std::sin(phase));
        }

        for (size_t n = 0; n < rf.size(); ++n)
        {
            double si = 0.0;
            double sq = 0.0;

            const size_t kmax =
                std::min<size_t>(
                    FIR_TAPS,
                    n + 1);

            for (size_t k = 0;
                 k < kmax;
                 ++k)
            {
                si +=
                    fir_coeff[k] *
                    mix_i[n - k];

                sq +=
                    fir_coeff[k] *
                    mix_q[n - k];
            }

            bb[n] =
                std::complex<float>(
                    (float)si,
                    (float)sq);
        }
    };

    // --------------------------------------------------------
    // Costas transitions
    //
    // EXACTLY as in D3a
    // --------------------------------------------------------

    struct Transition
    {
        int sample;
        int tone_before;
        int tone_after;
    };

    Transition transitions[22];
    int ntrans = 0;

    for (int k = 0;
         k < NSYM - 1;
         ++k)
    {
        if (tones[k] == tones[k + 1])
            continue;

        const int boundary =
            INTEGER_DELAY +
            (k + 1) * SPS;

        transitions[ntrans++] =
        {
            boundary,
            tones[k],
            tones[k + 1]
        };
    }

    assert(ntrans == 22);

    // --------------------------------------------------------
    // Print configuration
    // --------------------------------------------------------

    printf("Sample rate        : %d Hz\n",
           SAMPLE_RATE);

    printf("Carrier            : %.3f Hz\n",
           FREQ);

    printf("LPF cutoff         : %.1f Hz\n",
           LPF_CUTOFF);

    printf("FIR taps           : %d\n",
           FIR_TAPS);

    printf("FIR group delay    : %d samples\n",
           FIR_TAPS / 2);

    printf("\n");

    printf("Costas transitions : %d\n",
           ntrans);

    printf("Window             : %d samples\n",
           WIN);

    printf("Search             : +/- %d samples\n",
           SEARCH);

    printf("Trials             : %d\n",
           NTRIALS);

    // --------------------------------------------------------
    // Reference baseband
    //
    // EXACTLY as in D3a
    // --------------------------------------------------------

    std::vector<std::complex<float>>
        ref(ntrans * WIN);

    {
        std::vector<std::complex<float>> bb;

        downconvert(base, bb);

        const int fir_delay =
            FIR_TAPS / 2;

        for (int t = 0;
             t < ntrans;
             ++t)
        {
            const int center =
                transitions[t].sample +
                fir_delay;

            const int start =
                center - WIN / 2;

            for (int i = 0;
                 i < WIN;
                 ++i)
            {
                const int idx =
                    start + i;

                if (idx >= 0 &&
                    idx < (int)bb.size())
                {
                    ref[t * WIN + i] =
                        bb[idx];
                }
                else
                {
                    ref[t * WIN + i] =
                        std::complex<float>(
                            0.0f,
                            0.0f);
                }
            }
        }
    }

    // --------------------------------------------------------
    // RNG
    // --------------------------------------------------------

    std::mt19937 rng(1234567);

    std::uniform_real_distribution<double>
        frac_dist(0.0, 1.0);

    // --------------------------------------------------------
    // SNR loop
    // --------------------------------------------------------

    for (int isnr = 0;
         isnr < NSNR;
         ++isnr)
    {
        const double snr_db =
            snr_db_list[isnr];

        // Statistics

        double sum_error_integer = 0.0;
        double sum_error_parabola = 0.0;

        double sum_sq_integer = 0.0;
        double sum_sq_parabola = 0.0;

        double max_abs_integer = 0.0;
        double max_abs_parabola = 0.0;

        int wrong_integer = 0;

        // ----------------------------------------------------
        // Monte Carlo
        // ----------------------------------------------------

        for (int trial = 0;
             trial < NTRIALS;
             ++trial)
        {
            const double true_frac =
                frac_dist(rng);

            const double true_delay =
                INTEGER_DELAY +
                true_frac;

            // ------------------------------------------------
            // Fractionally delayed signal
            //
            // EXACTLY as D3a
            // ------------------------------------------------

            std::vector<float> delayed(nsamples);

            apply_fractional_delay(
                base.data(),
                delayed.data(),
                base.size(),
                true_frac);

            // ------------------------------------------------
            // Calculate signal power ONLY over active FT8
            // interval.
            // ------------------------------------------------

            const double signal_power =
                calculate_signal_power(
                    delayed.data(),
                    delayed.size(),
                    true_delay /
                        SAMPLE_RATE);

            // ------------------------------------------------
            // AWGN sigma
            //
            // SNR = signal_power / noise_power
            // ------------------------------------------------

            const double noise_power =
                signal_power /
                std::pow(
                    10.0,
                    snr_db / 10.0);

            const double sigma =
                std::sqrt(noise_power);

            // ------------------------------------------------
            // Add AWGN
            // ------------------------------------------------

            std::vector<float> noisy =
                delayed;

            add_awgn(
                noisy.data(),
                noisy.size(),
                sigma,
                rng);

            // ------------------------------------------------
            // RF -> complex baseband
            //
            // IMPORTANT:
            // noise is present before downconversion,
            // exactly as it would be in the RF signal.
            // ------------------------------------------------

            std::vector<std::complex<float>> bb;

            downconvert(
                noisy,
                bb);

            // ------------------------------------------------
            // Correlation score
            //
            // EXACTLY as D3a
            // ------------------------------------------------

            std::vector<double> score(
                2 * SEARCH + 1);

            const int fir_delay =
                FIR_TAPS / 2;

            for (int offset = -SEARCH;
                 offset <= SEARCH;
                 ++offset)
            {
                double total_score = 0.0;

                for (int t = 0;
                     t < ntrans;
                     ++t)
                {
                    const int center =
                        transitions[t].sample +
                        fir_delay +
                        offset;

                    const int start =
                        center - WIN / 2;

                    std::complex<double> corr(
                        0.0,
                        0.0);

                    double rx_power = 0.0;
                    double ref_power = 0.0;

                    for (int i = 0;
                         i < WIN;
                         ++i)
                    {
                        const int idx =
                            start + i;

                        if (idx < 0 ||
                            idx >= (int)bb.size())
                            continue;

                        const std::complex<float> rx =
                            bb[idx];

                        const std::complex<float> rr =
                            ref[t * WIN + i];

                        corr +=
                            std::complex<double>(
                                rx.real(),
                                rx.imag()) *
                            std::conj(
                                std::complex<double>(
                                    rr.real(),
                                    rr.imag()));

                        rx_power +=
                            (double)rx.real() *
                            rx.real() +
                            (double)rx.imag() *
                            rx.imag();

                        ref_power +=
                            (double)rr.real() *
                            rr.real() +
                            (double)rr.imag() *
                            rr.imag();
                    }

                    if (rx_power > 0.0 &&
                        ref_power > 0.0)
                    {
                        const double c2 =
                            corr.real() *
                            corr.real() +
                            corr.imag() *
                            corr.imag();

                        total_score +=
                            c2 /
                            (rx_power *
                             ref_power);
                    }
                }

                score[offset + SEARCH] =
                    total_score;
            }

            // ------------------------------------------------
            // Integer maximum
            // ------------------------------------------------

            int best_offset = 0;

            double best_score =
                score[SEARCH];

            for (int offset = -SEARCH;
                 offset <= SEARCH;
                 ++offset)
            {
                const double s =
                    score[offset + SEARCH];

                if (s > best_score)
                {
                    best_score = s;
                    best_offset = offset;
                }
            }

            const double integer_estimate =
                INTEGER_DELAY +
                best_offset;

            const double integer_error =
                integer_estimate -
                true_delay;

            sum_error_integer +=
                integer_error;

            sum_sq_integer +=
                integer_error *
                integer_error;

            max_abs_integer =
                std::max(
                    max_abs_integer,
                    std::abs(integer_error));

            const int expected_offset =
                (true_frac >= 0.5)
                    ? 1
                    : 0;

            if (best_offset !=
                expected_offset)
            {
                ++wrong_integer;
            }

            // ------------------------------------------------
            // Parabolic interpolation
            //
            // EXACTLY as D3a
            // ------------------------------------------------

            double parabola = 0.0;

            if (best_offset > -SEARCH &&
                best_offset < SEARCH)
            {
                const double sm =
                    score[
                        best_offset -
                        1 +
                        SEARCH];

                const double s0 =
                    score[
                        best_offset +
                        SEARCH];

                const double sp =
                    score[
                        best_offset +
                        1 +
                        SEARCH];

                const double denominator =
                    sm -
                    2.0 * s0 +
                    sp;

                if (std::abs(denominator) >
                    1e-20)
                {
                    parabola =
                        0.5 *
                        (sm - sp) /
                        denominator;
                }
            }

            const double parabolic_estimate =
                INTEGER_DELAY +
                best_offset +
                parabola;

            const double parabolic_error =
                parabolic_estimate -
                true_delay;

            sum_error_parabola +=
                parabolic_error;

            sum_sq_parabola +=
                parabolic_error *
                parabolic_error;

            max_abs_parabola =
                std::max(
                    max_abs_parabola,
                    std::abs(parabolic_error));
        }

        // ----------------------------------------------------
        // Final statistics
        // ----------------------------------------------------

        const double mean_integer =
            sum_error_integer /
            NTRIALS;

        const double mean_parabola =
            sum_error_parabola /
            NTRIALS;

        const double rms_integer =
            std::sqrt(
                sum_sq_integer /
                NTRIALS);

        const double rms_parabola =
            std::sqrt(
                sum_sq_parabola /
                NTRIALS);

        // ----------------------------------------------------
        // Results
        // ----------------------------------------------------

        printf("\n");
        printf("------------------------------------------------------------\n");
        printf("SNR = %+5.1f dB\n", snr_db);
        printf("------------------------------------------------------------\n");

        printf("\nInteger-grid estimator:\n");

        printf("  Mean error       = %10.6f samples\n",
               mean_integer);

        printf("  RMS error        = %10.6f samples\n",
               rms_integer);

        printf("  RMS error        = %10.6f ms\n",
               rms_integer *
               1000.0 /
               SAMPLE_RATE);

        printf("  Max abs error    = %10.6f samples\n",
               max_abs_integer);

        printf("  Wrong integer    = %d / %d\n",
               wrong_integer,
               NTRIALS);

        printf("\nParabolic estimator:\n");

        printf("  Mean error       = %10.6f samples\n",
               mean_parabola);

        printf("  RMS error        = %10.6f samples\n",
               rms_parabola);

        printf("  RMS error        = %10.6f ms\n",
               rms_parabola *
               1000.0 /
               SAMPLE_RATE);

        printf("  Max abs error    = %10.6f samples\n",
               max_abs_parabola);

        if (rms_parabola > 0.0)
        {
            printf("\n");
            printf("Improvement        = %10.3fx\n",
                   rms_integer /
                   rms_parabola);
        }
    }

    // --------------------------------------------------------
    // Resolution
    // --------------------------------------------------------

    printf("\n");
    printf("One sample         : %.9f ms\n",
           1000.0 / SAMPLE_RATE);

    printf("One tenth sample   : %.9f ms\n",
           100.0 / SAMPLE_RATE);

    printf("\n");
    printf("D3b completed.\n");
    printf("============================================================\n");
}

// ============================================================
// D3c - FT8 DELAY ESTIMATOR: WINDOW SIZE TEST
//
// Based EXACTLY on D3b.
//
// Tests the effect of the correlation window length:
//
//     WIN = 64
//     WIN = 128
//     WIN = 256
//     WIN = 512
//
// SNR:
//
//     +10 dB
//       0 dB
//
// Everything else is kept identical to D3b.
//
// Purpose:
//   Determine whether increasing the correlation window
//   improves delay accuracy.
// ============================================================

static void test_D3c()
{
    printf("\n");
    printf("============================================================\n");
    printf("D3c - FT8 DELAY ESTIMATOR: CORRELATION WINDOW TEST\n");
    printf("============================================================\n");

    constexpr int SAMPLE_RATE = 12000;
    constexpr int NSYM        = 79;
    constexpr int SPS         = 1920;

    constexpr double FREQ      = 1001.0;
    constexpr double AMPLITUDE = 1000.0;

    constexpr int INTEGER_DELAY = 24000;
    constexpr int SEARCH        = 15;

    constexpr int FIR_TAPS = 129;
    constexpr double LPF_CUTOFF = 150.0;

    constexpr int NTRIALS = 200;

    const int windows[] =
    {
         64,
        128,
        256,
        512
    };

    constexpr int NWINDOWS =
        sizeof(windows) / sizeof(windows[0]);

    const double snr_db_list[] =
    {
        10.0,
         0.0
    };

    constexpr int NSNR =
        sizeof(snr_db_list) / sizeof(snr_db_list[0]);

    // --------------------------------------------------------
    // FT8 tones
    // --------------------------------------------------------

    uint8_t tones[79];
    make_test_tones(tones);

    // --------------------------------------------------------
    // Generate integer-aligned reference waveform
    // --------------------------------------------------------

    const size_t nsamples =
        INTEGER_DELAY +
        NSYM * SPS +
        SAMPLE_RATE;

    std::vector<float> base(nsamples, 0.0f);

    generate_ft8_continuous(
        base.data(),
        base.size(),
        tones,
        FREQ,
        0.0,
        AMPLITUDE,
        (double)INTEGER_DELAY / SAMPLE_RATE);

    // --------------------------------------------------------
    // FIR
    // --------------------------------------------------------

    std::vector<float> fir_coeff(FIR_TAPS);

    {
        const int M = FIR_TAPS - 1;
        const double fc =
            LPF_CUTOFF / SAMPLE_RATE;

        for (int n = 0; n < FIR_TAPS; ++n)
        {
            const double x =
                n - M / 2.0;

            double h;

            if (std::abs(x) < 1e-12)
            {
                h = 2.0 * fc;
            }
            else
            {
                h =
                    std::sin(2.0 * PI * fc * x) /
                    (PI * x);
            }

            const double w =
                0.5 *
                (1.0 -
                 std::cos(2.0 * PI * n / M));

            fir_coeff[n] =
                (float)(h * w);
        }

        double sum = 0.0;

        for (float h : fir_coeff)
            sum += h;

        for (float &h : fir_coeff)
            h = (float)(h / sum);
    }

    // --------------------------------------------------------
    // RF -> complex baseband -> FIR
    //
    // EXACTLY as D3a / D3b
    // --------------------------------------------------------

    auto downconvert =
        [&](const std::vector<float> &rf,
            std::vector<std::complex<float>> &bb)
    {
        bb.assign(
            rf.size(),
            std::complex<float>(0.0f, 0.0f));

        const double w =
            2.0 * PI * FREQ / SAMPLE_RATE;

        std::vector<float> mix_i(rf.size());
        std::vector<float> mix_q(rf.size());

        for (size_t n = 0; n < rf.size(); ++n)
        {
            const double phase =
                w * n;

            mix_i[n] =
                (float)(rf[n] *
                        std::cos(phase));

            mix_q[n] =
                (float)(-rf[n] *
                        std::sin(phase));
        }

        for (size_t n = 0; n < rf.size(); ++n)
        {
            double si = 0.0;
            double sq = 0.0;

            const size_t kmax =
                std::min<size_t>(
                    FIR_TAPS,
                    n + 1);

            for (size_t k = 0;
                 k < kmax;
                 ++k)
            {
                si +=
                    fir_coeff[k] *
                    mix_i[n - k];

                sq +=
                    fir_coeff[k] *
                    mix_q[n - k];
            }

            bb[n] =
                std::complex<float>(
                    (float)si,
                    (float)sq);
        }
    };

    // --------------------------------------------------------
    // Costas transitions
    //
    // EXACTLY as D3a / D3b
    // --------------------------------------------------------

    struct Transition
    {
        int sample;
        int tone_before;
        int tone_after;
    };

    Transition transitions[22];
    int ntrans = 0;

    for (int k = 0;
         k < NSYM - 1;
         ++k)
    {
        if (tones[k] == tones[k + 1])
            continue;

        const int boundary =
            INTEGER_DELAY +
            (k + 1) * SPS;

        transitions[ntrans++] =
        {
            boundary,
            tones[k],
            tones[k + 1]
        };
    }

    assert(ntrans == 22);

    printf("Sample rate        : %d Hz\n",
           SAMPLE_RATE);

    printf("Carrier            : %.3f Hz\n",
           FREQ);

    printf("LPF cutoff         : %.1f Hz\n",
           LPF_CUTOFF);

    printf("FIR taps           : %d\n",
           FIR_TAPS);

    printf("FIR group delay    : %d samples\n",
           FIR_TAPS / 2);

    printf("\n");

    printf("Costas transitions : %d\n",
           ntrans);

    printf("Search             : +/- %d samples\n",
           SEARCH);

    printf("Trials             : %d\n",
           NTRIALS);

    // --------------------------------------------------------
    // RNG
    //
    // Reset for every window/SNR combination so that
    // comparisons use exactly the same random sequence.
    // --------------------------------------------------------

    // --------------------------------------------------------
    // Test each SNR
    // --------------------------------------------------------

    for (int isnr = 0;
         isnr < NSNR;
         ++isnr)
    {
        const double snr_db =
            snr_db_list[isnr];

        printf("\n");
        printf("============================================================\n");
        printf("SNR = %+5.1f dB\n", snr_db);
        printf("============================================================\n");

        printf("\n");
        printf(" WIN       Integer RMS       Parabolic RMS       Improvement\n");
        printf("----------------------------------------------------------------\n");

        for (int iw = 0;
             iw < NWINDOWS;
             ++iw)
        {
            const int WIN =
                windows[iw];

            // ------------------------------------------------
            // Reference baseband
            //
            // EXACTLY as D3a / D3b
            // ------------------------------------------------

            std::vector<std::complex<float>>
                ref(ntrans * WIN);

            {
                std::vector<std::complex<float>> bb;

                downconvert(base, bb);

                const int fir_delay =
                    FIR_TAPS / 2;

                for (int t = 0;
                     t < ntrans;
                     ++t)
                {
                    const int center =
                        transitions[t].sample +
                        fir_delay;

                    const int start =
                        center - WIN / 2;

                    for (int i = 0;
                         i < WIN;
                         ++i)
                    {
                        const int idx =
                            start + i;

                        if (idx >= 0 &&
                            idx < (int)bb.size())
                        {
                            ref[t * WIN + i] =
                                bb[idx];
                        }
                        else
                        {
                            ref[t * WIN + i] =
                                std::complex<float>(
                                    0.0f,
                                    0.0f);
                        }
                    }
                }
            }

            // ------------------------------------------------
            // RNG reset
            // ------------------------------------------------

            std::mt19937 rng(1234567);

            std::uniform_real_distribution<double>
                frac_dist(0.0, 1.0);

            // ------------------------------------------------
            // Statistics
            // ------------------------------------------------

            double sum_sq_integer = 0.0;
            double sum_sq_parabola = 0.0;

            double max_abs_integer = 0.0;
            double max_abs_parabola = 0.0;

            int wrong_integer = 0;

            // ------------------------------------------------
            // Monte Carlo
            // ------------------------------------------------

            for (int trial = 0;
                 trial < NTRIALS;
                 ++trial)
            {
                const double true_frac =
                    frac_dist(rng);

                const double true_delay =
                    INTEGER_DELAY +
                    true_frac;

                // --------------------------------------------
                // Fractional delay
                // --------------------------------------------

                std::vector<float> delayed(
                    nsamples);

                apply_fractional_delay(
                    base.data(),
                    delayed.data(),
                    base.size(),
                    true_frac);

                // --------------------------------------------
                // Signal power
                // --------------------------------------------

                const double signal_power =
                    calculate_signal_power(
                        delayed.data(),
                        delayed.size(),
                        true_delay /
                            SAMPLE_RATE);

                // --------------------------------------------
                // AWGN
                // --------------------------------------------

                const double noise_power =
                    signal_power /
                    std::pow(
                        10.0,
                        snr_db / 10.0);

                const double sigma =
                    std::sqrt(noise_power);

                std::vector<float> noisy =
                    delayed;

                add_awgn(
                    noisy.data(),
                    noisy.size(),
                    sigma,
                    rng);

                // --------------------------------------------
                // RF -> baseband
                // --------------------------------------------

                std::vector<std::complex<float>> bb;

                downconvert(
                    noisy,
                    bb);

                // --------------------------------------------
                // Correlation scores
                // --------------------------------------------

                std::vector<double> score(
                    2 * SEARCH + 1);

                const int fir_delay =
                    FIR_TAPS / 2;

                for (int offset = -SEARCH;
                     offset <= SEARCH;
                     ++offset)
                {
                    double total_score = 0.0;

                    for (int t = 0;
                         t < ntrans;
                         ++t)
                    {
                        const int center =
                            transitions[t].sample +
                            fir_delay +
                            offset;

                        const int start =
                            center - WIN / 2;

                        std::complex<double> corr(
                            0.0,
                            0.0);

                        double rx_power = 0.0;
                        double ref_power = 0.0;

                        for (int i = 0;
                             i < WIN;
                             ++i)
                        {
                            const int idx =
                                start + i;

                            if (idx < 0 ||
                                idx >= (int)bb.size())
                                continue;

                            const std::complex<float> rx =
                                bb[idx];

                            const std::complex<float> rr =
                                ref[t * WIN + i];

                            corr +=
                                std::complex<double>(
                                    rx.real(),
                                    rx.imag()) *
                                std::conj(
                                    std::complex<double>(
                                        rr.real(),
                                        rr.imag()));

                            rx_power +=
                                (double)rx.real() *
                                rx.real() +
                                (double)rx.imag() *
                                rx.imag();

                            ref_power +=
                                (double)rr.real() *
                                rr.real() +
                                (double)rr.imag() *
                                rr.imag();
                        }

                        if (rx_power > 0.0 &&
                            ref_power > 0.0)
                        {
                            const double c2 =
                                corr.real() *
                                corr.real() +
                                corr.imag() *
                                corr.imag();

                            total_score +=
                                c2 /
                                (rx_power *
                                 ref_power);
                        }
                    }

                    score[offset + SEARCH] =
                        total_score;
                }

                // --------------------------------------------
                // Integer maximum
                // --------------------------------------------

                int best_offset = 0;

                double best_score =
                    score[SEARCH];

                for (int offset = -SEARCH;
                     offset <= SEARCH;
                     ++offset)
                {
                    const double s =
                        score[offset + SEARCH];

                    if (s > best_score)
                    {
                        best_score = s;
                        best_offset = offset;
                    }
                }

                const double integer_estimate =
                    INTEGER_DELAY +
                    best_offset;

                const double integer_error =
                    integer_estimate -
                    true_delay;

                sum_sq_integer +=
                    integer_error *
                    integer_error;

                max_abs_integer =
                    std::max(
                        max_abs_integer,
                        std::abs(integer_error));

                const int expected_offset =
                    (true_frac >= 0.5)
                        ? 1
                        : 0;

                if (best_offset !=
                    expected_offset)
                {
                    ++wrong_integer;
                }

                // --------------------------------------------
                // Parabolic interpolation
                // --------------------------------------------

                double parabola = 0.0;

                if (best_offset > -SEARCH &&
                    best_offset < SEARCH)
                {
                    const double sm =
                        score[
                            best_offset -
                            1 +
                            SEARCH];

                    const double s0 =
                        score[
                            best_offset +
                            SEARCH];

                    const double sp =
                        score[
                            best_offset +
                            1 +
                            SEARCH];

                    const double denominator =
                        sm -
                        2.0 * s0 +
                        sp;

                    if (std::abs(denominator) >
                        1e-20)
                    {
                        parabola =
                            0.5 *
                            (sm - sp) /
                            denominator;
                    }
                }

                const double parabolic_estimate =
                    INTEGER_DELAY +
                    best_offset +
                    parabola;

                const double parabolic_error =
                    parabolic_estimate -
                    true_delay;

                sum_sq_parabola +=
                    parabolic_error *
                    parabolic_error;

                max_abs_parabola =
                    std::max(
                        max_abs_parabola,
                        std::abs(parabolic_error));
            }

            // ------------------------------------------------
            // Statistics
            // ------------------------------------------------

            const double rms_integer =
                std::sqrt(
                    sum_sq_integer /
                    NTRIALS);

            const double rms_parabola =
                std::sqrt(
                    sum_sq_parabola /
                    NTRIALS);

            const double improvement =
                (rms_parabola > 0.0)
                    ? rms_integer /
                      rms_parabola
                    : 0.0;

            printf(
                " %4d      %10.6f         %10.6f          %8.3fx\n",
                WIN,
                rms_integer,
                rms_parabola,
                improvement);

            printf(
                "          (%9.6f ms)       (%9.6f ms)"
                "   wrong=%3d/200\n",
                rms_integer *
                    1000.0 /
                    SAMPLE_RATE,
                rms_parabola *
                    1000.0 /
                    SAMPLE_RATE,
                wrong_integer);

            printf(
                "          max=%9.6f        max=%9.6f\n",
                max_abs_integer,
                max_abs_parabola);
        }
    }

    printf("\n");
    printf("One sample         : %.9f ms\n",
           1000.0 / SAMPLE_RATE);

    printf("\n");
    printf("D3c completed.\n");
    printf("============================================================\n");
}

int main()
{
    printf("FT8 frequency estimator Linux test\n");
    printf("rate = %d Hz\n", rate_);
    printf("block = %d samples\n",
           blocksize(rate_));

    test_frequency_estimator();
    test_frequency_estimator_noise();
    test_delay_estimator_noiseless();

    test_delay_score_debug();

    test_delay_estimator_noise();

    test_D1();
    test_D2();
    // test_D3();
    test_D3a();
    test_D3b();
    test_D3c();
    // Uncomment if you want to generate/check one waveform.
    //
    // test_generate_one();

    return 0;
}

