#include <cassert>
#include <vector>
#include <complex>
#include <cmath>
#include "subtract.h"
#include <cstdio>

// int rate_ = 12000;  // samples/second
double subtract_ramp = 0.11;

extern double elapsed_ms(const struct timespec *t0, const struct timespec *t1);


// return symbol length in samples at the given rate.
// insist on integer symbol lengths so that we can
// use whole FFT bins. e.g. at 12000 samples/second, symbol length is 1920 samples.
int
blocksize(int rate)
{
  // FT8 symbol length is 1920 at 12000 samples/second.
  int xblock = 1920 / (12000.0 / rate);
  assert(xblock == (int) xblock);
  int block = xblock;
  return block;
}

// Synthesize an FT8 time waveform based on the given tones, amplitudes, and phases, 
// and add it to the destination buffer.
// Synthesize an FT8 time waveform based on the given tones,
// amplitudes, and phases, and add it to the destination buffer.
//
// ESP32-S3 optimized:
//   - float arithmetic in the inner loops
//   - recursive oscillator in constant-frequency sections
//   - no cos() in the large steady-state section
//   - original transition equations preserved
void synthesize(
    float *dst,
    size_t nsamples,
    const uint8_t *tones,
    const std::vector<double>& amps,
    const std::vector<double>& phases,
    double hz0,
    double off_sec,
    double sign,
    int sample_rate)
{
    constexpr float TWO_PI =
        6.2831853071795864769f;

    const int block =
        blocksize(sample_rate);

    const int off0 =
        (int)lroundf(
            (float)off_sec *
            (float)sample_rate);

    int ramp =
        (int)lroundf(
            (float)block *
            (float)subtract_ramp);

    if (ramp < 1)
        ramp = 1;

    const float sign_f =
        (float)sign;

    const float fs =
        (float)sample_rate;


    // =========================================================
    // First symbol initial ramp
    // =========================================================

    {
        const float amp =
            (float)amps[0];

        const float phase =
            (float)phases[0];

        const float freq =
            (float)hz0 +
            6.25f * (float)tones[0];

        const float dtheta =
            TWO_PI * freq / fs;

        float c =
            cosf(phase);

        float s =
            sinf(phase);

        const float cd =
            cosf(dtheta);

        const float sd =
            sinf(dtheta);

        for (int jj = 0;
             jj < ramp;
             ++jj)
        {
            const int idx =
                off0 + jj;

            if (idx >= 0 &&
                idx < (int)nsamples)
            {
                float x =
                    amp * c;

                x *=
                    (float)jj /
                    (float)ramp;

                dst[idx] +=
                    sign_f * x;
            }

            /*
             * Recursive oscillator.
             */
            const float nc =
                c * cd -
                s * sd;

            const float ns =
                s * cd +
                c * sd;

            c = nc;
            s = ns;
        }
    }


    // =========================================================
    // All 79 symbols
    // =========================================================

    for (int si = 0;
         si < 79;
         ++si)
    {
        const float amp =
            (float)amps[si];

        const float phase =
            (float)phases[si];

        /*
         * Current symbol frequency.
         */
        const float freq =
            (float)hz0 +
            6.25f * (float)tones[si];

        /*
         * IMPORTANT:
         * dtheta must NOT be const because it is modified
         * during the transition section.
         */
        float dtheta =
            TWO_PI * freq / fs;


        // =====================================================
        // Steady part
        // =====================================================

        /*
         * Original code:
         *
         * theta = phase + jj*dtheta
         *
         * Therefore at jj=ramp:
         *
         * theta = phase + ramp*dtheta
         */
        const float theta_start =
            phase +
            (float)ramp * dtheta;

        float c =
            cosf(theta_start);

        float s =
            sinf(theta_start);

        const float cd =
            cosf(dtheta);

        const float sd =
            sinf(dtheta);

        for (int jj = ramp;
             jj < block - ramp;
             ++jj)
        {
            const int idx =
                off0 +
                si * block +
                jj;

            if (idx >= 0 &&
                idx < (int)nsamples)
            {
                dst[idx] +=
                    sign_f *
                    amp *
                    c;
            }

            /*
             * Advance oscillator by one sample.
             */
            const float nc =
                c * cd -
                s * sd;

            const float ns =
                s * cd +
                c * sd;

            c = nc;
            s = ns;
        }


        // =====================================================
        // Transition to next symbol
        // =====================================================

        /*
         * IMPORTANT:
         *
         * From here onward we retain the original algorithm.
         */
        float theta =
            phase +
            (float)(block - ramp) *
            dtheta;

        float freq1;
        float phase1;

        if (si + 1 >= 79)
        {
            freq1 =
                freq;

            phase1 =
                phase;
        }
        else
        {
            freq1 =
                (float)hz0 +
                6.25f *
                (float)tones[si + 1];

            phase1 =
                (float)phases[si + 1];
        }

        const float dtheta1 =
            TWO_PI *
            freq1 /
            fs;


        // -----------------------------------------------------
        // Frequency interpolation
        // -----------------------------------------------------

        const float inc =
            (dtheta1 - dtheta) /
            (2.0f * (float)ramp);


        // -----------------------------------------------------
        // Phase correction
        // -----------------------------------------------------

        const float actual =
            theta +
            dtheta *
            2.0f *
            (float)ramp +
            inc *
            4.0f *
            (float)ramp *
            (float)ramp /
            2.0f;

        float target =
            phase1 +
            dtheta1 *
            (float)ramp;

        while (
            fabsf(target - actual) >
            3.14159265358979323846f)
        {
            if (target < actual)
                target += TWO_PI;
            else
                target -= TWO_PI;
        }

        const float adj =
            target - actual;

        const float adj_step =
            adj /
            (2.0f * (float)ramp);


        int end =
            block + ramp;

        if (si == 78)
            end = block;


        // =====================================================
        // Transition waveform
        // =====================================================

        for (int jj = block - ramp;
             jj < end;
             ++jj)
        {
            const int idx =
                off0 +
                si * block +
                jj;

            if (idx >= 0 &&
                idx < (int)nsamples)
            {
                float x =
                    amp *
                    cosf(theta);

                /*
                 * Last symbol fade out.
                 */
                if (si == 78)
                {
                    x *=
                        1.0f -
                        (float)(jj - (block - ramp)) /
                        (float)ramp;
                }

                dst[idx] +=
                    sign_f * x;
            }

            /*
             * Original phase/frequency evolution.
             */
            theta += dtheta;
            dtheta += inc;
            theta += adj_step;
        }
    }
}


static void synthesize_float(
    float *dst,
    size_t nsamples,
    const uint8_t *tones,
    const float *amps,
    const float *phases,
    float hz0,
    float off_sec,
    float sign,
    int sample_rate)
{

    // printf(
    //     "SYNTH START: off=%.6f hz=%.3f nsamples=%u\n",
    //     off_sec,
    //     hz0,
    //     (unsigned)nsamples);


    constexpr float TWO_PI =
        6.2831853071795864769f;

    const int block =
        blocksize(sample_rate);

    const int off0 =
        (int)lroundf(
            off_sec * (float)sample_rate);

    int ramp =
        (int)lroundf(
            (float)block *
            (float)subtract_ramp);

    if (ramp < 1)
        ramp = 1;


    // =========================================================
    // First symbol initial ramp
    // =========================================================

    {
        const float amp =
            amps[0];

        const float phase =
            phases[0];

        const float freq =
            hz0 +
            6.25f * (float)tones[0];

        const float dtheta =
            TWO_PI * freq / (float)sample_rate;

        float c =
            cosf(phase);

        float s =
            sinf(phase);

        const float cd =
            cosf(dtheta);

        const float sd =
            sinf(dtheta);

        for (int jj = 0;
             jj < ramp;
             ++jj)
        {
            const int idx =
                off0 + jj;

            if (idx >= 0 &&
                idx < (int)nsamples)
            {
                float x =
                    amp * c;

                x *=
                    (float)jj /
                    (float)ramp;

                dst[idx] +=
                    sign * x;
            }

            const float nc =
                c * cd -
                s * sd;

            const float ns =
                s * cd +
                c * sd;

            c = nc;
            s = ns;
        }
    }


    // =========================================================
    // All 79 symbols
    // =========================================================
    // printf("Synthesize_float: Synthesizing and subtracting symbols...\n");

    for (int si = 0;
         si < 79;
         ++si)
    {
        const float amp =
            amps[si];

        const float phase =
            phases[si];

        const float freq =
            hz0 +
            6.25f * (float)tones[si];

        float dtheta =
            TWO_PI * freq /
            (float)sample_rate;


        // =====================================================
        // Steady part
        // =====================================================

        const float theta_start =
            phase +
            (float)ramp * dtheta;

        float c_steady =
            cosf(theta_start);

        float s_steady =
            sinf(theta_start);

        const float cd =
            cosf(dtheta);

        const float sd =
            sinf(dtheta);

        for (int jj = ramp;
             jj < block - ramp;
             ++jj)
        {
            const int idx =
                off0 +
                si * block +
                jj;

            if (idx >= 0 &&
                idx < (int)nsamples)
            {
                dst[idx] +=
                    sign *
                    amp *
                    c_steady;
            }

            const float nc =
                c_steady * cd -
                s_steady * sd;

            const float ns =
                s_steady * cd +
                c_steady * sd;

            c_steady = nc;
            s_steady = ns;
        }


        // =====================================================
        // Transition to next symbol
        // =====================================================

        float theta =
            phase +
            (float)(block - ramp) *
            dtheta;

        float freq1;
        float phase1;

        if (si == 78)
        {
            freq1 =
                freq;

            phase1 =
                phase;
        }
        else
        {
            freq1 =
                hz0 +
                6.25f *
                (float)tones[si + 1];

            phase1 =
                phases[si + 1];
        }

        const float dtheta1 =
            TWO_PI * freq1 /
            (float)sample_rate;

        const float inc =
            (dtheta1 - dtheta) /
            (2.0f * (float)ramp);


        // -----------------------------------------------------
        // Phase correction
        // -----------------------------------------------------

        const float actual =
            theta +
            dtheta *
            2.0f *
            (float)ramp +
            inc *
            4.0f *
            (float)ramp *
            (float)ramp /
            2.0f;

        float target =
            phase1 +
            dtheta1 *
            (float)ramp;

        // while (
        //     fabsf(target - actual) >
        //     3.14159265358979323846f)
        // {
        //     if (target < actual)
        //         target += TWO_PI;
        //     else
        //         target -= TWO_PI;
        // }

        float phase_error = target - actual;

        if (!std::isfinite(target) || !std::isfinite(actual) || !std::isfinite(phase_error)) {
            printf("ERROR synth nonfinite phase: "
                "si=%d target=%f actual=%f error=%f "
                "freq=%f freq1=%f amp=%f\n",
                si, target, actual, phase_error,
                freq, freq1, amp);
            return;
        }

        float adj = remainderf(phase_error, TWO_PI);
        float adj_step = adj / (2.0f * ramp);


        // -----------------------------------------------------
        // Number of transition samples
        // -----------------------------------------------------

        int end =
            block + ramp;

        if (si == 78)
            end =
                block;

        const int transition_samples =
            end - (block - ramp);


        // =====================================================
        // Recursive oscillator
        //
        // Original algorithm:
        //
        //   x = cos(theta)
        //   theta += dtheta
        //   dtheta += inc
        //   theta += adj_step
        //
        // Therefore:
        //
        //   delta = dtheta + adj_step
        //
        // and delta increases by 'inc' every sample.
        // =====================================================

        float delta =
            dtheta + adj_step;

        float c_trans =
            cosf(theta);

        float s_trans =
            sinf(theta);

        // Rotation corresponding to current delta.
        float rc =
            cosf(delta);

        float rs =
            sinf(delta);

        // Rotation corresponding to +inc.
        const float ric =
            cosf(inc);

        const float ris =
            sinf(inc);


        // =====================================================
        // Transition loop
        //
        // NO sinf()/cosf() here.
        // =====================================================

        for (int jj = 0;
             jj < transition_samples;
             ++jj)
        {
            const int idx =
                off0 +
                si * block +
                (block - ramp) +
                jj;

            if (idx >= 0 &&
                idx < (int)nsamples)
            {
                float x =
                    amp * c_trans;

                // Last-symbol fade out
                if (si == 78)
                {
                    x *=
                        1.0f -
                        (float)jj /
                        (float)ramp;
                }

                dst[idx] +=
                    sign * x;
            }


            // -------------------------------------------------
            // Advance waveform phase:
            //
            // theta_next =
            //     theta + delta
            // -------------------------------------------------

            const float nc =
                c_trans * rc -
                s_trans * rs;

            const float ns =
                s_trans * rc +
                c_trans * rs;

            c_trans =
                nc;

            s_trans =
                ns;


            // -------------------------------------------------
            // delta_next = delta + inc
            //
            // Rotate the oscillator increment by 'inc'.
            // -------------------------------------------------

            const float nrc =
                rc * ric -
                rs * ris;

            const float nrs =
                rs * ric +
                rc * ris;

            rc =
                nrc;

            rs =
                nrs;
        }
    }
}

void subtract(const uint8_t *tones,
              double hz0,
              double hz1,
              double off_sec,
              float *samples_,
              size_t num_samples,
              int sample_rate)
{
    constexpr int NTONES = 79;
    constexpr float TWO_PI =
        6.2831853071795864769f;

    const int block =
        blocksize(sample_rate);

    const int off0 =
        (int)lroundf(
            (float)off_sec *
            (float)sample_rate);

    const float fs =
        (float)sample_rate;

    float amps[NTONES];
    float phases[NTONES];

    // =========================================================
    // Estimate amplitudes and phases
    // =========================================================
    // printf("Subtract: estimating amplitudes and phases...\n");
     /*
     * Ignore the first and last 10% of every symbol.
     *
     * These regions contain the FT8 symbol ramps/transitions,
     * whereas the central 80% is steady-state.
     */
    const int first = block / 10;
    const int last  = block - block / 10;
    const int corr_len = last - first;

    for (int i = 0; i < NTONES; ++i)
    {
        const float freq =
            (float)hz0 +
            6.25f * (float)tones[i];

        const float dtheta =
            TWO_PI * freq / fs;

        /*
         * Start oscillator at theta = 0.
         *
         * We need exp(-j*dtheta*n):
         *
         *   real = cos(theta)
         *   imag = -sin(theta)
         */
        const float cd =
            cosf(dtheta);

        const float sd =
            sinf(dtheta);

        float c = 1.0f;
        float s = 0.0f;

        float ci = 0.0f;
        float cq = 0.0f;

        const int start =
            off0 + i * block;

        for (int n = 0; n < block; ++n)
        {
            const int idx =
                start + n;

            if (idx >= 0 &&
                idx < (int)num_samples &&
                n >= first &&
                n < last)
            {
                const float x =
                    samples_[idx];

                /*
                 * x * exp(-j theta)
                 */
                ci += x * c;
                cq -= x * s;
            }

            /*
             * Rotate oscillator.
             */
            const float nc =
                c * cd -
                s * sd;

            const float ns =
                s * cd +
                c * sd;

            c = nc;
            s = ns;
        }

        /*
         * FT8 real-signal amplitude.
         */
        const float amp =
            2.0f *
            sqrtf(
                ci * ci +
                cq * cq) /
            (float)corr_len;

        amps[i] =
            amp;

        phases[i] =
            atan2f(cq, ci);
    }


    // =========================================================
    // Reconstruct and subtract
    // =========================================================

    // printf("Reconstructing and subtracting...\n");
    synthesize_float(
        samples_,
        num_samples,
        tones,
        amps,
        phases,
        hz0,
        off_sec,
        -1.0,
        sample_rate);
    // printf("Subtraction complete.\n");
}

void diagnose_subtraction(
    float *samples,
    size_t num_samples,
    const uint8_t *tones,
    double freq,
    double delay,
    int sample_rate)
{
    constexpr int NTONES = 79;
    constexpr float TWO_PI =
        6.2831853071795864769f;

    const int block =
        blocksize(sample_rate);

    const int delay_samples =
        (int)lround(delay * sample_rate);

    float amps[NTONES];
    float phases[NTONES];

    const int first =
        block / 10;

    const int last =
        block - block / 10;

    const int corr_len =
        last - first;

    /*
     * =========================================================
     * 1. Estimate amplitude and phase
     * =========================================================
     *
     * This is exactly the same estimator used by subtract().
     */
    for (int i = 0; i < NTONES; ++i)
    {
        const float f =
            (float)freq +
            6.25f * (float)tones[i];

        const float dtheta =
            TWO_PI * f /
            (float)sample_rate;

        const float cd =
            cosf(dtheta);

        const float sd =
            sinf(dtheta);

        float c = 1.0f;
        float s = 0.0f;

        float ci = 0.0f;
        float cq = 0.0f;

        const int start =
            delay_samples +
            i * block;

        for (int n = 0;
             n < block;
             ++n)
        {
            const int idx =
                start + n;

            if (idx >= 0 &&
                idx < (int)num_samples &&
                n >= first &&
                n < last)
            {
                const float x =
                    samples[idx];

                ci += x * c;
                cq -= x * s;
            }

            const float nc =
                c * cd -
                s * sd;

            const float ns =
                s * cd +
                c * sd;

            c = nc;
            s = ns;
        }

        amps[i] =
            2.0f *
            sqrtf(ci * ci + cq * cq) /
            (float)corr_len;

        phases[i] =
            atan2f(cq, ci);
    }


    /*
     * =========================================================
     * 2. Generate candidate waveform
     * =========================================================
     */

    std::vector<float> reference(
        num_samples,
        0.0f);

    synthesize_float(
        reference.data(),
        reference.size(),
        tones,
        amps,
        phases,
        (float)freq,
        (float)delay,
        +1.0f,
        sample_rate);


    /*
     * =========================================================
     * 3. Measure BEFORE subtraction
     * =========================================================
     *
     * x = measured signal
     * s = candidate waveform
     *
     * We want to know how much of x lies in the direction of s.
     *
     * alpha = <x,s>/<s,s>
     *
     * If alpha = 1:
     *     candidate amplitude is correct.
     *
     * If alpha < 1:
     *     candidate is too large relative to actual signal.
     *
     * If alpha > 1:
     *     candidate is too small.
     *
     * Note: this assumes reference phase/amplitude are already
     * estimated from x, so alpha should normally be very close
     * to 1. Its deviation is therefore an excellent diagnostic.
     */

    const int waveform_start =
        std::max(0, delay_samples);

    const int waveform_end =
        std::min(
            delay_samples +
            NTONES * block,
            (int)num_samples);

    double xx = 0.0;
    double ss = 0.0;
    double xs = 0.0;

    for (int n = waveform_start;
         n < waveform_end;
         ++n)
    {
        const double x =
            samples[n];

        const double s =
            reference[n];

        xx += x * x;
        ss += s * s;
        xs += x * s;
    }

    /*
     * Projection coefficient.
     *
     * This is much more useful than RMS reduction.
     */
    const double alpha =
        (ss > 0.0)
        ? xs / ss
        : 0.0;

    /*
     * Correlation coefficient.
     */
    const double rho_before =
        (xx > 0.0 && ss > 0.0)
        ? xs / sqrt(xx * ss)
        : 0.0;


    /*
     * =========================================================
     * 4. Calculate the component orthogonal to the candidate
     * =========================================================
     *
     * x = alpha*s + e
     *
     * where e is orthogonal to s.
     *
     * Therefore:
     *
     * candidate_energy = alpha^2 * ss
     *
     * residual_error_energy =
     *     xx - alpha^2 * ss
     *
     * This tells us how much of the measured signal is NOT
     * explained by the candidate.
     */

    double candidate_energy =
        alpha * alpha * ss;

    double error_energy =
        xx - candidate_energy;

    if (error_energy < 0.0)
        error_energy = 0.0;


    /*
     * =========================================================
     * 5. Actual subtraction
     * =========================================================
     */

    std::vector<float> residual(
        samples,
        samples + num_samples);

    for (int n = waveform_start;
         n < waveform_end;
         ++n)
    {
        residual[n] -= reference[n];
    }


    /*
     * =========================================================
     * 6. Measure AFTER subtraction
     * =========================================================
     */

    double rr = 0.0;
    double rs = 0.0;

    for (int n = waveform_start;
         n < waveform_end;
         ++n)
    {
        const double r =
            residual[n];

        const double s =
            reference[n];

        rr += r * r;
        rs += r * s;
    }


    /*
     * Remaining projection of residual onto candidate.
     *
     * Ideally:
     *
     *     residual = x - s
     *
     * and if x contains exactly s:
     *
     *     residual projection = 0
     */
    const double alpha_after =
        (ss > 0.0)
        ? rs / ss
        : 0.0;


    /*
     * Correlation after subtraction.
     */
    const double rho_after =
        (rr > 0.0 && ss > 0.0)
        ? rs / sqrt(rr * ss)
        : 0.0;


    /*
     * =========================================================
     * 7. Candidate cancellation
     * =========================================================
     *
     * Before subtraction, the candidate component is:
     *
     *     alpha
     *
     * After subtraction:
     *
     *     alpha_after
     *
     * Therefore the cancellation ratio is:
     *
     *     |alpha_after / alpha|
     *
     * and cancellation in dB:
     *
     *     20 log10(...)
     *
     * This is the quantity we really care about.
     */

    double cancellation_ratio =
        0.0;

    if (fabs(alpha) > 1e-12)
    {
        cancellation_ratio =
            fabs(alpha_after / alpha);
    }

    double cancellation_db;

    if (cancellation_ratio > 0.0)
    {
        cancellation_db =
            20.0 *
            log10(cancellation_ratio);
    }
    else
    {
        cancellation_db =
            -INFINITY;
    }


    /*
     * =========================================================
     * 8. Subtraction error relative to the candidate
     * =========================================================
     *
     * Ideal:
     *
     *     residual candidate component = 0
     *
     * Define:
     *
     *     relative_error =
     *         |alpha_after| / |alpha|
     *
     * This is easier to interpret than total RMS.
     */

    const double relative_error =
        cancellation_ratio;


    /*
     * =========================================================
     * 9. Total energy reduction
     *
     * Keep this only as secondary information.
     * =========================================================
     */

    const double energy_reduction_db =
        (rr > 0.0 && xx > 0.0)
        ? 10.0 * log10(rr / xx)
        : -INFINITY;


    /*
     * =========================================================
     * 10. Candidate-to-residual ratio
     * =========================================================
     *
     * This estimates how much stronger the candidate component
     * was than the remaining signal.
     */

    const double candidate_rms =
        sqrt(candidate_energy /
             (double)(waveform_end -
                      waveform_start));

    const double residual_rms =
        sqrt(rr /
             (double)(waveform_end -
                      waveform_start));

    const double candidate_to_residual_db =
        (residual_rms > 0.0)
        ? 20.0 *
          log10(candidate_rms /
                residual_rms)
        : INFINITY;


    /*
     * =========================================================
     * REPORT
     * =========================================================
     */

    printf("\n");
    printf("============================================================\n");
    printf("SUBTRACTION DIAGNOSTIC\n");
    printf("============================================================\n");

    printf("Frequency              = %.6f Hz\n", freq);
    printf("Delay                  = %.6f s\n", delay);
    printf("Samples analyzed       = %d\n",
           waveform_end - waveform_start);

    printf("\n");
    printf("INPUT\n");
    printf("------------------------------------------------------------\n");

    printf("Signal energy          = %.9e\n", xx);
    printf("Reference energy       = %.9e\n", ss);

    printf("Projection alpha       = %.9f\n",
           alpha);

    printf("Correlation rho        = %.9f\n",
           rho_before);

    printf("Candidate energy       = %.9e\n",
           candidate_energy);

    printf("Unexplained energy     = %.9e\n",
           error_energy);

    printf("\n");
    printf("AFTER SUBTRACTION\n");
    printf("------------------------------------------------------------\n");

    printf("Residual energy        = %.9e\n",
           rr);

    printf("Residual projection    = %.9f\n",
           alpha_after);

    printf("Residual correlation   = %.9f\n",
           rho_after);

    printf("Energy reduction       = %.3f dB\n",
           energy_reduction_db);

    printf("\n");
    printf("CANDIDATE CANCELLATION\n");
    printf("------------------------------------------------------------\n");

    printf("Initial projection     = %.9f\n",
           alpha);

    printf("Final projection       = %.9f\n",
           alpha_after);

    printf("Remaining candidate    = %.6f %%\n",
           100.0 * relative_error);

    printf("Candidate cancellation = %.3f dB\n",
           cancellation_db);

    printf("Candidate/residual     = %.3f dB\n",
           candidate_to_residual_db);

    printf("============================================================\n");

    /*
     * =========================================================
     * INTERPRETATION
     * =========================================================
     */

    printf("\nINTERPRETATION\n");

    if (fabs(alpha_after) < 0.01 * fabs(alpha))
    {
        printf("  EXCELLENT: >99%% of candidate component removed.\n");
    }
    else if (fabs(alpha_after) < 0.05 * fabs(alpha))
    {
        printf("  GOOD: >95%% of candidate component removed.\n");
    }
    else if (fabs(alpha_after) < 0.10 * fabs(alpha))
    {
        printf("  FAIR: >90%% of candidate component removed.\n");
    }
    else if (fabs(alpha_after) < 0.25 * fabs(alpha))
    {
        printf("  POOR: significant candidate component remains.\n");
    }
    else
    {
        printf("  BAD: candidate component largely remains.\n");
    }

    /*
     * Detect possible over-subtraction.
     *
     * If alpha_after has the opposite sign to alpha, the
     * subtraction crossed zero and therefore overshot.
     */
    if (alpha * alpha_after < 0.0)
    {
        printf("  WARNING: subtraction appears to OVER-SHOOT.\n");
    }

    printf("============================================================\n\n");
}

float refine_ft8_delay(
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

    /*
     * =========================================================
     * Search range (unchanged: full +/-0.100 s window)
     * =========================================================
     */

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

    constexpr float TWO_PI =
        6.2831853071795864769f;

    /*
     * =========================================================
     * Per-delay score, accumulated one tone at a time.
     *
     * IMPORTANT (perf, not accuracy):
     * The original implementation interleaved all 79 tones for
     * every delay, so consecutive memory accesses jumped between
     * 79 widely-separated regions of `samples`. On a PSRAM-backed
     * buffer (as used on ESP32-S3) that access pattern thrashes
     * the external-memory cache and is drastically slower than on
     * a normal CPU/cache.
     *
     * Here we process one tone at a time across the whole delay
     * range, so for a fixed tone the samples read advance by one
     * element per delay step - a purely sequential scan that
     * PSRAM/cache hardware can prefetch/burst efficiently.
     *
     * This still touches exactly the same samples, with the same
     * 79-tone correlation and the same +/-0.300 s window, and now
     * scores every single delay exactly (no sparse-scan/local
     * refine approximation needed), so accuracy is at least as
     * good as before.
     * =========================================================
     */

    std::vector<float> score_total(ndelays, 0.0f);
    std::vector<int> used_count(ndelays, 0);

    for (int k = 0; k < ntones; ++k)
    {
        const float freq_hz =
            freq + 6.25f * (float)tones[k];

        const float w =
            TWO_PI * freq_hz /
            (float)sample_rate;

        const float cd =
            cosf(w);

        const float sd =
            sinf(w);

        /*
         * exp(-j*w*(block-1)), rotation applied to the incoming sample.
         */
        float end_c = 1.0f;
        float end_s = 0.0f;

        for (int n = 0;
             n < block - 1;
             ++n)
        {
            const float nc =
                end_c * cd + end_s * sd;

            const float ns =
                end_s * cd - end_c * sd;

            end_c = nc;
            end_s = ns;
        }

        /*
         * Correlation at d = 0 (delay = first_delay), from scratch.
         */
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

            score_total[0] +=
                ci * ci + cq * cq;

            ++used_count[0];
        }

        /*
         * Sequential slide for d = 1 .. ndelays-1: for this tone the
         * addresses touched only increase, one sample at a time.
         */
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

    /*
     * =========================================================
     * Exact argmax over the full 1-sample-resolution window.
     * =========================================================
     */

    int best_delay = first_delay;
    float best_score = -1.0f;

    for (int d = 0; d < ndelays; ++d)
    {
        if (used_count[d] == 0)
            continue;

        const float score =
            score_total[d] / (float)used_count[d];

        if (score > best_score)
        {
            best_score = score;
            best_delay = first_delay + d;
        }
    }

    return best_delay /
           (float)sample_rate;
}



float refine_ft8_delay_27(
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

    /*
     * ---------------------------------------------------------
     * Search range
     *
     * The original algorithm does:
     *
     *   coarse : delay0 +/- 300 ms
     *   medium : around coarse result
     *   fine   : around medium result
     *
     * Therefore ALL candidate delays are contained inside
     * delay0 +/- 300 ms.
     *
     * this complete range once. (now +/- 100 ms instead of 300 ms)
     * this complete range once.
     * ---------------------------------------------------------
     */

    const int first_delay =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.100f) * sample_rate));

    const int last_delay =
        std::min(
            num_samples - 1,
            (int)std::lround(
                (delay0 + 0.100f) * sample_rate));

    const int ndelays =
        last_delay - first_delay + 1;

    /*
     * Scores for every integer-sample delay.
     *
     * At 12 kHz and +/-300 ms this is only about 7201 floats.
     */
    std::vector<float> scores(ndelays, 0.0);

    /*
     * ---------------------------------------------------------
     * Precalculate oscillator coefficients for all 79 tones.
     * ---------------------------------------------------------
     */

    struct Osc
    {
        float cd;
        float sd;

        /*
         * exp(+j*w)
         */
        float pc;
        float ps;

        /*
         * exp(-j*w*(block-1))
         */
        float end_c;
        float end_s;
    };

    Osc osc[ntones];

    constexpr float TWO_PI =
        2.0 * M_PI;

    for (int k = 0; k < ntones; ++k)
    {
        const float freq_hz =
            freq + 6.25f * tones[k];

        const float w =
            TWO_PI * freq_hz / sample_rate;

        const float cd =
            std::cos(w);

        const float sd =
            std::sin(w);

        osc[k].cd = cd;
        osc[k].sd = sd;

        /*
         * exp(+j*w)
         */
        osc[k].pc = cd;
        osc[k].ps = sd;

        /*
         * exp(-j*w*(block-1))
         *
         * Calculate recursively rather than calling sin/cos
         * again.
         */
        float ec = 1.0f;
        float es = 0.0f;

        for (int n = 0; n < block - 1; ++n)
        {
            const float nc =
                ec * cd + es * sd;

            const float ns =
                es * cd - ec * sd;

            ec = nc;
            es = ns;
        }

        osc[k].end_c = ec;
        osc[k].end_s = es;
    }

    /*
     * ---------------------------------------------------------
     * Complex correlations.
     *
     * ci[k] + j*cq[k]
     *
     * corresponds exactly to the old:
     *
     *   ci += x*c
     *   cq -= x*s
     *
     * ---------------------------------------------------------
     */

    float ci[ntones];
    float cq[ntones];

    for (int k = 0; k < ntones; ++k)
    {
        ci[k] = 0.0;
        cq[k] = 0.0;

        const int start =
            first_delay + k * block;

        if (start < 0 ||
            start >= num_samples)
            continue;

        const int count =
            std::min(
                block,
                num_samples - start);

        /*
         * Initial correlation.
         *
         * This is the only full 1920-sample calculation for
         * this symbol.
         */
        float c = 1.0f;
        float s = 0.0f;

        const float cd = osc[k].cd;
        const float sd = osc[k].sd;

        for (int n = 0; n < count; ++n)
        {
            const float x =
                samples[start + n];

            ci[k] += x * c;
            cq[k] -= x * s;

            const float nc =
                c * cd - s * sd;

            const float ns =
                s * cd + c * sd;

            c = nc;
            s = ns;
        }
    }

    /*
     * ---------------------------------------------------------
     * Calculate score at current delay.
     * ---------------------------------------------------------
     */

    auto calculate_score =
        [&]() -> float
    {
        float score = 0.0;
        int used_symbols = 0;

        for (int k = 0; k < ntones; ++k)
        {
            const int start =
                first_delay + k * block;

            if (start >= num_samples)
                break;

            score +=
                ci[k] * ci[k] +
                cq[k] * cq[k];

            ++used_symbols;
        }

        if (used_symbols == 0)
            return 0.0;

        return score / used_symbols;
    };

    /*
     * Score for first delay.
     */
    scores[0] =
        calculate_score();

    /*
     * ---------------------------------------------------------
     * SLIDING CORRELATION
     *
     * Move delay one sample at a time.
     *
     * For:
     *
     * C(d) = sum x[d+n] exp(-j*w*n)
     *
     * we have:
     *
     * C(d+1) =
     *     (C(d) - x[d]) exp(+j*w)
     *     + x[d+B] exp(-j*w*(B-1))
     *
     * ---------------------------------------------------------
     */

    for (int d = 1; d < ndelays; ++d)
    {
        const int delay =
            first_delay + d;

        float score = 0.0;
        int used_symbols = 0;

        for (int k = 0; k < ntones; ++k)
        {
            const int start =
                delay + k * block;

            if (start >= num_samples)
                break;

            /*
             * Old first sample of this window.
             */
            const float x_old =
                samples[start - 1];

            /*
             * New sample entering the window.
             *
             * If the complete block is available this is
             * exactly start + block - 1.
             */
            const int new_index =
                start + block - 1;

            if (new_index >= num_samples)
                break;

            const float x_new =
                samples[new_index];

            /*
             * -------------------------------------------------
             * Remove old sample and rotate the remaining
             * correlation by exp(+j*w).
             *
             * z = ci + j*cq
             *
             * z' = (z - x_old) * exp(+j*w)
             * -------------------------------------------------
             */

            const float r =
                ci[k] - x_old;

            const float im =
                cq[k];

            const float rotated_re =
                r * osc[k].pc -
                im * osc[k].ps;

            const float rotated_im =
                r * osc[k].ps +
                im * osc[k].pc;

            /*
             * Add new sample with phase:
             *
             * exp(-j*w*(block-1))
             */
            ci[k] =
                rotated_re +
                x_new * osc[k].end_c;

            cq[k] =
                rotated_im +
                x_new * osc[k].end_s;

            score +=
                ci[k] * ci[k] +
                cq[k] * cq[k];

            ++used_symbols;
        }

        if (used_symbols > 0)
        {
            scores[d] =
                score / used_symbols;
        }
    }

    /*
     * ---------------------------------------------------------
     * Helper matching the original search().
     *
     * We already have every integer-sample score, so searching
     * now costs essentially nothing.
     * ---------------------------------------------------------
     */

    auto search =
        [&](float center,
            float half_range,
            float step_ms) -> float
    {
        const int first =
            std::max(
                first_delay,
                (int)std::lround(
                    (center - half_range) *
                    sample_rate));

        const int last =
            std::min(
                last_delay,
                (int)std::lround(
                    (center + half_range) *
                    sample_rate));

        const int step =
            std::max(
                1,
                (int)std::lround(
                    step_ms *
                    sample_rate /
                    1000.0));

        float best_score = -1.0;
        int best_delay = first;

        for (int delay = first;
             delay <= last;
             delay += step)
        {
            const int index =
                delay - first_delay;

            const float score =
                scores[index];

            if (score > best_score)
            {
                best_score = score;
                best_delay = delay;
            }
        }

        return best_delay /
               (float)sample_rate;
    };

    /*
     * ---------------------------------------------------------
     * PASS 1: coarse
     * ---------------------------------------------------------
     */

    const float delay1 =
        search(
            delay0,
            0.300f,
            10.0f);

    // printf(
    //     "fine delay coarse = %.6f\n",
    //     delay1);

    /*
     * ---------------------------------------------------------
     * PASS 2: medium
     * ---------------------------------------------------------
     */

    const float delay2 =
        search(
            delay1,
            0.020f,
            1.0f);

    // printf(
    //     "fine delay medium = %.6f\n",
    //     delay2);

    /*
     * ---------------------------------------------------------
     * PASS 3: fine
     *
     * 1 sample = 83.333 us
     * ---------------------------------------------------------
     */

    const float delay3 =
        search(
            delay2,
            0.001f,
            1000.0f / sample_rate);

    // printf(
    //     "fine delay final  = %.6f\n",
    //     delay3);

    return delay3;
}



struct ft8_freq_moments_t
{
    float re[11];
    float im[11];
};

struct ft8_weight_t
{
    float w[8];
};

struct ft8_osc_t
{
    float re;
    float im;
};


static ft8_weight_t weights[1536] __attribute__((aligned(16)));
static ft8_osc_t osc[8][1536] __attribute__((aligned(16)));

static bool weights_initialized = false;


static void ft8_init_moment_weights()
{
    const int block = blocksize(12000);

    const int first = block / 10;
    const int last  = block - block / 10;
    const int nsamp = last - first;

    const float norm =
        1.0f / (float)(nsamp - 1);

    float u = 0.0f;

    for (int i = 0; i < nsamp; i++)
    {
        float up = u;

        for (int p = 0; p < 8; p++)
        {
            weights[i].w[p] = up;
            up *= u;
        }

        u += norm;
    }

    weights_initialized = true;
}

static void ft8_prepare_frequency_moments(
    const float * __restrict samples,
    int num_samples,
    const uint8_t * __restrict tones,
    float delay,
    float freq_ref,
    ft8_freq_moments_t * __restrict moments)
{
    constexpr float sample_rate = 12000.0f;
    constexpr int ntones = 79;
    constexpr int NT = 8;

    constexpr float TWO_PI =
        6.2831853071795864769f;

    const int block = blocksize(12000);

    const int delay_samples =
        (int)lroundf(delay * sample_rate);

    const int first = block / 10;
    const int last  = block - block / 10;
    const int nsamp = last - first;

    /*
     * ---------------------------------------------------------
     * Initialize polynomial weights only once.
     * ---------------------------------------------------------
     */

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    if (!weights_initialized)
        ft8_init_moment_weights();
    clock_gettime(CLOCK_MONOTONIC, &t1);
    printf("FT8 moment weights initialization time: %.3f ms\n", elapsed_ms(&t0, &t1));

    /*
     * ---------------------------------------------------------
     * Generate the 8 oscillators.
     *
     * These depend on freq_ref, so they must be regenerated
     * when freq_ref changes.
     * ---------------------------------------------------------
     */
    
    for (int t = 0; t < NT; t++)
    {
        const float f =
            freq_ref + 6.25f * (float)t;

        const float dtheta =
            TWO_PI * f / sample_rate;

        const float cd = cosf(dtheta);
        const float sd = sinf(dtheta);

        const float theta0 =
            dtheta * (float)first;

        float c = cosf(theta0);
        float s = sinf(theta0);

        ft8_osc_t * __restrict o = osc[t];

        for (int i = 0; i < nsamp; i++)
        {
            o[i].re = c;
            o[i].im = -s;

            const float nc =
                c * cd - s * sd;

            const float ns =
                s * cd + c * sd;

            c = nc;
            s = ns;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t0);
    printf("Oscillator generation time: %.3f ms\n", elapsed_ms(&t1, &t0));
    /*
     * ---------------------------------------------------------
     * Process symbols.
     * ---------------------------------------------------------
     */

    for (int k = 0; k < ntones; k++)
    {
        /*
         * Clear all 11 elements because the structure has
         * 11 slots, even though this routine currently uses
         * only moments 0..8.
         */
        for (int p = 0; p < 11; p++)
        {
            moments[k].re[p] = 0.0f;
            moments[k].im[p] = 0.0f;
        }

        const int start =
            delay_samples + k * block;

        if (start < 0 || start + block > num_samples)
            continue;

        const int tone = tones[k];

        const ft8_osc_t * __restrict o =
            osc[tone];

        const ft8_weight_t * __restrict w =
            weights;

        const float * __restrict x =
            samples + start + first;

        float * __restrict r =
            moments[k].re;

        float * __restrict im =
            moments[k].im;


        /*
         * -----------------------------------------------------
         * Hot loop.
         * -----------------------------------------------------
         */
            float r0 = 0.0f, r1 = 0.0f, r2 = 0.0f;
            float r3 = 0.0f, r4 = 0.0f, r5 = 0.0f;
            float r6 = 0.0f, r7 = 0.0f, r8 = 0.0f;

            float i0 = 0.0f, i1 = 0.0f, i2 = 0.0f;
            float i3 = 0.0f, i4 = 0.0f, i5 = 0.0f;
            float i6 = 0.0f, i7 = 0.0f, i8 = 0.0f;

            for (int i = 0; i < nsamp; i++)
            {
                const float sample = x[i];

                const float yr = sample * o[i].re;
                const float yi = sample * o[i].im;

                const float *wp = weights[i].w;

                 // moment 0: weight = 1
                r0 += yr;
                i0 += yi;
    
                r1 += yr * wp[1];
                r2 += yr * wp[2];
                r3 += yr * wp[3];
                r4 += yr * wp[4];
                r5 += yr * wp[5];
                r6 += yr * wp[6];
                r7 += yr * wp[7];
                r8 += yr * wp[8];

                i1 += yi * wp[1];
                i2 += yi * wp[2];
                i3 += yi * wp[3];
                i4 += yi * wp[4];
                i5 += yi * wp[5];
                i6 += yi * wp[6];
                i7 += yi * wp[7];
                i8 += yi * wp[8];
            }

            r[0] = r0;
            r[1] = r1;
            r[2] = r2;
            r[3] = r3;
            r[4] = r4;
            r[5] = r5;
            r[6] = r6;
            r[7] = r7;
            r[8] = r8;

            im[0] = i0;
            im[1] = i1;
            im[2] = i2;
            im[3] = i3;
            im[4] = i4;
            im[5] = i5;
            im[6] = i6;
            im[7] = i7;
            im[8] = i8;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    printf("FT8 frequency moments computation time: %.3f ms\n", elapsed_ms(&t0, &t1));
}


static void ft8_prepare_frequency_moments_ori(
    const float *samples,
    int num_samples,
    const uint8_t *tones,
    float delay,
    float freq_ref,
    ft8_freq_moments_t *moments)
{
    constexpr float sample_rate = 12000.0f;
    constexpr int ntones = 79;
    constexpr int ORDER = 8;

    const int block = blocksize(12000);

    const int delay_samples =
        (int)lroundf(delay * sample_rate);

    constexpr float TWO_PI =
        6.2831853071795864769f;

    const int first = block / 10;
    const int last  = block - block / 10;

    const float norm =
        1.0f / (float)(last - first - 1);

    for (int k = 0; k < ntones; k++)
    {
        /*
         * Clear moments.
         */
        for (int p = 0; p <= ORDER; p++)
        {
            moments[k].re[p] = 0.0f;
            moments[k].im[p] = 0.0f;
        }

        const int start =
            delay_samples + k * block;

        if (start < 0 || start + block > num_samples)
            continue;

        /*
         * Reference frequency for this FT8 tone.
         */
        const float f =
            freq_ref + 6.25f * (float)tones[k];

        const float dtheta =
            TWO_PI * f / sample_rate;

        /*
         * IMPORTANT:
         * dtheta is constant for this symbol.
         *
         * Calculate the oscillator rotation only once.
         */
        const float cd = cosf(dtheta);
        const float sd = sinf(dtheta);

        /*
         * Oscillator starts at n = first.
         */
        const float theta0 =
            dtheta * (float)first;

        float c = cosf(theta0);
        float s = sinf(theta0);

        for (int n = first; n < last; n++)
        {
            const float x =
                samples[start + n];

            /*
             * Baseband sample:
             *
             * y = x * exp(-j*w*n)
             */
            const float yr = x * c;
            const float yi = -x * s;

            /*
             * Normalized local time.
             */
            const float u =
                (float)(n - first) * norm;

            /*
             * Accumulate powers of u.
             */
            float up = 1.0f;

            for (int p = 0; p <= ORDER; p++)
            {
                moments[k].re[p] += yr * up;
                moments[k].im[p] += yi * up;

                up *= u;
            }

            /*
             * Oscillator recurrence.
             *
             * No sinf()/cosf() here anymore.
             */
            const float nc =
                c * cd - s * sd;

            const float ns =
                s * cd + c * sd;

            c = nc;
            s = ns;
        }
    }
}



static float ft8_frequency_score_moments(
    const ft8_freq_moments_t *moments,
    float delta_freq)
{
    constexpr int ntones = 79;
    constexpr int ORDER = 8;

    constexpr float sample_rate = 12000.0f;
    constexpr float TWO_PI =
        6.2831853071795864769f;

    const float a =
        TWO_PI * delta_freq /
        sample_rate * 1535.0f;

    float cr[ORDER + 1];
    float ci[ORDER + 1];

    cr[0] = 1.0f;
    ci[0] = 0.0f;

    for (int p = 1; p <= ORDER; p++)
    {
        cr[p] =
            ( ci[p - 1] * a) / (float)p;

        ci[p] =
            (-cr[p - 1] * a) / (float)p;
    }

    float score = 0.0f;

    for (int k = 0; k < ntones; k++)
    {
        float sum_re = 0.0f;
        float sum_im = 0.0f;

        for (int p = 0; p <= ORDER; p++)
        {
            const float mr =
                moments[k].re[p];

            const float mi =
                moments[k].im[p];

            sum_re +=
                mr * cr[p] -
                mi * ci[p];

            sum_im +=
                mr * ci[p] +
                mi * cr[p];
        }

        score +=
            sum_re * sum_re +
            sum_im * sum_im;
    }

    return score;
}

/*
 * Refine the frequency estimate for an FT8 signal.
 *
 * samples: the input signal samples
 * num_samples: the number of samples
 * tones: the decoded tones
 * delay: the estimated signal delay
 * freq_coarse: the coarse frequency estimate
 *
 * Returns the refined frequency estimate.
*/
float refine_ft8_frequency_ori(
    const float *samples,
    int num_samples,
    const uint8_t *tones,
    float delay,
    float freq_coarse)
{
    constexpr int ntones = 79;

    /*
     * Moment storage.
     *
     * ~7 KB.
     */
    ft8_freq_moments_t moments[ntones];

    /*
     * Prepare signal at the decoder frequency.
     * This is the only expensive pass.
     */
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    ft8_prepare_frequency_moments(
        samples,
        num_samples,
        tones,
        delay,
        freq_coarse,
        moments);
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    printf("ft8_prepare_frequency_moments took %f seconds\n",
           (end_time.tv_sec - start_time.tv_sec) +
           (end_time.tv_nsec - start_time.tv_nsec) / 1e9);

    /*
     * ---------------------------------------------------------
     * Stage 1
     *
     * +/- 2 Hz
     * 0.20 Hz
     *
     * 21 evaluations
     * ---------------------------------------------------------
     */

    constexpr float search1 = 2.0f;
    constexpr float step1   = 0.20f;

    float best_freq =
        freq_coarse;

    float best_score =
        -1.0f;

    for (float f = freq_coarse - search1;
         f <= freq_coarse + search1 + 0.0001f;
         f += step1)
    {
        const float delta =
            f - freq_coarse;

        const float score =
            ft8_frequency_score_moments(
                moments,
                delta);

        if (score > best_score)
        {
            best_score = score;
            best_freq = f;
        }
    }

    /*
     * ---------------------------------------------------------
     * Stage 2
     *
     * Three-point parabolic interpolation.
     * ---------------------------------------------------------
     */

    constexpr float step2 = 0.20f;

    const float f_minus =
        best_freq - step2;

    const float f_zero =
        best_freq;

    const float f_plus =
        best_freq + step2;

    const float s_minus =
        ft8_frequency_score_moments(
            moments,
            f_minus - freq_coarse);

    const float s_zero =
        ft8_frequency_score_moments(
            moments,
            f_zero - freq_coarse);

    const float s_plus =
        ft8_frequency_score_moments(
            moments,
            f_plus - freq_coarse);

    const float denominator =
        s_minus -
        2.0f * s_zero +
        s_plus;

    float final_freq =
        best_freq;

    if (denominator < 0.0f)
    {
        const float offset =
            0.5f * step2 *
            (s_minus - s_plus) /
            denominator;

        if (offset >= -step2 &&
            offset <= step2)
        {
            final_freq =
                best_freq + offset;
        }
    }

    return final_freq;
}


static float ft8_frequency_score_costas(
    const float *samples,
    int num_samples,
    float delay,
    float freq)
{
    constexpr float sample_rate = 12000.0f;
    constexpr float TWO_PI = 6.2831853071795864769f;

    constexpr int NT = 7;

    // FT8 Costas pattern
    constexpr uint8_t costas[7] =
    {
        3, 1, 4, 0, 6, 5, 2
    };

    const int block = blocksize(12000);

    const int delay_samples =
        (int)lroundf(delay * sample_rate);

    const int first = block / 10;
    const int last  = block - block / 10;
    const int nsamp = last - first;

    float total_re = 0.0f;
    float total_im = 0.0f;

    /*
     * Three Costas arrays:
     *
     *   0..6
     *   36..42
     *   72..78
     */
    constexpr int starts[3] =
    {
        0,
        36,
        72
    };

    for (int c = 0; c < 3; c++)
    {
        const int symbol0 = starts[c];

        for (int j = 0; j < NT; j++)
        {
            const int k = symbol0 + j;

            const int start =
                delay_samples + k * block;

            if (start < 0 ||
                start + block > num_samples)
                continue;

            const int tone = costas[j];

            const float f =
                freq + 6.25f * (float)tone;

            const float dtheta =
                TWO_PI * f / sample_rate;

            const float cd = cosf(dtheta);
            const float sd = sinf(dtheta);

            const float theta0 =
                dtheta * (float)first;

            float cr = cosf(theta0);
            float sr = sinf(theta0);

            const float *x =
                samples + start + first;

            float sum_re = 0.0f;
            float sum_im = 0.0f;

            for (int i = 0; i < nsamp; i++)
            {
                const float sample = x[i];

                sum_re += sample * cr;
                sum_im -= sample * sr;

                const float nr =
                    cr * cd - sr * sd;

                const float ni =
                    sr * cd + cr * sd;

                cr = nr;
                sr = ni;
            }

            total_re += sum_re;
            total_im += sum_im;
        }
    }

    return total_re * total_re +
           total_im * total_im;
}

float refine_ft8_frequency_costas(
    const float *samples,
    int num_samples,
    const uint8_t *tones,
    float delay,
    float freq_coarse)
{
    (void)tones; // Costas tones are known

    constexpr float search1 = 2.0f;
    constexpr float step1   = 0.20f;

    float best_freq  = freq_coarse;
    float best_score = -1.0f;

    /*
     * ---------------------------------------------------------
     * Stage 1
     *
     * Search +/- 2 Hz in 0.20 Hz steps.
     *
     * 21 frequency evaluations.
     * ---------------------------------------------------------
     */

    for (float f = freq_coarse - search1;
         f <= freq_coarse + search1 + 0.0001f;
         f += step1)
    {
        const float score =
            ft8_frequency_score_costas(
                samples,
                num_samples,
                delay,
                f);

        if (score > best_score)
        {
            best_score = score;
            best_freq = f;
        }
    }

    /*
     * ---------------------------------------------------------
     * Stage 2
     *
     * Three-point parabolic interpolation.
     * ---------------------------------------------------------
     */

    constexpr float step2 = 0.20f;

    const float s_minus =
        ft8_frequency_score_costas(
            samples,
            num_samples,
            delay,
            best_freq - step2);

    const float s_zero =
        ft8_frequency_score_costas(
            samples,
            num_samples,
            delay,
            best_freq);

    const float s_plus =
        ft8_frequency_score_costas(
            samples,
            num_samples,
            delay,
            best_freq + step2);

    const float denominator =
        s_minus -
        2.0f * s_zero +
        s_plus;

    float final_freq = best_freq;

    if (denominator < 0.0f)
    {
        const float offset =
            0.5f * step2 *
            (s_minus - s_plus) /
            denominator;

        if (offset >= -step2 &&
            offset <= step2)
        {
            final_freq =
                best_freq + offset;
        }
    }

    return final_freq;
}


static float refine_ft8_frequency_costas_phase(
    const float *samples,
    int num_samples,
    float delay,
    float freq_coarse)
{
    constexpr float sample_rate = 12000.0f;
    constexpr float TWO_PI =
        6.2831853071795864769f;

    constexpr float TONE_SPACING = 6.25f;

    constexpr int NCOSTAS = 21;

    constexpr uint8_t costas[7] =
    {
        3, 1, 4, 0, 6, 5, 2
    };

    constexpr int symbol_numbers[21] =
    {
         0,  1,  2,  3,  4,  5,  6,
        36, 37, 38, 39, 40, 41, 42,
        72, 73, 74, 75, 76, 77, 78
    };

    const int block =
        blocksize((int)sample_rate);

    const int delay_samples =
        (int)lroundf(delay * sample_rate);

    /*
     * Use the central 80% of each symbol.
     */
    const int first = block / 10;
    const int last  = block - block / 10;
    const int nsamp = last - first;

    float phase[NCOSTAS];
    float weight[NCOSTAS];

    /*
     * ------------------------------------------------------------
     * Measure Costas phase for all 21 symbols.
     * ------------------------------------------------------------
     *
     * IMPORTANT:
     *
     * Do NOT add theta0 to the measured phase.
     *
     * The correlation reference already removes the known
     * coarse frequency (including the Costas tone). What remains
     * is essentially the phase caused by the frequency error.
     */
    for (int n = 0; n < NCOSTAS; ++n)
    {
        const int k =
            symbol_numbers[n];

        const int start =
            delay_samples + k * block;

        if (start < 0 ||
            start + block > num_samples)
        {
            phase[n]  = 0.0f;
            weight[n] = 0.0f;
            continue;
        }

        const int costas_index =
            (k < 7) ? k :
            (k < 43) ? (k - 36) :
                       (k - 72);

        const int tone =
            costas[costas_index];

        const float f =
            freq_coarse +
            TONE_SPACING * (float)tone;

        const float dtheta =
            TWO_PI * f / sample_rate;

        /*
         * Reference oscillator starts at the first sample
         * actually used by the correlation.
         */
        float c = cosf(dtheta * (float)first);
        float s = sinf(dtheta * (float)first);

        const float cd = cosf(dtheta);
        const float sd = sinf(dtheta);

        const float *x =
            samples + start + first;

        float sum_re = 0.0f;
        float sum_im = 0.0f;

        for (int i = 0; i < nsamp; ++i)
        {
            const float sample = x[i];

            sum_re += sample * c;
            sum_im -= sample * s;

            const float nc =
                c * cd - s * sd;

            const float ns =
                s * cd + c * sd;

            c = nc;
            s = ns;
        }

        phase[n] =
            atan2f(sum_im, sum_re);

        weight[n] =
            sum_re * sum_re +
            sum_im * sum_im;
    }

    /*
     * ------------------------------------------------------------
     * Estimate frequency independently from each Costas group.
     * ------------------------------------------------------------
     *
     * This is the crucial part.
     *
     * We NEVER unwrap:
     *
     *     symbol 6 -> symbol 36
     *
     * or
     *
     *     symbol 42 -> symbol 72
     *
     * because those intervals are 4.8 seconds long.
     *
     * At only 1 Hz offset, that corresponds to almost 5 complete
     * phase rotations.
     */
    double frequency_sum    = 0.0;
    double frequency_weight = 0.0;

    for (int group = 0; group < 3; ++group)
    {
        const int base = group * 7;

        float unwrapped[7];

        unwrapped[0] =
            phase[base];

        /*
         * Unwrap only the seven consecutive Costas symbols.
         */
        for (int i = 1; i < 7; ++i)
        {
            float d =
                phase[base + i] -
                phase[base + i - 1];

            while (d > (float)M_PI)
                d -= TWO_PI;

            while (d < -(float)M_PI)
                d += TWO_PI;

            unwrapped[i] =
                unwrapped[i - 1] + d;
        }

        /*
         * Weighted linear regression:
         *
         *       phase = a + slope * time
         *
         * frequency offset = slope / (2*pi)
         */
        double sum_w  = 0.0;
        double sum_t  = 0.0;
        double sum_p  = 0.0;
        double sum_tt = 0.0;
        double sum_tp = 0.0;

        for (int i = 0; i < 7; ++i)
        {
            const int n =
                base + i;

            if (weight[n] <= 0.0f)
                continue;

            const double t =
                (double)symbol_numbers[n] *
                ((double)block / sample_rate);

            const double p =
                (double)unwrapped[i];

            const double w =
                (double)weight[n];

            sum_w  += w;
            sum_t  += w * t;
            sum_p  += w * p;
            sum_tt += w * t * t;
            sum_tp += w * t * p;
        }

        const double denominator =
            sum_w * sum_tt -
            sum_t * sum_t;

        if (denominator <= 0.0)
            continue;

        const double slope =
            (sum_w * sum_tp -
             sum_t * sum_p) /
            denominator;

        const double df =
            slope / TWO_PI;

        /*
         * Weight the group estimate according to correlation
         * strength.
         */
        frequency_sum += df * sum_w;
        frequency_weight += sum_w;
    }

    if (frequency_weight <= 0.0)
        return freq_coarse;

    const float delta_freq =
        (float)(frequency_sum / frequency_weight);

    /*
     * Sanity check.
     */
    if (!std::isfinite(delta_freq) ||
        delta_freq < -3.0f ||
        delta_freq > 3.0f)
    {
        return freq_coarse;
    }

    return freq_coarse + delta_freq;
}



float refine_ft8_frequency(
    const float *samples,
    int num_samples,
    const uint8_t *tones,
    float delay,
    float freq_coarse)
{
    (void)tones;

    struct timespec start_time, end_time;

    clock_gettime(CLOCK_MONOTONIC, &start_time);

    const float result =
        refine_ft8_frequency_costas_phase(
            samples,
            num_samples,
            delay,
            freq_coarse);

    clock_gettime(CLOCK_MONOTONIC, &end_time);

    // printf("Costas phase frequency refinement took %f msec\n",
    //        (end_time.tv_sec - start_time.tv_sec) * 1000.0 +
    //        (end_time.tv_nsec - start_time.tv_nsec) / 1e6);

    return result;
}
