#include <cassert>
#include <vector>
#include <complex>
#include <cmath>
#include "subtract.h"

// int rate_ = 12000;  // samples/second
double subtract_ramp = 0.11;

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

    /*
     * Use float here.
     *
     * ESP32-S3 has hardware float support but double is much
     * more expensive.
     */
    float amps[NTONES];
    float phases[NTONES];


    // =========================================================
    // Estimate amplitudes and phases
    // =========================================================

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
                idx < (int)num_samples)
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
            (float)block;

        amps[i] =
            amp;

        phases[i] =
            atan2f(cq, ci);
    }


    // =========================================================
    // Reconstruct and subtract
    // =========================================================

    
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
}



complex_amp_t estimate_amplitude_phase(
    float *samples,
    int num_samples,
    uint8_t *tones,
    float delay,
    float freq)
{
    constexpr int sample_rate = 12000;
    constexpr int ntones = 79;

    const int block = blocksize(sample_rate);
    const int waveform_samples = ntones * block;

    std::vector<double> amps(ntones, 1.0);
    std::vector<double> phases_i(ntones, 0.0);
    std::vector<double> phases_q(ntones, M_PI / 2.0);

    std::vector<float> reference_i(waveform_samples, 0.0f);
    std::vector<float> reference_q(waveform_samples, 0.0f);

    synthesize(
        reference_i.data(),
        waveform_samples,
        tones,
        amps,
        phases_i,
        freq,
        0.0,
        +1.0,
        sample_rate);

    synthesize(
        reference_q.data(),
        waveform_samples,
        tones,
        amps,
        phases_q,
        freq,
        0.0,
        +1.0,
        sample_rate);

    const int delay_samples =
        (int)std::lround(delay * sample_rate);

    double ci = 0.0;
    double cq = 0.0;

    double energy_i = 0.0;
    double energy_q = 0.0;

    for (int n = 0; n < waveform_samples; ++n)
    {
        const int index = delay_samples + n;

        if (index < 0 || index >= num_samples)
            continue;

        const double x = samples[index];

        ci += x * reference_i[n];
        cq += x * reference_q[n];

        energy_i += reference_i[n] * reference_i[n];
        energy_q += reference_q[n] * reference_q[n];
    }

    //
    // Both references should have essentially the same energy.
    //
    const double energy =
        0.5 * (energy_i + energy_q);

    //
    // Normalize correlation.
    //
    const double i = ci / energy;
    const double q = cq / energy;

    complex_amp_t result;

    result.ci = i;
    result.cq = q;

    result.amplitude =
        std::sqrt(i * i + q * q);

    result.phase =
        std::atan2(q, i);

    return result;
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
     * Search range
     * =========================================================
     */

    const int first_delay =
        std::max(
            0,
            (int)std::lround(
                (delay0 - 0.300f) * sample_rate));

    const int last_delay =
        std::min(
            num_samples - 1,
            (int)std::lround(
                (delay0 + 0.300f) * sample_rate));

    const int ndelays =
        last_delay - first_delay + 1;


    /*
     * =========================================================
     * Oscillator information for the 79 tones
     * =========================================================
     */

    struct Osc
    {
        float pc;
        float ps;

        float end_c;
        float end_s;

        float cd;
        float sd;
    };

    Osc osc[ntones];

    constexpr float TWO_PI =
        6.2831853071795864769f;

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

        osc[k].cd = cd;
        osc[k].sd = sd;

        /*
         * exp(+j*w)
         */
        osc[k].pc = cd;
        osc[k].ps = sd;

        /*
         * exp(-j*w*(block-1))
         */
        float ec = 1.0f;
        float es = 0.0f;

        for (int n = 0;
             n < block - 1;
             ++n)
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
     * =========================================================
     * Complex correlations
     * =========================================================
     */

    float ci[ntones];
    float cq[ntones];

    for (int k = 0; k < ntones; ++k)
    {
        ci[k] = 0.0f;
        cq[k] = 0.0f;

        const int start =
            first_delay + k * block;

        if (start < 0 ||
            start >= num_samples)
            continue;

        const int count =
            std::min(
                block,
                num_samples - start);

        float c = 1.0f;
        float s = 0.0f;

        const float cd =
            osc[k].cd;

        const float sd =
            osc[k].sd;

        for (int n = 0;
             n < count;
             ++n)
        {
            const float x =
                samples[start + n];

            ci[k] += x * c;
            cq[k] -= x * s;

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


    /*
     * =========================================================
     * Calculate score
     *
     * IMPORTANT:
     * This is now only called during the sparse scan and
     * during the final local exact scan.
     * =========================================================
     */

    auto calculate_score =
        [&]() -> float
    {
        float score = 0.0f;
        int used_symbols = 0;

        for (int k = 0;
             k < ntones;
             ++k)
        {
            const int start =
                first_delay +
                k * block;

            if (start >= num_samples)
                break;

            score +=
                ci[k] * ci[k] +
                cq[k] * cq[k];

            ++used_symbols;
        }

        if (used_symbols == 0)
            return 0.0f;

        return score /
               (float)used_symbols;
    };


    /*
     * =========================================================
     * FIRST SCORE
     * =========================================================
     */

    float best_score =
        calculate_score();

    int best_delay =
        first_delay;


    /*
     * =========================================================
     * SPARSE SLIDING-CORRELATION SCAN
     *
     * We still update the correlation for EVERY sample.
     *
     * But we calculate the expensive 79-tone score only
     * every SCORE_STEP samples.
     *
     * 10 samples = 0.833 ms
     * =========================================================
     */

    constexpr int SCORE_STEP = 10;

    for (int d = 1;
         d < ndelays;
         ++d)
    {
        const int delay =
            first_delay + d;

        /*
         * Update all 79 correlations.
         */
        for (int k = 0;
             k < ntones;
             ++k)
        {
            const int start =
                delay +
                k * block;

            if (start >= num_samples)
                break;

            const int new_index =
                start +
                block -
                1;

            if (new_index >= num_samples)
                break;

            /*
             * Remove old sample.
             */
            const float x_old =
                samples[start - 1];

            /*
             * Add new sample.
             */
            const float x_new =
                samples[new_index];

            /*
             * Rotate existing correlation by +w.
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
             * Add incoming sample.
             */
            ci[k] =
                rotated_re +
                x_new * osc[k].end_c;

            cq[k] =
                rotated_im +
                x_new * osc[k].end_s;
        }


        /*
         * -----------------------------------------------------
         * Sparse score evaluation.
         * -----------------------------------------------------
         */

        if ((d % SCORE_STEP) == 0)
        {
            const float score =
                calculate_score();

            if (score > best_score)
            {
                best_score =
                    score;

                best_delay =
                    delay;
            }
        }
    }


    /*
     * =========================================================
     * FINAL EXACT SEARCH
     *
     * The sparse scan found the peak to within SCORE_STEP
     * samples.
     *
     * Now search every sample around it.
     *
     * Since we no longer have the correlations corresponding
     * to best_delay, we reconstruct them for the local window.
     *
     * Search:
     *
     *     best_delay +/- SCORE_STEP
     *
     * at full 1-sample resolution.
     * =========================================================
     */

    const int refine_first =
        std::max(
            first_delay,
            best_delay - SCORE_STEP);

    const int refine_last =
        std::min(
            last_delay,
            best_delay + SCORE_STEP);


    /*
     * =========================================================
     * Recalculate correlation at refine_first.
     * =========================================================
     */

    for (int k = 0;
         k < ntones;
         ++k)
    {
        ci[k] = 0.0f;
        cq[k] = 0.0f;

        const int start =
            refine_first +
            k * block;

        if (start < 0 ||
            start >= num_samples)
            continue;

        const int count =
            std::min(
                block,
                num_samples - start);

        float c = 1.0f;
        float s = 0.0f;

        const float cd =
            osc[k].cd;

        const float sd =
            osc[k].sd;

        for (int n = 0;
             n < count;
             ++n)
        {
            const float x =
                samples[start + n];

            ci[k] += x * c;
            cq[k] -= x * s;

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


    /*
     * =========================================================
     * EXACT LOCAL SEARCH
     * =========================================================
     */

    float final_score =
        calculate_score();

    int final_delay =
        refine_first;

    for (int delay =
             refine_first + 1;
         delay <= refine_last;
         ++delay)
    {
        /*
         * Slide correlation by one sample.
         */
        for (int k = 0;
             k < ntones;
             ++k)
        {
            const int start =
                delay +
                k * block;

            if (start >= num_samples)
                break;

            const int new_index =
                start +
                block -
                1;

            if (new_index >= num_samples)
                break;

            const float x_old =
                samples[start - 1];

            const float x_new =
                samples[new_index];

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

            ci[k] =
                rotated_re +
                x_new * osc[k].end_c;

            cq[k] =
                rotated_im +
                x_new * osc[k].end_s;
        }

        const float score =
            calculate_score();

        if (score > final_score)
        {
            final_score =
                score;

            final_delay =
                delay;
        }
    }


    /*
     * =========================================================
     * Result
     * =========================================================
     */

    return final_delay /
           (float)sample_rate;
}

struct ft8_freq_moments_t
{
    float re[11];
    float im[11];
};

static void ft8_prepare_frequency_moments(
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

float refine_ft8_frequency(
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
     *
     * This is the only expensive pass.
     */
    ft8_prepare_frequency_moments(
        samples,
        num_samples,
        tones,
        delay,
        freq_coarse,
        moments);

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