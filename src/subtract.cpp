#include <cassert>
#include <vector>
#include <complex>
#include <cmath>

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

void synthesize(
    float *dst,
    size_t nsamples,
    // const std::vector<int>& re79,
    const uint8_t * tones,
    const std::vector<double>& amps,
    const std::vector<double>& phases,
    double hz0,
    double off_sec,
    double sign,
    int sample_rate)
{
    int block = blocksize(sample_rate);
    int off0 = std::round(off_sec * sample_rate);

    int ramp = std::round(block * subtract_ramp);

    if (ramp < 1)
        ramp = 1;

    //
    // First symbol initial ramp
    //
    {
        double amp = amps[0];
        double phase = phases[0];
        double freq = hz0 + 6.25 * tones[0];

        double dtheta =
            2.0 * M_PI * freq / sample_rate;

        for (int jj = 0; jj < ramp; jj++)
        {
            double theta =
                phase + jj * dtheta;

            double x =
                amp * std::cos(theta);

            x *= jj / (double)ramp;

            int idx = off0 + jj;

            if (idx >= 0 && idx < (int)nsamples)
                dst[idx] += sign * x;
        }
    }

    //
    // All 79 symbols
    //
    for (int si = 0; si < 79; si++)
    {
        double amp = amps[si];
        double phase = phases[si];

        double freq =
            hz0 + 6.25 * tones[si];

        double dtheta =
            2.0 * M_PI * freq / sample_rate;

        //
        // Steady part
        //
        for (int jj = ramp; jj < block-ramp; jj++)
        {
            double theta =
                phase + jj * dtheta;

            double x =
                amp * std::cos(theta);

            int idx =
                off0 + si * block + jj;

            if (idx >= 0 && idx < (int)nsamples)
                dst[idx] += sign * x;
        }

        //
        // Transition to next symbol
        //
        double theta =
            phase + (block-ramp) * dtheta;

        double freq1;
        double phase1;

        if (si + 1 >= 79)
        {
            freq1 = freq;
            phase1 = phase;
        }
        else
        {
            freq1 =
                hz0 + 6.25 * tones[si+1];

            phase1 =
                phases[si+1];
        }

        double dtheta1 =
            2.0 * M_PI * freq1 / sample_rate;

        //
        // Frequency interpolation
        //
        double inc =
            (dtheta1 - dtheta) /
            (2.0 * ramp);

        //
        // Phase correction
        //
        double actual =
            theta +
            dtheta * 2.0 * ramp +
            inc * 4.0 * ramp * ramp / 2.0;

        double target =
            phase1 + dtheta1 * ramp;

        while (std::fabs(target - actual) > M_PI)
        {
            if (target < actual)
                target += 2.0 * M_PI;
            else
                target -= 2.0 * M_PI;
        }

        double adj =
            target - actual;

        int end = block + ramp;

        if (si == 78)
            end = block;

        for (int jj = block-ramp; jj < end; jj++)
        {
            int idx =
                off0 + si * block + jj;

            if (idx >= 0 && idx < (int)nsamples)
            {
                double x =
                    amp * std::cos(theta);

                //
                // Last symbol fade out
                //
                if (si == 78)
                {
                    x *= 1.0 -
                        (jj - (block-ramp)) /
                        (double)ramp;
                }

                dst[idx] += sign * x;
            }

            theta += dtheta;
            dtheta += inc;
            theta += adj / (2.0 * ramp);
        }
    }
}

void subtract(const uint8_t * tones,
              double hz0,
              double hz1,
              double off_sec,
              float * samples_,
              size_t num_samples,
              int sample_rate)
{
    int block = blocksize(sample_rate);
    int off0 = std::round(off_sec * sample_rate);

    std::vector<double> phases(79);
    std::vector<double> amps(79);

    //
    // Estimate amplitudes and phases
    //
    for (int i = 0; i < 79; i++)
    {
        
        double freq =
            hz0 + 6.25 * tones[i];

        double dtheta =
            2.0 * M_PI * freq / sample_rate;

        std::complex<double> c(0.0, 0.0);

        for (int n = 0; n < block; n++)
        {
            int idx = off0 + i * block + n;

            if (idx < 0 || idx >= num_samples)
                continue;

            double theta =
                -dtheta * n;

            std::complex<float> lo(
                std::cos(theta),
                std::sin(theta));

            c += samples_[idx] * lo;
        }

        amps[i] =
            2.0 * std::abs(c) / block;

        phases[i] =
            std::arg(c);

        // printf("symbol %2d: amp=%10.4f phase=%12.8f\n",
        //    i,
        //    amps[i],
        //    phases[i]);
    }

    double amp_err2 = 0.0;
    double phase_err2 = 0.0;

    for (int i = 0; i < 79; ++i)
    {
        double da = amps[i] - 1000.0;

        double dp = phases[i];

        while (dp > M_PI)
            dp -= 2.0 * M_PI;

        while (dp < -M_PI)
            dp += 2.0 * M_PI;

        amp_err2 += da * da;
        phase_err2 += dp * dp;
    }

    double amp_rmse =
        std::sqrt(amp_err2 / 79.0);

    double phase_rmse =
        std::sqrt(phase_err2 / 79.0);

    printf("amp RMSE   = %.6f\n", amp_rmse);
    printf("phase RMSE = %.9f rad (%.6f deg)\n",
        phase_rmse,
        phase_rmse * 180.0 / M_PI);

    //
    // Reconstruct and subtract
    //
    synthesize(
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


// void
// subtract_ori(const std::vector<int> re79,
//          double hz0,
//          double hz1,
//          double off_sec)
// {
//     int block = blocksize(rate_);
//     int off0 = round(off_sec * rate_);

//     std::vector<double> phases(79);
//     std::vector<double> amps(79);
//     //
//     // Estimate the 79 amplitudes and phases using coherent correlation
//     // at the exact FT8 tone frequency.
//     //
//     for(int i = 0; i < 79; i++)
//     {
//         double freq = hz0 + 6.25 * re79[i];
//         double dtheta = 2.0 * M_PI * freq / rate_;
//         std::complex<double> c(0.0, 0.0);
//         for(int n = 0; n < block; n++)
//         {
//             double theta = -dtheta * n;
//             std::complex<float> lo(cos(theta), sin(theta));
//             c += samples_[off0 + i*block + n] * lo;
//         }

//         amps[i] = 2.0 * abs(c) / block;
//         phases[i] = arg(c);
//     }

//     int ramp = round(block * subtract_ramp);
//     if(ramp < 1)
//         ramp = 1;

//     //
//     // First symbol initial ramp
//     //
//     {
//         double amp = amps[0];
//         double phase = phases[0];
//         double freq = hz0 + 6.25 * re79[0];
//         double dtheta = 2.0 * M_PI * freq / rate_;

//         for(int jj = 0; jj < ramp; jj++)
//         {
//             double theta = phase + jj*dtheta;
//             double x = amp*cos(theta);
//             x *= jj/(double)ramp;
//             int idx = off0 + jj;
//             samples_[idx] -= x;
//         }
//     }

//     //
//     // All the 79 symbols
//     //
//     for(int si = 0; si < 79; si++)
//     {
//         double amp = amps[si];
//         double phase = phases[si];
//         double freq = hz0 + 6.25 * re79[si];
//         double dtheta = 2.0 * M_PI * freq / rate_;

//         //
//         // steady part
//         //
//         for(int jj = ramp; jj < block-ramp; jj++)
//         {
//             double theta = phase + jj*dtheta;
//             double x = amp*cos(theta);
//             int idx = off0 + si*block + jj;
//             samples_[idx] -= x;
//         }

//         //
//         // transition to next symbol
//         //
//         double theta = phase + (block-ramp)*dtheta;
//         double freq1;
//         double phase1;

//         if(si+1 >= 79)
//         {
//             freq1 = freq;
//             phase1 = phase;
//         }
//         else
//         {
//             freq1 = hz0 + 6.25 * re79[si+1];
//             phase1 = phases[si+1];
//         }

//         double dtheta1 =
//             2.0 * M_PI * freq1 / rate_;

//         //
//         // frequency interpolation
//         //
//         double inc =
//             (dtheta1 - dtheta)/(2.0*ramp);

//         //
//         // phase correction
//         //
//         double actual =
//             theta +
//             dtheta*2.0*ramp +
//             inc*4.0*ramp*ramp/2.0;

//         double target =
//             phase1 + dtheta1*ramp;

//         while(fabs(target-actual) > M_PI)
//         {
//             if(target < actual)
//                 target += 2*M_PI;
//             else
//                 target -= 2*M_PI;
//         }

//         double adj =
//             target - actual;

//         int end = block + ramp;

//         if(si == 78)
//             end = block;

//         for(int jj = block-ramp; jj < end; jj++)
//         {
//             int idx = off0 + si*block + jj;
//             double x = amp*cos(theta);
//             //
//             // last symbol fade out
//             //
//             if(si == 78)
//             {
//                 x *= 1.0 -
//                     (jj-(block-ramp))/(double)ramp;
//             }

//             samples_[idx] -= x;

//             theta += dtheta;
//             dtheta += inc;
//             theta += adj/(2.0*ramp);
//         }
//     }
// }
