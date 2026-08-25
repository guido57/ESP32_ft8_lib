#include <cstdio>
#include <vector>
#include <cmath>
#include "ft8/encode.h"
#include "ft8/decode.h"
#include "ft8/message.h"
#include "subtract.h"
#include "common/wave.h"
#include "cstring"

extern "C" void encode174(const uint8_t* message, uint8_t* codeword);


static double rms(const float *p, size_t n)
{
    double e = 0.0;

    for (size_t i = 0; i < n; ++i)
        e += double(p[i]) * double(p[i]);

    return std::sqrt(e / n);
}

#define FT8_CRC_WIDTH      (14)
#define TOPBIT (1u << (FT8_CRC_WIDTH - 1))

// Compute 14-bit CRC for a sequence of given number of bits
// Adapted from https://barrgroup.com/Embedded-Systems/How-To/CRC-Calculation-C-Code
// [IN] message  - byte sequence (MSB first)
// [IN] num_bits - number of bits in the sequence
uint16_t ft8lib_compute_crc(const uint8_t message[], int num_bits)
{
    uint16_t remainder = 0;
    int idx_byte = 0;

    // Perform modulo-2 division, a bit at a time.
    for (int idx_bit = 0; idx_bit < num_bits; ++idx_bit)
    {
        if (idx_bit % 8 == 0)
        {
            // Bring the next byte into the remainder.
            remainder ^= (message[idx_byte] << (FT8_CRC_WIDTH - 8));
            ++idx_byte;
        }

        // Try to divide the current data bit.
        if (remainder & TOPBIT)
        {
            remainder = (remainder << 1) ^ FT8_CRC_POLYNOMIAL;
        }
        else
        {
            remainder = (remainder << 1);
        }
    }

    return remainder & ((TOPBIT << 1) - 1u);
}



void ft8lib_add_crc(const uint8_t payload[], uint8_t a91[])
{
    // Copy 77 bits of payload data
    for (int i = 0; i < 10; i++)
        a91[i] = payload[i];

    // Clear 3 bits after the payload to make 82 bits
    a91[9] &= 0xF8u;
    a91[10] = 0;

    // Calculate CRC of 82 bits (77 + 5 zeros)
    // 'The CRC is calculated on the source-encoded message, zero-extended from 77 to 82 bits'
    uint16_t checksum = ft8lib_compute_crc(a91, 96 - 14);

    // Store the CRC at the end of 77 bit message
    a91[9] |= (uint8_t)(checksum >> 11);
    a91[10] = (uint8_t)(checksum >> 3);
    a91[11] = (uint8_t)(checksum << 5);
}


void ft8lib_encode(const uint8_t* payload, uint8_t* tones)
{
    uint8_t a91[FTX_LDPC_K_BYTES]; // Store 77 bits of payload + 14 bits CRC

    // Compute and add CRC at the end of the message
    // a91 contains 77 bits of payload + 14 bits of CRC
    ft8lib_add_crc(payload, a91);

    uint8_t codeword[FTX_LDPC_N_BYTES];
    encode174(a91, codeword);

    // Message structure: S7 D29 S7 D29 S7
    // Total symbols: 79 (FT8_NN)

    uint8_t mask = 0x80u; // Mask to extract 1 bit from codeword
    int i_byte = 0;       // Index of the current byte of the codeword
    for (int i_tone = 0; i_tone < FT8_NN; ++i_tone)
    {
        if ((i_tone >= 0) && (i_tone < 7))
        {
            tones[i_tone] = kFT8_Costas_pattern[i_tone];
        }
        else if ((i_tone >= 36) && (i_tone < 43))
        {
            tones[i_tone] = kFT8_Costas_pattern[i_tone - 36];
        }
        else if ((i_tone >= 72) && (i_tone < 79))
        {
            tones[i_tone] = kFT8_Costas_pattern[i_tone - 72];
        }
        else
        {
            // Extract 3 bits from codeword at i-th position
            uint8_t bits3 = 0;

            if (codeword[i_byte] & mask)
                bits3 |= 4;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }
            if (codeword[i_byte] & mask)
                bits3 |= 2;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }
            if (codeword[i_byte] & mask)
                bits3 |= 1;
            if (0 == (mask >>= 1))
            {
                mask = 0x80u;
                i_byte++;
            }

            tones[i_tone] = kFT8_Gray_map[bits3];
        }
    }
}


static void test_synthesize_ft8()
{
    const char *msg = "CQ IW5ALZ JN53";
    constexpr int sample_rate = 12000;
    constexpr int nsamples = 15 * sample_rate;
    const double hz0 = 1500.0;
    const double offset = 2.5; // 2500 ms
    const double amplitude = 0.5;
    const char *wav_path = "test_real_ft8.wav";

    printf("\nSynthesize the FT8 message %s at sample rate %d with time_offset %.3f into %s\n", msg, sample_rate, offset, wav_path);

    
    std::vector<float> samples(nsamples, 0.0f);

    // --------------------------------------------------
    // Encode real FT8 message
    // --------------------------------------------------

    ftx_message_t m;

    if (ftx_message_encode(&m, nullptr, msg) != FTX_MESSAGE_RC_OK)
    {
        printf("ERROR: cannot encode message: %s\n", msg);
        return;
    }

    uint8_t tones[79] = {
    3, 1, 4, 0, 6, 5, 2,
    0, 3, 1, 7, 4, 5, 2, 6, 4, 5, 0, 5, 4, 7, 6, 7, 0, 4, 6, 0, 6, 0, 2, 1, 4, 3, 2, 0, 5,
    3, 1, 4, 0, 6, 5, 2,
    6, 4, 0, 4, 0, 1, 3, 6, 5, 0, 5, 4, 5, 4, 5, 0, 7, 0, 6, 4, 0, 4, 1, 1, 4, 0, 0, 4, 2,
    3, 1, 4, 0, 6, 5, 2
    };
    // ft8lib_encode(m.payload, tones);

    printf("Message : %s\n", msg);
    printf("Rate    : %d Hz\n", sample_rate);
    printf("Frequency: %.3f Hz\n", hz0);
    printf("Offset  : %.3f ms\n", offset * 1000.0);
    printf("Amplitude: %.3f\n", amplitude);

    printf("Tones:\n");

    for (int i = 0; i < 79; ++i)
    {
        printf("%d%s",
               tones[i],
               (i == 78) ? "\n" : " ");
    }

    // --------------------------------------------------
    // Synthesize
    // --------------------------------------------------

    std::vector<double> amps(79, amplitude);
    std::vector<double> phases(79, 0.0);

    synthesize(
        samples.data(),
        samples.size(),
        tones,
        amps,
        phases,
        hz0,
        offset,
        +1.0,
        sample_rate);

    // --------------------------------------------------
    // Save WAV
    // --------------------------------------------------
    
    save_wav(
        samples.data(),
        nsamples,
        sample_rate,
        wav_path);

    printf("WAV saved: %s\n", wav_path);

    // --------------------------------------------------
    // Basic signal check
    // --------------------------------------------------

    double rms_signal = 0.0;

    int off0 = std::lround(offset * sample_rate);

    if (off0 >= 0 && off0 + 79 * 1920 <= nsamples)
    {
        rms_signal =
            rms(samples.data() + off0, 79 * 1920);
    }

    printf("Signal RMS: %.6f\n", rms_signal);

    printf("PASS\n");
}

int main()
{
    
    test_synthesize_ft8();

    printf("\nDONE\n");

    return 0;
}