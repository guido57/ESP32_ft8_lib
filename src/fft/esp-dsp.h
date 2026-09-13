#pragma once

#include "dsps_fft2r.h"
#include "dsps_fft4r.h"
#include "esp_heap_caps.h"
#include <string.h>

// ============================================================
// Maximum supported REAL FFT size
// ============================================================

#define MAX_FFT_SIZE       4096
#define MAX_FFT_HALF_SIZE  (MAX_FFT_SIZE / 2)

// ============================================================
// FFT2R TABLE
// ============================================================
//
// The internal complex FFT has MAX_FFT_SIZE / 2 points.
//
// dsps_fft2r_init_fc32() needs:
//
//     2 * max_fft_size_complex floats
//
// For a maximum real FFT of 4096:
//
//     complex FFT = 2048 points
//     table = 4096 floats = 16 KB
// ============================================================

#define FFT2_TABLE_FLOATS  (MAX_FFT_HALF_SIZE * 2)
#define FFT2_TABLE_BYTES   (FFT2_TABLE_FLOATS * sizeof(float))

// ============================================================
// FFT4R TABLE
// ============================================================
//
// dsps_cplx2real_fc32_ae32_() uses the FFT4R table.
//
// For dsps_fft4r_init_fc32(table, 4096), the table requires:
//
//     4 * 4096 floats
//
//     = 16384 floats
//     = 65536 bytes
//     = 64 KB
//
// IMPORTANT:
// This is intentionally based on MAX_FFT_SIZE, not
// MAX_FFT_HALF_SIZE.
// ============================================================

#define FFT4_TABLE_FLOATS  (MAX_FFT_SIZE * 4)
#define FFT4_TABLE_BYTES   (FFT4_TABLE_FLOATS * sizeof(float))

// ============================================================
// Maximum output size
// ============================================================
//
// Real FFT:
//
//     N real samples
//          |
//          v
//     N/2 + 1 complex bins
//
// For N = 4096:
//
//     2049 complex
//     = 4098 floats
// ============================================================

#define FFT_OUTPUT_FLOATS  (MAX_FFT_HALF_SIZE * 2 + 2)


// ============================================================
// FFT TABLES
// ============================================================

static float *fft2_table = NULL;
static float *fft4_table = NULL;

static bool esp_dsp_fftr_initialized = false;


// ============================================================
// INITIALIZATION
// ============================================================
//
// Initializes the tables for the maximum supported FFT size.
//
// After this, esp_dsp_fftr() can be called with:
//
//     nfft = 2048
//
// or:
//
//     nfft = 4096
//
// The same tables are reused.
// ============================================================

bool esp_dsp_fftr_init()
{
    if (esp_dsp_fftr_initialized) {
        return true;
    }

    // --------------------------------------------------------
    // FFT2R table
    // --------------------------------------------------------

    fft2_table = (float *)heap_caps_aligned_alloc(
        16,
        FFT2_TABLE_BYTES,
        MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT
    );

    if (fft2_table == NULL) {
        return false;
    }

    esp_err_t ret =
        dsps_fft2r_init_fc32(
            fft2_table,
            MAX_FFT_HALF_SIZE
        );

    if (ret != ESP_OK) {
        heap_caps_free(fft2_table);
        fft2_table = NULL;
        return false;
    }

    // --------------------------------------------------------
    // FFT4R table
    // --------------------------------------------------------

    fft4_table = (float *)heap_caps_aligned_alloc(
        16,
        FFT4_TABLE_BYTES,
        MALLOC_CAP_INTERNAL | MALLOC_CAP_8BIT
    );

    if (fft4_table == NULL) {
        heap_caps_free(fft2_table);
        fft2_table = NULL;
        return false;
    }

    ret =
        dsps_fft4r_init_fc32(
            fft4_table,
            MAX_FFT_SIZE
        );

    if (ret != ESP_OK) {
        heap_caps_free(fft4_table);
        fft4_table = NULL;

        heap_caps_free(fft2_table);
        fft2_table = NULL;

        return false;
    }

    esp_dsp_fftr_initialized = true;

    return true;
}


// ============================================================
// REAL FFT
// ============================================================
//
// Supported:
//
//     nfft = 2048
//     nfft = 4096
//
// Input:
//
//     timedata[0 ... nfft-1]
//
// Output:
//
//     freqdata[0 ... nfft+1]
//
// as:
//
//     Re(X[0]),   Im(X[0])
//     Re(X[1]),   Im(X[1])
//     ...
//     Re(X[nfft/2]), Im(X[nfft/2])
//
// Therefore:
//
//     nfft = 2048 -> 1025 complex bins -> 2050 floats
//     nfft = 4096 -> 2049 complex bins -> 4098 floats
//
// The caller must provide a sufficiently large freqdata buffer.
// ============================================================

void esp_dsp_fftr(
    void *cfg,
    const float *timedata,
    float *freqdata,
    int nfft)
{
    (void)cfg;

    if (!esp_dsp_fftr_initialized) {
        return;
    }

    // --------------------------------------------------------
    // Validate FFT size
    // --------------------------------------------------------

    if (nfft != 2048 && nfft != 4096) {
        return;
    }

    const int half = nfft / 2;

    // --------------------------------------------------------
    // Copy real input.
    //
    // The first nfft floats are subsequently interpreted as
    // half complex samples:
    //
    //   Re[0], Im[0], Re[1], Im[1], ...
    // --------------------------------------------------------

    memcpy(
        freqdata,
        timedata,
        nfft * sizeof(float)
    );

    // --------------------------------------------------------
    // Complex FFT
    //
    // nfft real samples -> nfft/2 complex samples
    // --------------------------------------------------------

    /*
     * Select the implementation for the target.  In particular, the ESP32-S3
     * macro resolves to the AES3 SIMD kernel; calling the AE32 symbol directly
     * bypasses that faster implementation.
     */
    dsps_fft2r_fc32(
        freqdata,
        half
    );

    // --------------------------------------------------------
    // Bit reversal
    // --------------------------------------------------------

    dsps_bit_rev2r_fc32(
        freqdata,
        half
    );

    // --------------------------------------------------------
    // Complex -> real reconstruction
    //
    // The FFT4 table was initialized for MAX_FFT_SIZE.
    // Keep the table_size corresponding to that table.
    // --------------------------------------------------------

    const int fft4_table_size = MAX_FFT_SIZE * 2;

    dsps_cplx2real_fc32_ae32_(
        freqdata,
        half,
        fft4_table,
        fft4_table_size
    );
}
