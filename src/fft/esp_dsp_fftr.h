#ifndef _ESP_DSP_FFTR_H_
#define _ESP_DSP_FFTR_H_    

#if defined (ARDUINO_ARCH_ESP32)
#include "dsps_fft2r.h"
#include "_kiss_fft_guts.h"
#include <esp_heap_caps.h>
#include <string.h>
/**
 * Wrapper esp-dsp con la stessa ABI di kiss_fftr
 * @param cfg      Puntatore di configurazione (può essere NULL se non usato da esp-dsp)
 * @param timedata Array di input contenente esattamente 4096 campioni reali
 * @param freqdata Array di output complesso di destinazione (almeno 2049 elementi)
 */
#include "dsps_fft2r.h"
#include <string.h>

/**
 * Hardware-accelerated Real-to-Complex FFT wrapper using esp-dsp.
 * Perfectly mimics kiss_fftr ABI for the FT8 decoder at 12800 Hz (4096 points).
 *
 * @param cfg      Configuration pointer (Can be NULL, unused by esp-dsp)
 * @param timedata Input array containing exactly 4096 real samples
 * @param freqdata Output destination array (Must be allocated for at least 2049 complex bins)
 */
void esp_dsp_fftr(void* cfg, const kiss_fft_scalar* timedata, kiss_fft_cpx* freqdata) 
{
    #define FFT_SIZE 4096

    static __attribute__((aligned(16))) float* local_cplx_buffer = NULL;

    if (local_cplx_buffer == NULL) {
        // Allocazione protetta sulla RAM interna veloce (32KB)
        local_cplx_buffer = (float*)heap_caps_malloc(FFT_SIZE * 2 * sizeof(float), MALLOC_CAP_INTERNAL);
        if (local_cplx_buffer == NULL) {
            printf("ERROR: Insufficient internal memory for 4096 FFT buffer!\n");
            return;
        }
    }

    // Caricamento sequenziale dei 4096 campioni reali
    for (int i = 0; i < FFT_SIZE; i++) {
        local_cplx_buffer[2 * i]     = timedata[i]; // Parte Reale
        local_cplx_buffer[2 * i + 1] = 0.0f;        // Parte Immaginaria
    }

    // Esecuzione FFT complessa a 4096 punti
    dsps_fft2r_fc32(local_cplx_buffer, FFT_SIZE);
    dsps_bit_rev_fc32(local_cplx_buffer, FFT_SIZE);

    // Copia lineare dei primi 2049 bin complessi (da 0 a 6400 Hz)
    memcpy(freqdata, local_cplx_buffer, (FFT_SIZE / 2 + 1) * sizeof(kiss_fft_cpx));
}


 #endif

#endif // _ESP_DSP_FFTR_H_
