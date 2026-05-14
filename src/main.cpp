// ft8_lib cross-platform entrypoint
// - ESP32 Arduino: scans LittleFS and decodes .wav files
// - Linux host: decodes a wav path provided as argv[1]
//
// Flash the filesystem image first:
//   pio run --target uploadfs --environment esp32s3
//
// Then flash the firmware:
//   pio run --target upload --environment esp32s3

#include <time.h>
#include <stdint.h>

#ifdef ARDUINO_ARCH_ESP32

#include <Arduino.h>
#include <LittleFS.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <esp_heap_caps.h>

#include "ft8/decoder_api.h"

static float* alloc_sample_buffer(size_t count)
{
    size_t bytes = count * sizeof(float);
#ifdef ARDUINO_ARCH_ESP32
    // Prefer PSRAM for large FT8 slot buffers; internal DRAM is usually too small.
    if (psramFound()) {
        float* p = (float*)heap_caps_malloc(bytes, MALLOC_CAP_SPIRAM | MALLOC_CAP_8BIT);
        if (p) {
            return p;
        }
    }
#endif
    return (float*)malloc(bytes);
}

// ---------------------------------------------------------------------------
// WAV loader for LittleFS
// Returns heap-allocated float samples (caller must free), or nullptr on error.
// Fills *out_num_samples, *out_num_channels, *out_sample_rate.
// ---------------------------------------------------------------------------
static float* load_wav_littlefs(const char* path,
                                 int* out_num_samples,
                                 int* out_num_channels,
                                 int* out_sample_rate)
{
    File f = LittleFS.open(path, "r");
    if (!f) {
        Serial.printf("[ft8] Cannot open %s\n", path);
        return nullptr;
    }

    // --- Minimal RIFF/WAV header parse ---
    uint8_t hdr[44];
    if (f.read(hdr, 44) != 44) { f.close(); return nullptr; }

    // "RIFF" at 0, "WAVE" at 8, "fmt " at 12
    if (memcmp(hdr, "RIFF", 4) != 0 || memcmp(hdr+8, "WAVE", 4) != 0 ||
        memcmp(hdr+12, "fmt ", 4) != 0) {
        Serial.printf("[ft8] %s: not a WAV file\n", path);
        f.close();
        return nullptr;
    }

    uint16_t audio_fmt   = hdr[20] | (hdr[21] << 8);  // 1 = PCM
    uint16_t num_ch      = hdr[22] | (hdr[23] << 8);
    uint32_t sample_rate = hdr[24] | (hdr[25]<<8) | (hdr[26]<<16) | (hdr[27]<<24);
    uint16_t bits        = hdr[34] | (hdr[35] << 8);

    if (audio_fmt != 1 || bits != 16) {
        Serial.printf("[ft8] %s: only 16-bit PCM WAV supported (fmt=%u bits=%u)\n",
                      path, audio_fmt, bits);
        f.close();
        return nullptr;
    }

    // Skip to "data" chunk (it may not be exactly at offset 36)
    uint32_t data_size = 0;
    {
        uint8_t chunk[8];
        // We already read 44 bytes; fmt chunk size is at hdr[16..19]
        uint32_t fmt_size = hdr[16]|(hdr[17]<<8)|(hdr[18]<<16)|(hdr[19]<<24);
        // Seek to end of fmt chunk: 12 (RIFF+size+WAVE) + 8 (fmt id+size) + fmt_size
        size_t pos = 12 + 8 + fmt_size;
        f.seek(pos);
        while (f.available() >= 8) {
            f.read(chunk, 8);
            uint32_t csz = chunk[4]|(chunk[5]<<8)|(chunk[6]<<16)|(chunk[7]<<24);
            if (memcmp(chunk, "data", 4) == 0) { data_size = csz; break; }
            f.seek(f.position() + csz);
        }
    }

    if (data_size == 0) {
        Serial.printf("[ft8] %s: no data chunk\n", path);
        f.close();
        return nullptr;
    }

    int total_samples = data_size / (bits/8);   // interleaved across channels
    int mono_samples  = total_samples / num_ch;

    float* buf = alloc_sample_buffer((size_t)mono_samples);
    if (!buf) {
        Serial.printf("[ft8] OOM: %d samples (%lu bytes)\n", mono_samples,
                      (unsigned long)(mono_samples * (int)sizeof(float)));
        f.close();
        return nullptr;
    }

    // Read 16-bit PCM, mix down to mono float
    const float scale = 1.0f / 32768.0f;
    int16_t frame[8];  // up to 8 channels, more than enough
    for (int i = 0; i < mono_samples; i++) {
        f.read((uint8_t*)frame, num_ch * 2);
        float s = 0;
        for (int ch = 0; ch < (int)num_ch; ch++) s += frame[ch];
        buf[i] = s / num_ch * scale;
    }

    f.close();
    *out_num_samples  = mono_samples;
    *out_num_channels = num_ch;
    *out_sample_rate  = (int)sample_rate;
    return buf;
}

// ---------------------------------------------------------------------------
// Decode one WAV file
// ---------------------------------------------------------------------------
static void decode_file(const char* path, float base_freq_mhz, bool is_ft8)
{
    Serial.printf("\n[ft8] Decoding %s (base %.3f MHz, %s)\n",
                  path, base_freq_mhz, is_ft8 ? "FT8" : "FT4");

    int num_samples = 0, num_channels = 0, sample_rate = 0;
    float* signal = load_wav_littlefs(path, &num_samples, &num_channels, &sample_rate);
    if (!signal) return;

    Serial.printf("[ft8] Loaded: %d samples @ %d Hz, %d ch, free heap: %lu bytes\n",
                  num_samples, sample_rate, num_channels,
                  (unsigned long)ESP.getFreeHeap());

    ft8_decode_context_t ctx = {};
    ctx.is_ft8         = is_ft8;
    ctx.base_freq_mhz  = base_freq_mhz;
    // ctx.utc and ctx.utc_frac_sec left zero — fill from NTP or RTC if available

    unsigned long t0 = millis();
    int rc = ft8_decode_slot(signal, sample_rate, num_samples, &ctx);
    unsigned long elapsed = millis() - t0;

    free(signal);

    Serial.printf("[ft8] Done in %lu ms, rc=%d, free heap: %lu bytes\n",
                  elapsed, rc, (unsigned long)ESP.getFreeHeap());
}

static void decode_task(void* /*arg*/)
{
    if (!LittleFS.begin(true)) {
        Serial.println("[ft8] LittleFS mount failed - reflash filesystem");
        vTaskDelete(nullptr);
        return;
    }
    Serial.println("[ft8] LittleFS mounted");

    // List available WAV files and decode each one
    File root = LittleFS.open("/");
    File entry;
    while ((entry = root.openNextFile())) {
        String name = String("/") + entry.name();
        entry.close();
        if (name.endsWith(".wav") || name.endsWith(".WAV")) {
            // Default: FT8, 14.074 MHz (20 m). Adjust as needed.
            decode_file(name.c_str(), 14.074f, true);
        }
    }

    Serial.println("\n[ft8] All files decoded.");
    vTaskDelete(nullptr);
}

// ---------------------------------------------------------------------------
// Arduino entry points
// ---------------------------------------------------------------------------
void setup()
{
    Serial.begin(115200);
    while (!Serial) delay(10);
    delay(1500);
    Serial.println("\n[ft8] FT8 decoder starting");
    Serial.printf("[ft8] Free heap: %lu bytes\n", (unsigned long)ESP.getFreeHeap());
#ifdef ARDUINO_ARCH_ESP32
    Serial.printf("[ft8] PSRAM found: %s, free PSRAM: %lu bytes\n",
                  psramFound() ? "yes" : "no",
                  (unsigned long)ESP.getFreePsram());
#endif

    const uint32_t decode_stack_bytes = 32768;
    BaseType_t rc = xTaskCreatePinnedToCore(
        decode_task,
        "ft8_decode",
        decode_stack_bytes,
        nullptr,
        1,
        nullptr,
        APP_CPU_NUM
    );
    if (rc != pdPASS) {
        Serial.printf("[ft8] Failed to create decode task (rc=%ld)\n", (long)rc);
    }
}

void loop()
{
    // Nothing to do — single-shot decode in setup()
    delay(10000);
}

#else

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

#include "common/wave.h"
#include "ft8/decoder_api.h"

int main(int argc, char** argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <wavfile> [base_freq_mhz]\n", argv[0]);
        return 2;
    }

    const char* wav_path = argv[1];
    float base_freq_mhz = 14.074f;
    if (argc >= 3) {
        base_freq_mhz = (float)atof(argv[2]);
    }

    int fd = open(wav_path, O_RDONLY);
    if (fd < 0) {
        perror("open");
        return 1;
    }

    float* signal = NULL;
    int num_samples = 0;
    int num_channels = 0;
    int sample_rate = 0;
    if (load_wav(&signal, &num_samples, &num_channels, &sample_rate, wav_path, fd) != 0) {
        close(fd);
        fprintf(stderr, "Failed to load wav file: %s\n", wav_path);
        return 1;
    }

    ft8_decode_context_t ctx = {0};
    ctx.is_ft8 = true;
    ctx.base_freq_mhz = base_freq_mhz;

    struct timespec t0 = {0};
    struct timespec t1 = {0};
    clock_gettime(CLOCK_MONOTONIC, &t0);
    int rc = ft8_decode_slot(signal, sample_rate, num_samples, &ctx);
    clock_gettime(CLOCK_MONOTONIC, &t1);

    long elapsed_ms = (t1.tv_sec - t0.tv_sec) * 1000L +
                      (t1.tv_nsec - t0.tv_nsec) / 1000000L;
    free(signal);

    fprintf(stdout, "decode time: %ld ms\n", elapsed_ms);
    fprintf(stdout, "process_buffer returned %d\n", rc);
    return rc == 0 ? 0 : 1;
}

#endif
