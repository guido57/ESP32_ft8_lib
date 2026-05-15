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

#if defined(ARDUINO_ARCH_ESP32)
#include <Arduino.h>
#include <LittleFS.h>
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <esp_heap_caps.h>
#else
#include <cstdlib>   // malloc, free
#include <cstring>   // memcmp
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#include "ft8/decoder_api.h"

#if defined(ARDUINO_ARCH_ESP32)

// ============================================================================
// ESP32 Arduino Build
// ============================================================================

static float* alloc_sample_buffer(size_t count)
{
    size_t bytes = count * sizeof(float);
    // Prefer PSRAM for large FT8 slot buffers; internal DRAM is usually too small.
    if (psramFound()) {
        float* p = (float*)heap_caps_malloc(bytes, MALLOC_CAP_SPIRAM | MALLOC_CAP_8BIT);
        if (p) {
            return p;
        }
    }
    return (float*)malloc(bytes);
}

// WAV loader for LittleFS
// Returns heap-allocated float samples (caller must free), or nullptr on error.
static float* load_wav_littlefs(const char* path, int* out_num_samples, int* out_num_channels, int* out_sample_rate)
{
    File f = LittleFS.open(path, "r");
    if (!f) {
        Serial.printf("[ft8] Cannot open %s\n", path);
        return nullptr;
    }
    
    // Read and parse WAV header (RIFF container)
    uint8_t riff_header[12];
    if (f.read(riff_header, 12) != 12 || memcmp(riff_header, "RIFF", 4) != 0 || memcmp(riff_header + 8, "WAVE", 4) != 0) {
        Serial.println("[ft8] Not a valid WAVE file");
        f.close();
        return nullptr;
    }
    
    // Read chunks until we find "fmt " and "data"
    uint8_t chunk_header[8];
    uint32_t chunk_size;
    uint16_t audio_format, num_channels, bits_per_sample;
    uint32_t sample_rate;
    bool found_fmt = false, found_data = false;
    
    while (f.read(chunk_header, 8) == 8) {
        chunk_size = (chunk_header[7] << 24) | (chunk_header[6] << 16) | (chunk_header[5] << 8) | chunk_header[4];
        
        if (memcmp(chunk_header, "fmt ", 4) == 0) {
            uint8_t fmt_data[16];
            if (f.read(fmt_data, 16) != 16) break;
            audio_format = (fmt_data[1] << 8) | fmt_data[0];
            num_channels = (fmt_data[3] << 8) | fmt_data[2];
            sample_rate = (fmt_data[7] << 24) | (fmt_data[6] << 16) | (fmt_data[5] << 8) | fmt_data[4];
            bits_per_sample = (fmt_data[15] << 8) | fmt_data[14];
            found_fmt = true;
            if (chunk_size > 16) f.seek(f.position() + chunk_size - 16);
        } else if (memcmp(chunk_header, "data", 4) == 0) {
            found_data = true;
            break;
        } else {
            f.seek(f.position() + chunk_size);
        }
    }
    
    if (!found_fmt || !found_data) {
        Serial.println("[ft8] Missing fmt or data chunk");
        f.close();
        return nullptr;
    }
    
    if (audio_format != 1 || bits_per_sample != 16) {
        Serial.printf("[ft8] Unsupported format: %u-bit, format %u\n", bits_per_sample, audio_format);
        f.close();
        return nullptr;
    }
    
    uint32_t num_samples = chunk_size / (num_channels * (bits_per_sample / 8));
    float* samples = (float*)malloc(num_samples * sizeof(float));
    if (!samples) {
        Serial.println("[ft8] Failed to allocate sample buffer");
        f.close();
        return nullptr;
    }
    
    int16_t raw_sample;
    for (uint32_t i = 0; i < num_samples; i++) {
        if (f.read((uint8_t*)&raw_sample, 2) != 2) {
            Serial.println("[ft8] Unexpected EOF reading samples");
            free(samples);
            f.close();
            return nullptr;
        }
        samples[i] = (float)raw_sample / 32768.0f;
        
        // Skip remaining channels if multichannel
        for (int ch = 1; ch < num_channels; ch++) {
            if (f.read((uint8_t*)&raw_sample, 2) != 2) {
                Serial.println("[ft8] Unexpected EOF reading channels");
                free(samples);
                f.close();
                return nullptr;
            }
        }
    }
    
    f.close();
    *out_num_samples = num_samples;
    *out_num_channels = num_channels;
    *out_sample_rate = sample_rate;
    return samples;
}

static void decode_file(const char* path, float base_freq_mhz, bool is_ft8)
{
    Serial.printf("\n[ft8] Decoding %s (base %.3f MHz, %s)\n", path, base_freq_mhz, is_ft8 ? "FT8" : "FT4");
    int num_samples = 0, num_channels = 0, sample_rate = 0;
    float* signal = load_wav_littlefs(path, &num_samples, &num_channels, &sample_rate);
    if (!signal) return;
    
    Serial.printf("[ft8] Loaded: %d samples @ %d Hz, %d ch, free heap: %lu bytes\n", 
                  num_samples, sample_rate, num_channels, (unsigned long)ESP.getFreeHeap());
    
    ft8_decode_context_t ctx = {};
    ctx.is_ft8 = is_ft8;
    ctx.base_freq_mhz = base_freq_mhz;
    
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
    File root = LittleFS.open("/");
    File entry;
    while ((entry = root.openNextFile())) {
        String name = String("/") + entry.name();
        entry.close();
        if (name.endsWith(".wav") || name.endsWith(".WAV")) {
            decode_file(name.c_str(), 14.074f, true);
        }
    }
    
    Serial.println("\n[ft8] All files decoded.");
    vTaskDelete(nullptr);
}

void setup()
{
    Serial.begin(115200);
    while (!Serial) delay(10);
    delay(1500);
    
    Serial.println("\n[ft8] FT8 decoder starting");
    Serial.printf("[ft8] Free heap: %lu bytes\n", (unsigned long)ESP.getFreeHeap());
    Serial.printf("[ft8] PSRAM found: %s, free PSRAM: %lu bytes\n", 
                  psramFound() ? "yes" : "no", (unsigned long)ESP.getFreePsram());
    
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
    delay(10000);
}

#else

// ============================================================================
// Native Linux Build
// ============================================================================

#include "common/wave.h"

static float* alloc_sample_buffer(size_t count)
{
    size_t bytes = count * sizeof(float);
    return (float*)malloc(bytes);
}

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
    
    long elapsed_ms = (t1.tv_sec - t0.tv_sec) * 1000L + (t1.tv_nsec - t0.tv_nsec) / 1000000L;
    
    free(signal);
    close(fd);
    
    fprintf(stdout, "decode time: %ld ms\n", elapsed_ms);
    fprintf(stdout, "process_buffer returned %d\n", rc);
    
    return rc == 0 ? 0 : 1;
}

#endif
