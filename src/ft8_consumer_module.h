#pragma once

#include <stdbool.h>
#include <stdint.h>

#include <freertos/FreeRTOS.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int sample_rate;
    float base_freq_mhz;
    int sample_queue_depth;
    int finalize_queue_depth;
    int append_batch_size;
    int consumer_task_stack;
    int finalize_task_stack;
    UBaseType_t consumer_task_priority;
    UBaseType_t finalize_task_priority;
    BaseType_t consumer_task_core;
    BaseType_t finalize_task_core;
} ft8_consumer_module_config_t;

// Capture one UTC-aligned FT8 slot in PSRAM, then decode it on the finalizer
// task.  This preserves the full waveform required by wide subtraction.
bool ft8_consumer_module_init(const ft8_consumer_module_config_t* cfg);
bool ft8_consumer_module_start(void);
int ft8_consumer_module_enqueue_i16(const int16_t* samples, int count, TickType_t timeout_ticks);

#ifdef __cplusplus
}
#endif
