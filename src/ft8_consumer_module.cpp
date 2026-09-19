#include "ft8_consumer_module.h"

#include <atomic>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include <Arduino.h>
#include <freertos/queue.h>
#include <freertos/task.h>

#include "decoder_api.h"
#include "common/debug.h"

typedef struct {
    ft8_stream_decoder_t* stream;
    TickType_t queued_at;
    bool is_early_pass;
} finalize_job_t;

static QueueHandle_t s_sample_queue = nullptr;
static QueueHandle_t s_finalize_queue = nullptr;
static ft8_consumer_module_config_t s_cfg = {};
static std::atomic<uint32_t> s_input_samples{0};
static std::atomic<uint32_t> s_dropped_samples{0};
static constexpr float kEarlyPassSeconds = 11.75f;

static double utc_now(void)
{
    struct timeval tv = {};
    gettimeofday(&tv, nullptr);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

static void decoder_consumer_task(void*)
{
    const int slot_samples = (int)(15.0 * s_cfg.sample_rate + 0.5);
    const int early_samples =
        (int)lroundf(kEarlyPassSeconds * s_cfg.sample_rate);
    int16_t* batch = static_cast<int16_t*>(malloc((size_t)s_cfg.append_batch_size * sizeof(*batch)));
    if (batch == nullptr) {
        LOG(LOG_ERROR, "[ft8] slot capture batch allocation failed\n");
        vTaskDelete(nullptr);
        return;
    }

    double next_slot = ceil(utc_now() / 15.0) * 15.0;
    for (;;) {
        while (utc_now() < next_slot)
            vTaskDelay(pdMS_TO_TICKS(2));

        // Do not let samples queued before the UTC boundary contaminate this
        // slot. The producer is non-blocking, so a full queue simply drops.
        xQueueReset(s_sample_queue);
        const uint32_t input_at_start =
            s_input_samples.load(std::memory_order_relaxed);
        const uint32_t dropped_at_start =
            s_dropped_samples.load(std::memory_order_relaxed);

        time_t slot_seconds = (time_t)next_slot;
        ft8_decode_context_t ctx = {};
        ctx.is_ft8 = true;
        ctx.base_freq_mhz = s_cfg.base_freq_mhz;
        gmtime_r(&slot_seconds, &ctx.utc);
        ctx.utc_frac_sec = next_slot - (double)slot_seconds;
        ctx.cand_to_subtract = -1;

        ft8_stream_decoder_t* stream = ft8_stream_open(s_cfg.sample_rate, &ctx);
        if (stream == nullptr) {
            LOG(LOG_ERROR, "[ft8] slot capture allocation failed\n");
            next_slot += 15.0;
            continue;
        }

        const double slot_end = next_slot + 15.0;
        int captured = 0;
        bool early_attempted = false;
        while (utc_now() < slot_end && captured < slot_samples) {
            int n = 0;
            TickType_t wait = pdMS_TO_TICKS(5);
            if (xQueueReceive(s_sample_queue, &batch[n], wait) != pdTRUE)
                continue;
            ++n;
            while (n < s_cfg.append_batch_size &&
                   xQueueReceive(s_sample_queue, &batch[n], 0) == pdTRUE)
                ++n;

            const int accepted = ft8_stream_append_i16(stream, batch, n);
            if (accepted < 0) {
                captured = 0;
                break;
            }
            captured += accepted;
            if (!early_attempted && captured >= early_samples) {
                early_attempted = true;
                ft8_stream_decoder_t* early =
                    ft8_stream_snapshot_early(stream);
                if (early == nullptr) {
                    LOG(LOG_WARN, "[ft8] early pass snapshot allocation failed\n");
                } else {
                    finalize_job_t early_job = {
                        early, xTaskGetTickCount(), true
                    };
                    if (xQueueSend(s_finalize_queue, &early_job, 0) == pdTRUE) {
                        LOG(LOG_INFO,"[ft8] early pass queued samples=%d duration=%.3fs\n",
                                    captured,
                                    (double)captured / s_cfg.sample_rate);
                    } else {
                        LOG(LOG_WARN,"[ft8] finalizer busy; early pass dropped\n");
                        ft8_stream_close(early);
                    }
                }
            }
            if ((captured & 0x3ff) == 0)
                taskYIELD();
        }

        if (captured > 0) {
            const uint32_t input_samples =
                s_input_samples.load(std::memory_order_relaxed) - input_at_start;
            const uint32_t dropped_samples =
                s_dropped_samples.load(std::memory_order_relaxed) - dropped_at_start;
            LOG(LOG_INFO,"[ft8] captured samples=%d/%d input=%lu dropped=%lu\n",
                          captured, slot_samples,
                          (unsigned long)input_samples,
                          (unsigned long)dropped_samples);
            finalize_job_t job = {
                stream, xTaskGetTickCount(), false
            };
            if (xQueueSend(s_finalize_queue, &job, 0) != pdTRUE) {
                LOG(LOG_WARN,"[ft8] finalizer busy; slot dropped\n");
                ft8_stream_close(stream);
            }
        } else {
            ft8_stream_close(stream);
        }

        next_slot += 15.0;
        // Finishing a capture a few milliseconds after the UTC boundary is
        // normal: the last queue receive can straddle that boundary.  Do not
        // turn that tiny jitter into a complete 15-second slot skip.  FT8
        // transmissions begin about 0.5 s into a slot, so a prompt capture
        // after such jitter still has the full signal.  Only resynchronise
        // when the capture task was genuinely unable to start the next slot.
        const double late_sec = utc_now() - next_slot;
        if (late_sec > 0.250) {
            LOG(LOG_WARN,"[ft8] capture late=%.0fms; resynchronizing\n",
                          late_sec * 1000.0);
            next_slot = ceil(utc_now() / 15.0) * 15.0;
        }
    }
}

static void finalize_worker_task(void*)
{
    for (;;) {
        finalize_job_t job = {};
        if (xQueueReceive(s_finalize_queue, &job, portMAX_DELAY) != pdTRUE)
            continue;
        const TickType_t waited_ticks = xTaskGetTickCount() - job.queued_at;
        const uint32_t waited_ms = (uint32_t)(waited_ticks * portTICK_PERIOD_MS);
        if (job.is_early_pass) {
            if (waited_ms > 1000) {
                // An early pass is useful only if it starts promptly. The
                // untouched full stream will retain the normal two passes.
                LOG(LOG_WARN,"[ft8] early pass stale wait=%lums; full decode retained\n",
                              (unsigned long)waited_ms);
            } else {
                const int early_decoded = ft8_stream_finalize(job.stream);
                if (early_decoded < 0) {
                    LOG(LOG_ERROR,"[ft8] early pass failed rc=%d; full decode retained\n",
                                  early_decoded);
                } else {
                    LOG(LOG_INFO,"[ft8] early pass complete decoded=%d; full waveform unchanged\n",
                                  early_decoded);
                }
            }
            ft8_stream_close(job.stream);
            continue;
        }

        if (waited_ms > 1000)
            LOG(LOG_WARN,"[ft8] full decode wait=%lums; preserving all configured passes\n",
                          (unsigned long)waited_ms);
        const int rc = ft8_stream_finalize(job.stream);
        if (rc < 0)
            LOG(LOG_ERROR,"[ft8] slot decode failed rc=%d\n", rc);
        ft8_stream_close(job.stream);
    }
}

bool ft8_consumer_module_init(const ft8_consumer_module_config_t* cfg)
{
    if (cfg == nullptr || cfg->sample_rate <= 0 || cfg->sample_queue_depth <= 0 ||
        cfg->finalize_queue_depth <= 0 || cfg->append_batch_size <= 0)
        return false;
    s_cfg = *cfg;
    s_sample_queue = xQueueCreate((UBaseType_t)cfg->sample_queue_depth, sizeof(int16_t));
    s_finalize_queue = xQueueCreate((UBaseType_t)cfg->finalize_queue_depth, sizeof(finalize_job_t));
    return s_sample_queue != nullptr && s_finalize_queue != nullptr;
}

bool ft8_consumer_module_start(void)
{
    if (s_sample_queue == nullptr || s_finalize_queue == nullptr)
        return false;
    return xTaskCreatePinnedToCore(decoder_consumer_task, "ft8_capture",
                                   (uint32_t)s_cfg.consumer_task_stack, nullptr,
                                   s_cfg.consumer_task_priority, nullptr,
                                   s_cfg.consumer_task_core) == pdPASS &&
           xTaskCreatePinnedToCore(finalize_worker_task, "ft8_finalize",
                                   (uint32_t)s_cfg.finalize_task_stack, nullptr,
                                   s_cfg.finalize_task_priority, nullptr,
                                   s_cfg.finalize_task_core) == pdPASS;
}

int ft8_consumer_module_enqueue_i16(const int16_t* samples, int count, TickType_t timeout_ticks)
{
    if (s_sample_queue == nullptr || samples == nullptr || count < 0)
        return -1;
    int i = 0;
    for (; i < count; ++i) {
        s_input_samples.fetch_add(1, std::memory_order_relaxed);
        if (xQueueSend(s_sample_queue, &samples[i], timeout_ticks) != pdTRUE) {
            s_dropped_samples.fetch_add(1, std::memory_order_relaxed);
            break;
        }
    }
    return i;
}
