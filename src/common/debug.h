#pragma once

#include <stdio.h>

#define LOG_DEBUG   0
#define LOG_INFO    1
#define LOG_WARN    2
#define LOG_ERROR   3
#define LOG_FATAL   4

#ifndef LOG_LEVEL
#define LOG_LEVEL LOG_INFO
#endif

#ifndef LOG_USE_ANSI_COLORS
#define LOG_USE_ANSI_COLORS 0
#endif

static inline const char* log_level_name(int level)
{
    switch (level) {
    case LOG_DEBUG: return "DBG";
    case LOG_INFO:  return "INF";
    case LOG_WARN:  return "WRN";
    case LOG_ERROR: return "ERR";
    case LOG_FATAL: return "FTL";
    default:        return "UNK";
    }
}

static inline const char* log_level_color(int level)
{
#if LOG_USE_ANSI_COLORS
    switch (level) {
    case LOG_DEBUG: return "\033[90m"; // grey
    case LOG_INFO:  return "\033[36m"; // cyan
    case LOG_WARN:  return "\033[33m"; // yellow
    case LOG_ERROR: return "\033[31m"; // red
    case LOG_FATAL: return "\033[1;31m"; // bright red
    default:        return "\033[35m"; // magenta
    }
#else
    (void)level;
    return "";
#endif
}

static inline const char* log_color_reset(void)
{
#if LOG_USE_ANSI_COLORS
    return "\033[0m";
#else
    return "";
#endif
}

static inline const char* log_decode_color(void)
{
#if LOG_USE_ANSI_COLORS
    return "\033[32m"; // green
#else
    return "";
#endif
}

// Emit each record with one stdio call.  Capture and decode run in separate
// FreeRTOS tasks, so splitting the prefix and payload can interleave lines.
#define LOG(level, format, ...)                                               \
    do {                                                                      \
        if ((level) >= LOG_LEVEL) {                                           \
            fprintf(stderr, "%s[%s]%s " format,                              \
                    log_level_color(level), log_level_name(level),            \
                    log_color_reset(), ##__VA_ARGS__);                        \
        }                                                                     \
    } while (0)

// Decoded FT8 messages are the useful payload of a decode cycle.  Keep their
// INFO-level filtering semantics while making the entire record green in an
// ANSI-capable monitor.
#define LOG_DECODED(format, ...)                                              \
    do {                                                                      \
        if (LOG_INFO >= LOG_LEVEL) {                                          \
            fprintf(stderr, "%s[INF] " format,                              \
                    log_decode_color(), ##__VA_ARGS__);                      \
        }                                                                     \
    } while (0)
