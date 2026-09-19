#ifndef _DEBUG_H_INCLUDED_
#define _DEBUG_H_INCLUDED_

#define LOG_DEBUG 0
#define LOG_INFO  1
#define LOG_WARN  2
#define LOG_ERROR 3
#define LOG_FATAL 4

#ifndef LOG_LEVEL
#define LOG_LEVEL LOG_INFO
#endif

#ifndef LOG_USE_ANSI_COLORS
#define LOG_USE_ANSI_COLORS 0
#endif

#ifndef LOG_PRINTF
#include <stdio.h>
#define LOG_PRINTF(...) fprintf(stderr, __VA_ARGS__)
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
    case LOG_DEBUG: return "\033[90m";
    case LOG_INFO:  return "\033[36m";
    case LOG_WARN:  return "\033[33m";
    case LOG_ERROR: return "\033[31m";
    case LOG_FATAL: return "\033[1;31m";
    default:        return "\033[35m";
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

#define LOG(level, format, ...)                                               \
    do {                                                                      \
        if ((level) >= LOG_LEVEL) {                                           \
            LOG_PRINTF("%s[%s]%s " format,                                  \
                       log_level_color(level), log_level_name(level),         \
                       log_color_reset(), ##__VA_ARGS__);                     \
        }                                                                     \
    } while (0)

#endif // _DEBUG_H_INCLUDED_
