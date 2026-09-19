#ifndef FT8_NATIVE_OSD_H
#define FT8_NATIVE_OSD_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Native-only ordered-statistics list decoder.  A successful return provides
// a complete, re-encoded 174-bit LDPC codeword.  The caller must still apply
// the normal FT8 CRC and message-unpacking checks.
int native_osd_recover(const float log174[], uint8_t recovered174[], float* score);

#ifdef __cplusplus
}
#endif

#endif
