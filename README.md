# FT8/FT4 Decoder - ESP32 + Linux

Upstream reference (original work): Karlis Goba (YL3JG), repository: https://github.com/kgoba/ft8_lib

This repository is a PlatformIO-focused adaptation of that FT8/FT4 decoder, set up to run with one shared application entrypoint on:

- ESP32-S3 (Arduino framework)
- Linux host (PlatformIO native environment)

## What this repo is for

- Decode FT8/FT4 from WAV files on ESP32-S3 using LittleFS storage.
- Run the same decode application logic on Linux for fast testing.
- Keep the decoder core in C (`src/decode_ft8.c`, `ft8/`, `fft/`, `common/`) and use thin platform-specific glue in `src/main.cpp`.

## Project layout (important parts)

- `src/main.cpp`: Cross-platform app entrypoint
	- ESP32 path: scans LittleFS and decodes WAV files
	- Linux path: decodes the WAV path passed on command line
- `src/decode_ft8.c`: Main decoder processing implementation (`process_buffer`)
- `ft8/`, `fft/`, `common/`: Decoder/FFT/support code
- `tests/`: WAV test assets (kept intentionally)
- `data/`: Files packaged into LittleFS image for ESP32
- `platformio.ini`: Build environments and flags
- `utils/pio_sources.py`: Adds non-`src/` core sources to PlatformIO build
- `utils/pio_littlefs.py`: Controls WAV staging behavior for LittleFS image

## Build and run

### Linux host (native)

Build:

```bash
pio run --environment native
```

Run:

```bash
./.pio/build/native/program tests/CQIW5ALZ.wav 14.074
```

The output includes decoded messages plus decode execution time.

### ESP32-S3 (Arduino)

Build firmware:

```bash
pio run --environment esp32s3
```

Upload LittleFS image (`data/` contents):

```bash
pio run --environment esp32s3 --target uploadfs --upload-port /dev/ttyACM1
```

Upload firmware:

```bash
pio run --environment esp32s3 --target upload --upload-port /dev/ttyACM1
```

Monitor serial output:

```bash
pio device monitor --port /dev/ttyACM1 --baud 115200
```

## LittleFS WAV staging behavior

`utils/pio_littlefs.py` works as follows:

- If `FT8_WAV_FILE` is set, that WAV is copied into `data/`.
- Else if `data/` already contains `.wav`, they are kept as-is.
- Else one WAV is auto-staged from available test assets.

Example forcing a specific file:

```bash
FT8_WAV_FILE=tests/191111_110130.wav pio run --environment esp32s3 --target uploadfs --upload-port /dev/ttyACM1
```

## Notes

- Decoding is demonstrated end-to-end on ESP32-S3 and Linux.
- If output timestamps show `1900/...`, UTC metadata is not yet being populated in the decode context.

## Credits

- Karlis Goba (YL3JG) for the original `ft8_lib` project and core implementation direction.
- Mark Borgerding for KissFFT: https://github.com/mborgerding/kissfft
- WSJT-X authors for FT8/FT4 protocol ecosystem.
