#!/usr/bin/env python3

import subprocess

program = ".pio/build/native/program"
wav = "tests/20m_busy/test_01.wav"

for i in range(1400, 1601):
    value = i / 1000.0

    print(f"\n{'='*60}", flush=True)
    print(f"Running with parameter: {value:.3f}", flush=True)
    print(f"{'='*60}", flush=True)

    result = subprocess.run(
        [program, wav, f"{value:.3f}"],
        capture_output=True,
        text=True
    )

    print(result.stdout, end="")

    if result.stderr:
        print("STDERR:")
        print(result.stderr, end="")

    print(f"\nReturn code: {result.returncode}", flush=True)

    if result.returncode != 0:
        print(f"*** FAILED at {value:.3f} ***", flush=True)s