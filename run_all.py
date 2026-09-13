#!/usr/bin/env python3

import subprocess

PROGRAM = ".pio/build/native/program"
DATA_DIR = "tests/20m_busy"

for i in range(1, 39):
    wav_file = f"{DATA_DIR}/test_{i:02d}.wav"

    cmd = [
        PROGRAM,
        wav_file,
        "-1",
        "0.0",
        "0.0",
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    print(f"\n===== {wav_file} =====")

    for line in result.stdout.splitlines():
        if "Pass 1: decoded" in line:
            print(line)