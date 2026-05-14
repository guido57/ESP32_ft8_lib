"""
PlatformIO pre-script: utils/pio_littlefs.py
Stages WAV files into data/ for LittleFS packaging.

Behavior:
- If FT8_WAV_FILE is set: copy that WAV into data/ (overwrite if exists).
- Else if data/ already contains .wav files: keep them untouched.
- Else auto-stage one WAV (prefer CQIW5ALZ.wav, otherwise first tests/*.wav).
"""

Import("env")
import os
import shutil

project_dir = env.subst("$PROJECT_DIR")
data_dir = os.path.join(project_dir, "data")
tests_dir = os.path.join(project_dir, "tests")

os.makedirs(data_dir, exist_ok=True)

selected = None
override = os.environ.get("FT8_WAV_FILE")
if override:
    selected = override if os.path.isabs(override) else os.path.join(project_dir, override)
else:
    existing_wavs = [
        name for name in os.listdir(data_dir)
        if name.lower().endswith(".wav")
    ]
    if existing_wavs:
        print("[littlefs] keeping existing data/*.wav:", ", ".join(sorted(existing_wavs)))
    else:
        preferred = os.path.join(project_dir, "CQIW5ALZ.wav")
        if os.path.exists(preferred):
            selected = preferred
        else:
            candidates = sorted(
                os.path.join(tests_dir, name)
                for name in os.listdir(tests_dir)
                if name.lower().endswith(".wav")
            )
            if candidates:
                selected = candidates[0]

if selected and os.path.exists(selected):
    dst = os.path.join(data_dir, os.path.basename(selected))
    shutil.copy2(selected, dst)
    print(f"[littlefs] staged {os.path.basename(selected)}")
elif selected:
    print(f"[littlefs] selected WAV not found: {selected}")
else:
    print("[littlefs] no WAV selected/staged")
