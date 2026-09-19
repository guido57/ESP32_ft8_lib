#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Data: percentage relative to WSJT-X = 100%

pyft8_perc = [
    65,70,60,65,70,75,68,91,53,91,
    56,73,64,110,73,107,74,76,64,91,
    58,73,59,84,63,91,66,71,62,68,
    68,70,66,64,59,67,57,67
]

ft8lib_perc = [
    69,70,60,65,67,75,61,91,56,86,
    59,68,68,107,73,107,70,81,64,81,
    58,73,66,84,67,82,63,67,69,68,
    72,73,66,61,59,67,57,62
]

ft8mon_perc = [
    122,114,150,129,90,96,94,100,96,110,
    117,105,114,106,115,107,96,90,100,95,
    103,118,115,113,96,96,100,95,96,100,
    100,96,100,96,103,91,136,105
]


# WSJT-X reference decoded messages

wsjtx = [
    27,22,16,21,30,27,31,17,28,20,
    29,21,29,16,26,15,27,20,32,22,
    32,22,26,23,27,23,30,22,26,29,
    25,27,28,25,30,22,22,21
]


# ft8mon1 decoded messages (absolute values)

ft8mon1 = [
    33,25,25,27,29,27,29,19,27,21,
    33,21,32,15,29,18,26,19,31,23,
    33,24,32,25,29,23,30,23,25,29,
    26,26,29,27,33,23,32,22
]


ft8mon2 = [
    28,25,26,26,28,30,27,18,27,20,
    31,22,29,15,29,17,25,19,31,18,
    31,26,31,24,25,23,30,22,25,29,
    25,28,28,26,31,21,31,22
]

ft8lib_subtract_time_osr_2 = [
    26,20,18,19,25,26,23,21,22,20,
    24,16,25,17,25,17,21,20,26,19,
    26,24,22,20,22,20,25,16,22,22,
    23,26,25,19,24,18,22,15
]

# FT8_lib subtraction, frequency/time OSR = 2, LDPC iterations = 40
ft8lib_subtract_wide = [
    27,20,16,19,24,23,22,21,22,19,
    24,16,25,16,24,18,21,19,25,19,
    24,24,22,20,21,20,24,17,21,21,
    23,23,25,17,22,18,23,16
]

ft8lib_subtract_time_osr_4 = [
    25, 20, 19, 21, 24, 29, 24, 21, 24, 20,
    26, 18, 24, 18, 25, 17, 22, 20, 26, 20,
    26, 24, 21, 20, 22, 23, 25, 21, 23, 22,
    23, 26, 25, 18, 24, 18, 22, 16
]

ft8lib_subtract_DDC = [
    25, 19, 14, 17, 22, 24, 20, 19, 20, 18,
    20, 16, 20, 16, 23, 17, 22, 17, 25, 17,
    25, 21, 17, 21, 20, 20, 22, 18, 19, 19,
    20, 20, 19, 17, 22, 17, 18, 13
]



# Calculate ft8mon1 percentage relative to WSJT-X

ft8mon1_perc = 100 * np.array(ft8mon1) / np.array(wsjtx)
ft8mon2_perc = 100 * np.array(ft8mon2) / np.array(wsjtx)
ft8lib_subtract_time_osr_2_perc = 100 * np.array(ft8lib_subtract_time_osr_2) / np.array(wsjtx)
ft8lib_subtract_wide_perc = 100 * np.array(ft8lib_subtract_wide) / np.array(wsjtx)
ft8lib_subtract_time_osr_4_perc = 100 * np.array(ft8lib_subtract_time_osr_4) / np.array(wsjtx)
ft8lib_subtract_DDC_perc = 100 * np.array(ft8lib_subtract_DDC) / np.array(wsjtx)



# X axis: test files 1...38

files = np.arange(1, 39)


# Bar width

width = 0.11


plt.figure(figsize=(16, 6))


plt.bar(files - 4*width,
        pyft8_perc,
        width,
        label=f"PyFT8 avg={np.mean(pyft8_perc):.1f}%")


plt.bar(files - 3*width,
        ft8lib_perc,
        width,
        label=f"FT8_lib avg={np.mean(ft8lib_perc):.1f}%")

plt.bar(files - 2*width,
        ft8lib_subtract_time_osr_2_perc,
        width,
        label=f"FT8_lib_subtract_time_osr_2 avg={np.mean(ft8lib_subtract_time_osr_2_perc):.1f}%")

plt.bar(files - width,
        ft8lib_subtract_wide_perc,
        width,
        label=f"FT8_lib_subtract_wide (LDPC=40) avg={np.mean(ft8lib_subtract_wide_perc):.1f}%")

plt.bar(files,
        ft8lib_subtract_time_osr_4_perc,
        width,
        label=f"FT8_lib_subtract_time_osr_4 avg={np.mean(ft8lib_subtract_time_osr_4_perc):.1f}%")


plt.bar(files + width,
        ft8lib_subtract_DDC_perc,
        width,
        label=f"FT8_lib_subtract_DDC avg={np.mean(ft8lib_subtract_DDC_perc):.1f}%")

plt.bar(files + 2*width,
        ft8mon_perc,
        width,
        label=f"ft8mon hz_frac_n=4 off_frac_n=4 avg={np.mean(ft8mon_perc):.1f}%")


plt.bar(files + 3*width,
        ft8mon1_perc,
        width,
        label=f"ft8mon1 hz_frac_n=1 off_frac_n=4 avg={np.mean(ft8mon1_perc):.1f}%")

plt.bar(files + 4*width,
        ft8mon2_perc,
        width,
        label=f"ft8mon2 hz_frac_n=1 off_frac_n=1 avg={np.mean(ft8mon2_perc):.1f}%")

# WSJT-X reference line

plt.axhline(100,
            linestyle="--",
            label="WSJT-X = 100%")


plt.xlabel("Test file")
plt.ylabel("Decoded messages relative to WSJT-X (%)")

plt.title("FT8 Decoder Performance per Test File (38 WAV files)")


plt.xticks(files)

plt.grid(axis="y")

plt.legend()


plt.tight_layout()


plt.savefig("decoder_per_file.png", dpi=150)


plt.show()


# Print averages

print()
print("Average performance:")
print(f"PyFT8   : {np.mean(pyft8_perc):.1f}%")
print(f"FT8_lib  : {np.mean(ft8lib_perc):.1f}%")
print(f"FT8_lib DDC: {np.mean(ft8lib_subtract_DDC_perc):.1f}%")
print(f"FT8_lib subtract wide (LDPC=40): {np.mean(ft8lib_subtract_wide_perc):.1f}%")
print(f"ft8mon   : {np.mean(ft8mon_perc):.1f}%")
print(f"ft8mon2  : {np.mean(ft8mon2_perc):.1f}%")
print(f"ft8mon1  : {np.mean(ft8mon1_perc):.1f}%")
