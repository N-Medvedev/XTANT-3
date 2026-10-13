import os
import sys
import re
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Find the photon data file
# ============================================================

photon_files = sorted(glob.glob("OUTPUT_*_Photon*.dat"))

if len(photon_files) == 0:
    print("Error: No file matching 'OUTPUT_*_Photon*.dat' found.")
    sys.exit(1)

if len(photon_files) > 1:
    print("Error: Multiple photon data files found:")
    for f in photon_files:
        print(f"   {f}")
    sys.exit(1)

data_file = photon_files[0]

print(f"Using data file: {data_file}")

# ============================================================
# Determine matching parameter file
# ============================================================

match = re.match(
    r"OUTPUT_(.+?)_Photon.*\.dat$",
    os.path.basename(data_file)
)

if not match:
    print(f"Error: Cannot determine material from '{data_file}'.")
    sys.exit(1)

material = match.group(1)

param_file = f"!OUTPUT_{material}_Parameters.txt"

if not os.path.exists(param_file):
    print(f"Error: Matching parameter file not found:")
    print(f"       {param_file}")
    sys.exit(1)

print(f"Using parameter file: {param_file}")

# ============================================================
# Read atomic density
# ============================================================

density = None

with open(param_file, "r") as f:
    for line in f:
        match = re.search(
            r"atomic density.*?:\s*([0-9.]+[Ee][+-]?[0-9]+)",
            line,
            re.IGNORECASE,
        )

        if match:
            density = float(match.group(1))
            break

if density is None:
    print(f"Error: Atomic density not found in '{param_file}'.")
    sys.exit(1)

print(f"Atomic density = {density:.6E} [1/cm^3]")

# ============================================================
# Read attenuation data
# ============================================================

df = pd.read_csv(
    data_file,
    sep=r"\s+",
    header=None,
    comment="#",
    skipinitialspace=True,
    skip_blank_lines=True,
)

if df.shape[1] < 2:
    print(f"Error: '{data_file}' contains insufficient columns.")
    sys.exit(1)

print(
    f"Read {df.shape[0]} rows and {df.shape[1]} columns. "
    f"Using last column ({df.shape[1]-1}) as attenuation length."
    f"Insert correct damage threshold(s) into the script!!"
)

# ============================================================
# Extract data
# ============================================================

energy = df.iloc[:, 0]
attenuation_length = df.iloc[:, -1]  # last column

# ============================================================
# Physical constants and conversion factors
# ============================================================

coef = density        # [1/cm^3]
coef_e = 1.6e-19      # eV -> J

factor = (
    np.array([4.0, 2.0, 1.7, 1.1])      # INSERT CORRECT DAMAGE THRESHOLD(S) FOR THIS MATERIAL !!!
    * coef
    * 1.0e-8
    * coef_e
)

labels = ["NVE BO", "NPH BO", "NVE, DC", "NPH DC"]

# ============================================================
# Create figure
# ============================================================

plt.figure(figsize=(5, 4))

plt.xlabel("Photon energy (eV)", fontsize=14)
plt.ylabel("Threshold fluence (J/cm$^2$)", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.xscale("log")
plt.yscale("log")

plt.grid(False)

colors = plt.cm.tab10(np.linspace(1, 0, 10))
linestyles = ["-", "-.", "--", ":"]

# ============================================================
# Plot curves
# ============================================================

for i in range(len(factor)):
    plt.plot(
        energy,
        attenuation_length * factor[i],
        color=colors[i],
        linestyle=linestyles[i],
        label=labels[i],
        linewidth=2,
    )

# ============================================================
# Axes limits
# ============================================================

xmax = energy.max()

plt.xlim(10.0, min(xmax, 50000.0))
plt.ylim(1.0e-2, 1.0e4)

plt.title("Damage threshold fluence")
plt.legend(loc="best", fontsize=12)

# ============================================================
# Save figure
# ============================================================

output_file = "OUTPUT_threshold_fluence.png"

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

print(f"Plot saved: {output_file}")