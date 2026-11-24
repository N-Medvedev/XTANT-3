#!/usr/bin/env python3
"""
Go through folders like OUTPUT_MgO_hw=7000.0_t=6.0_F=1.50_-400, read
OUTPUT_diffraction_peaks.dat from each, align the data in time, average the
peak intensities over the folders, and plot them.

The peak list is read from the header line starting with '#Time', e.g.
    #Time          (200)       (220)       (311)       (222)
It may differ between simulations: peaks are matched by their label, e.g.
"(200)", and each peak is averaged only over the folders that contain it.

Alignment: by the absolute Time column (default). Each run is interpolated
onto the common grid of all time points; at each time the average runs only
over the runs that cover that time. With --relative, every run is shifted to
start at t = 0.

Output (in <root>/):
  peaks_average.dat   (Time + one column per peak, with 2theta/q in the header)
  peaks_average.png

Usage:
    python XTANT_peaks_average.py [root_dir] [-p "OUTPUT_MgO_*"] [--relative] [--log]
"""
import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_FILE = "OUTPUT_diffraction_peaks_DW.dat"
LABEL_RE = re.compile(r"\(([^)]+)\)")

parser = argparse.ArgumentParser()
parser.add_argument("root", nargs="?", default=".", help="directory with the folders")
parser.add_argument("-p", "--pattern", default="OUTPUT_*", help="folder name pattern")
parser.add_argument("--relative", action="store_true",
                    help="shift each run so that it starts at t = 0")
parser.add_argument("--log", action="store_true", help="logarithmic y axis")
args = parser.parse_args()

root = Path(args.root)


def read_header(fpath):
    """Return peak labels (from the '#Time' line) and {label: 2theta string}."""
    labels, two_theta = [], {}
    with open(fpath) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.lstrip("#").lstrip().startswith("Time"):
                labels = LABEL_RE.findall(line)
            elif "2theta" in line:
                vals = LABEL_RE.findall(line.split(":", 1)[1])
                two_theta = dict(zip(labels, vals))
    return labels, two_theta


runs = []          # list of (time, {label: intensity array})
two_theta_all = {}  # label -> 2theta string (from the first file with that peak)

for folder in sorted(root.glob(args.pattern)):
    fpath = folder / DATA_FILE
    if not folder.is_dir() or not fpath.is_file():
        continue
    labels, two_theta = read_header(fpath)
    data = np.atleast_2d(np.loadtxt(fpath, comments="#"))
    if not labels or data.shape[1] - 1 < len(labels):
        print(f"WARNING: {fpath}: header/data mismatch, skipping")
        continue
    data = data[np.argsort(data[:, 0])]
    t = data[:, 0] - (data[0, 0] if args.relative else 0.0)
    runs.append((t, {lab: data[:, 1 + i] for i, lab in enumerate(labels)}))
    for lab in labels:
        two_theta_all.setdefault(lab, two_theta.get(lab, "?"))
    print(f"Read {fpath}  ({len(t)} points, peaks: {' '.join(labels)})")

if not runs:
    raise SystemExit("No matching folders / data found.")

grid = np.unique(np.concatenate([r[0] for r in runs]))
all_labels = sorted({lab for r in runs for lab in r[1]},
                    key=lambda s: [int(x) if x.lstrip("-").isdigit() else x
                                   for x in re.split(r"(-?\d+)", s)])

fig, ax = plt.subplots(figsize=(6, 4.5))
columns = []
for lab in all_labels:
    curves = [np.interp(grid, t, peaks[lab], left=np.nan, right=np.nan)
              for t, peaks in runs if lab in peaks]
    avg = np.nanmean(curves, axis=0)
    columns.append(avg)
    ax.plot(grid, avg, lw=1.5, label=f"({lab})")

header = ("Time   " + "   ".join(f"({lab})" for lab in all_labels) + "\n[fs]   [arb.units]\n"
          "2theta [deg]: " + "   ".join(f"({two_theta_all[lab]})" for lab in all_labels))

np.savetxt(root / "OUT_DW_diffraction_peaks_average.dat", np.column_stack([grid] + columns),
           fmt="%.10E", header=header)

ax.set_xlabel("Time (fs)", fontsize=12)
ax.set_ylabel("DW peak intensity (arb. units)", fontsize=12)

ax.set_xlim(-100.0,500.0)
ax.set_ylim(None,None)

ax.set_title(f"DW diffraction peaks (average of {len(runs)} runs)")
if args.log:
    ax.set_yscale("log")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(root / "OUT_DW_diffraction_peaks_average.png", dpi=300, bbox_inches="tight")
print(f"Saved OUT_DW_diffraction_peaks_average.png and OUT_DW_diffraction_peaks_average.dat ({len(runs)} runs, {len(all_labels)} peaks)")
