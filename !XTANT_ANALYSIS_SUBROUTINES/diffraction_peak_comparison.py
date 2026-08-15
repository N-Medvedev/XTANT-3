"""
diffraction_peak_comparison.py

Compares the "ideal" Miller-index diffraction peaks (OUTPUT_diffraction_peaks.dat)
against the full simulated powder pattern (OUTPUT_diffraction_powder.dat), by
searching a small angular window around each nominal 2theta position in the
powder pattern to locate where the peak actually sits at every time step
(it may drift/broaden away from the ideal position as the simulation evolves,
e.g. on melting).

Inputs
------
OUTPUT_diffraction_peaks.dat
    4 comment lines:
      #Time      (hkl1) (hkl2) ...
      #[fs]      [arb.units]
      #2theta [deg]: (val1) (val2) ...
      #q [1/A]:      (val1) (val2) ...
    Then rows: time, intensity_hkl1, intensity_hkl2, ...

OUTPUT_diffraction_powder.dat
    Block format, one block per time step:
      #    <time>
      <angle>  <Total>  <In-In>  <In-O>  <O-O>
      ... (angle = 1 .. 180 deg)
      <blank line(s)>
      #    <next time>
      ...

Outputs
-------
Two PNGs:
  1. Intensity comparison: nominal peak-file intensity vs. the intensity
     found (via windowed local-max search) in the powder pattern, per peak,
     vs. time. Each curve normalized to its own max so the two data sources
     (arbitrary units vs. raw powder counts) can be compared by shape.
  2. Peak-position comparison: nominal 2theta (constant, dashed) vs. the
     angular position of the local maximum found in the powder pattern,
     per peak, vs. time - to see peak drift/broadening.

Usage
-----
Run the script directly; all parameters are in the CONFIG block below.
"""

from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt

# ============================== CONFIG ==============================
SCRIPT_DIR = Path(__file__).resolve().parent

PEAKS_FILE = SCRIPT_DIR / "OUTPUT_diffraction_peaks.dat"
POWDER_FILE = SCRIPT_DIR / "OUTPUT_diffraction_powder.dat"

OUTPUT_PNG_INTENSITY = SCRIPT_DIR / "diffraction_intensity_comparison.png"
OUTPUT_PNG_POSITION = SCRIPT_DIR / "diffraction_peak_positions.png"

# Half-width [deg] of the search window around each nominal 2theta position,
# used to locate the actual peak in the powder pattern.
SEARCH_WINDOW_DEG = 2.0

# Which powder column to use for peak identification / comparison.
# One of: 'Total', 'In-In', 'In-O', 'O-O'
POWDER_COLUMN = "Total"

# Time tolerance [fs] used when matching a powder-file block to a peaks-file row.
TIME_MATCH_TOL_FS = 1e-3

NORMALIZE_INTENSITY_CURVES = True
# ======================================================================


def read_peaks_file(path):
    """
    Parse OUTPUT_diffraction_peaks.dat.

    Returns
    -------
    names   : list[str]     Miller-index labels, e.g. '(211)'
    angles  : ndarray       nominal 2theta [deg], one per peak
    times   : ndarray       time grid [fs]
    inten   : ndarray       shape (ntimes, npeaks), intensities [arb.units]
    """
    with open(path, "r") as f:
        header_lines = [f.readline() for _ in range(4)]

    name_line = header_lines[0]
    angle_line = header_lines[2]

    names = re.findall(r"\(([^)]+)\)", name_line)
    names = [f"({n})" for n in names]

    angle_strs = re.findall(r"\(([^)]+)\)", angle_line)
    angles = np.array([float(a) for a in angle_strs])

    if len(names) != len(angles):
        raise ValueError(
            f"Mismatch between number of peak names ({len(names)}) and "
            f"number of 2theta values ({len(angles)}) parsed from the header."
        )

    data = np.loadtxt(path, skiprows=4)
    times = data[:, 0]
    inten = data[:, 1:1 + len(names)]

    return names, angles, times, inten


def read_powder_file(path):
    """
    Parse the block-structured powder file.

    Each block is preceded by one or more comment lines ('#...'). Only the
    ones that parse as a single float (the time stamp, e.g. '#    -300.0')
    mark the start of a new block; other comment lines (e.g. the repeated
    '#Angle  Total  In-In  In-O  O-O' column header) are simply skipped.

    Returns
    -------
    times  : ndarray               time [fs], one per block
    blocks : list[ndarray]         each shape (nangles, 5):
                                    columns = Angle, Total, In-In, In-O, O-O
    """
    times = []
    blocks = []
    current_time = None
    current_rows = []

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if s == "":
                continue
            if s.startswith("#"):
                payload = s.lstrip("#").strip()
                try:
                    t = float(payload)
                except ValueError:
                    # Not a time stamp (e.g. the column-header comment line) - skip it.
                    continue
                if current_time is not None and current_rows:
                    times.append(current_time)
                    blocks.append(np.array(current_rows, dtype=float))
                current_time = t
                current_rows = []
            else:
                current_rows.append([float(x) for x in s.split()])

        if current_time is not None and current_rows:
            times.append(current_time)
            blocks.append(np.array(current_rows, dtype=float))

    return np.array(times), blocks


POWDER_COLUMNS = ["Angle", "Total", "In-In", "In-O", "O-O"]


def find_peak_in_window(block, angle_nominal, window_deg, col_index):
    """
    Within angle_nominal +/- window_deg, find the local maximum of the
    chosen intensity column in one powder block.

    Returns
    -------
    pos  : float   angle [deg] of the local maximum
    val  : float   intensity value at that maximum
    """
    ang = block[:, 0]
    mask = (ang >= angle_nominal - window_deg) & (ang <= angle_nominal + window_deg)
    if not np.any(mask):
        # Nominal angle is outside the tabulated range - fall back to nearest point.
        idx = np.argmin(np.abs(ang - angle_nominal))
        return ang[idx], block[idx, col_index]

    sub_ang = ang[mask]
    sub_val = block[mask, col_index]
    imax = np.argmax(sub_val)
    return sub_ang[imax], sub_val[imax]


def main():
    names, nominal_angles, peak_times, peak_inten = read_peaks_file(PEAKS_FILE)
    powder_times, powder_blocks = read_powder_file(POWDER_FILE)
    col_index = POWDER_COLUMNS.index(POWDER_COLUMN)

    npeaks = len(names)
    print(f"Peaks file: {len(peak_times)} time steps, {npeaks} peaks: {', '.join(names)}")
    print(f"Powder file: {len(powder_times)} time steps.")

    # Match each powder block to a row in the peaks file by nearest time.
    matched_peak_idx = []
    matched_powder_idx = []
    for ip, t in enumerate(powder_times):
        j = np.argmin(np.abs(peak_times - t))
        if abs(peak_times[j] - t) <= TIME_MATCH_TOL_FS:
            matched_peak_idx.append(j)
            matched_powder_idx.append(ip)
    if len(matched_powder_idx) == 0:
        raise RuntimeError(
            "No time steps could be matched between the two files - "
            "check TIME_MATCH_TOL_FS and that both files come from the same run."
        )
    matched_peak_idx = np.array(matched_peak_idx)
    matched_powder_idx = np.array(matched_powder_idx)
    common_times = peak_times[matched_peak_idx]
    print(f"Matched {len(common_times)} common time steps.")

    # Extract, per peak, the powder-identified position & intensity at every matched time.
    powder_positions = np.empty((len(common_times), npeaks))
    powder_intensities = np.empty((len(common_times), npeaks))
    for k, ip in enumerate(matched_powder_idx):
        block = powder_blocks[ip]
        for p in range(npeaks):
            pos, val = find_peak_in_window(block, nominal_angles[p], SEARCH_WINDOW_DEG, col_index)
            powder_positions[k, p] = pos
            powder_intensities[k, p] = val

    peaks_intensities_matched = peak_inten[matched_peak_idx, :]

    # ---------------- Plot 1: intensity comparison ----------------
    ncols = 3
    nrows = int(np.ceil(npeaks / ncols))
    fig1, axes1 = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    for p in range(npeaks):
        ax = axes1[p // ncols][p % ncols]
        y_peaks = peaks_intensities_matched[:, p].astype(float)
        y_powder = powder_intensities[:, p].astype(float)
        if NORMALIZE_INTENSITY_CURVES:
            if y_peaks.max() > 0:
                y_peaks = y_peaks / y_peaks.max()
            if y_powder.max() > 0:
                y_powder = y_powder / y_powder.max()
        ax.plot(common_times, y_peaks, label="peaks file", lw=1.5)
        ax.plot(common_times, y_powder, label=f"powder ({POWDER_COLUMN})", lw=1.5, ls="--")
        ax.set_title(f"{names[p]}  (2\u03b8={nominal_angles[p]:.2f}\u00b0)", fontsize=10)
        ax.set_xlabel("Time [fs]")
        ax.set_ylabel("Intensity" + (" (norm.)" if NORMALIZE_INTENSITY_CURVES else ""))
        ax.legend(fontsize=7)
    for p in range(npeaks, nrows * ncols):
        axes1[p // ncols][p % ncols].axis("off")
    fig1.suptitle("Diffraction peak intensity: peaks-file vs. powder pattern")
    fig1.tight_layout()
    fig1.savefig(OUTPUT_PNG_INTENSITY, dpi=200)
    print(f"Saved: {OUTPUT_PNG_INTENSITY}")

    # ---------------- Plot 2: peak position comparison ----------------
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    for p in range(npeaks):
        ax = axes2[p // ncols][p % ncols]
        ax.axhline(nominal_angles[p], color="k", ls="--", lw=1, label="nominal (hkl)")
        ax.plot(common_times, powder_positions[:, p], color="tab:red", lw=1.5,
                label="powder local max")
        ax.set_title(f"{names[p]}", fontsize=10)
        ax.set_xlabel("Time [fs]")
        ax.set_ylabel("2\u03b8 [deg]")
        ax.legend(fontsize=7)
    for p in range(npeaks, nrows * ncols):
        axes2[p // ncols][p % ncols].axis("off")
    fig2.suptitle(f"Diffraction peak position drift (search window \u00b1{SEARCH_WINDOW_DEG}\u00b0)")
    fig2.tight_layout()
    fig2.savefig(OUTPUT_PNG_POSITION, dpi=200)
    print(f"Saved: {OUTPUT_PNG_POSITION}")


if __name__ == "__main__":
    main()
