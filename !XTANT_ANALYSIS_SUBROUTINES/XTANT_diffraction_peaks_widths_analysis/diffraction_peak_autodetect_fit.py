"""
diffraction_peak_fit.py
------------------------
Extract integrated Bragg-peak intensities (and widths) from XTANT-3-style
powder diffraction output, using the known Miller-index peak positions as
fitting guides.

Handles:
  * Multi-block "OUTPUT_diffraction_powder*.dat" files (one block per time
    stamp, columns: 2theta, Total, <element-pair columns...>).
  * "OUTPUT_diffraction_peaks*.dat" header, which gives the nominal 2theta
    (and q) position for each named Miller-index reflection.
  * Fitting each reflection (or group of overlapping reflections) with a
    pseudo-Voigt profile + local linear background, and extracting the
    ANALYTIC integrated area (not a noisy trapezoidal sum).
  * Averaging over several independent runs (different initial velocities)
    either before fitting (average raw curves -> one fit) or after fitting
    (fit each run -> average the fitted areas, get a run-to-run std as an
    error bar).

This is meant to be read top to bottom and adapted -- it is not a black box.
"""

from __future__ import annotations

import os
import re
import glob
from dataclasses import dataclass, field
from typing import Sequence
from datetime import datetime


import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


# --------------------------------------------------------------------------
# 1. Parsing
# --------------------------------------------------------------------------

def parse_powder_file(path: str):
    """Parse a multi-block OUTPUT_diffraction_powder*.dat file.

    Returns
    -------
    col_names : list[str]
        Column names from the '#Angle ...' header line, e.g.
        ['Angle', 'Total', 'In-In', 'In-O', 'O-O'].
    blocks : list[tuple[float, np.ndarray]]
        (time_fs, data) pairs, data has shape (n_angles, n_columns) with
        columns in the same order as col_names.
    """
    col_names = None
    blocks = []
    current_time = None
    current_rows = []

    with open(path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                content = line[1:].strip()
                try:
                    t = float(content)
                except ValueError:
                    # column-header line, e.g. "Angle  Total  In-In  In-O  O-O"
                    col_names = content.split()
                    continue
                # new block: flush the previous one first
                if current_time is not None and current_rows:
                    blocks.append((current_time, np.array(current_rows)))
                current_time = t
                current_rows = []
            else:
                current_rows.append([float(x) for x in line.split()])

    if current_time is not None and current_rows:
        blocks.append((current_time, np.array(current_rows)))

    if col_names is None:
        raise ValueError(f"Could not find the '#Angle ...' header line in {path}")

    return col_names, blocks


@dataclass
class PeakInfo:
    """Nominal positions of the Miller-index reflections, from the peaks file
    header (or supplied by hand)."""
    hkl: list[str]
    twotheta: np.ndarray
    q: np.ndarray = None


def parse_peaks_header(path: str) -> PeakInfo:
    """Read the 3-4 comment lines of an OUTPUT_diffraction_peaks*.dat file
    to get the nominal 2theta (and q) position of each named reflection."""
    with open(path, "r") as f:
        header_lines = []
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break

    hkl_line = next(l for l in header_lines if "Time" in l)
    tth_line = next(l for l in header_lines if "2theta" in l)
    hkl = re.findall(r"\(([^)]+)\)", hkl_line)
    twotheta = np.array([float(x) for x in re.findall(r"\(([-\d.eE+]+)\)", tth_line)])

    q = None
    q_line = next((l for l in header_lines if l.strip().lower().startswith("#q")), None)
    if q_line is not None:
        q = np.array([float(x) for x in re.findall(r"\(([-\d.eE+]+)\)", q_line)])

    assert len(hkl) == len(twotheta), "hkl labels and 2theta positions don't match in count"
    return PeakInfo(hkl=hkl, twotheta=twotheta, q=q)


def parse_peaks_data(path: str, peak_info: PeakInfo) -> pd.DataFrame:
    """Read the time-resolved I(hkl) table from OUTPUT_diffraction_peaks*.dat
    (this is XTANT's own normalized per-hkl intensity -- useful to compare
    against what this script's own fit extracts from the powder pattern)."""
    data = np.loadtxt(path, comments="#")
    cols = ["time_fs"] + list(peak_info.hkl)
    return pd.DataFrame(data, columns=cols)




def auto_peak_info(powder_path: str,
                   intensity_col="Total",
                   theta_min=10.0,
                   theta_max=None,
                   baseline_lam=3e4,
                   baseline_p=0.001,
                   prominence=0.05) -> PeakInfo:
    """
    Detect peak positions automatically from the initial block of a powder file.
    Excludes the forward-scattering divergence below theta_min.
    """
    col_names, blocks = parse_powder_file(powder_path)
    acol = col_names.index("Angle")
    icol = col_names.index(intensity_col)

    # take first block only
    t0, data0 = blocks[0]
    two_theta_full = data0[:, acol]
    intensity_full = data0[:, icol]

    # baseline correction, excluding forward scattering < theta_min
    two_theta, intensity, _ = subtract_baseline(
        two_theta_full, intensity_full,
        theta_min=theta_min, theta_max=theta_max,
        lam=baseline_lam, p=baseline_p)

    # detect peaks
    peak_indices, _ = find_peaks(intensity, prominence=prominence*np.max(intensity))
    centers = two_theta[peak_indices]

    # build PeakInfo (no hkl labels, just generic IDs)
    hkl_labels = [f"auto{i}" for i in range(len(centers))]
    return PeakInfo(hkl=hkl_labels, twotheta=np.array(centers))



# --------------------------------------------------------------------------
# 1b. Baseline (diffuse-background) removal
# --------------------------------------------------------------------------

def als_baseline(y, lam=3e4, p=0.001, niter=15):
    """Asymmetric Least Squares baseline (Eilers & Boelens, 2005). Follows
    the smooth diffuse/continuum background while ignoring the (positive)
    Bragg peaks sitting on top of it.

    lam : smoothness -- larger = stiffer/smoother baseline. For a 400-atom
        cluster's powder pattern, something in the 1e3-1e5 range over a
        ~40-50 deg window is typical; always check the plot.
    p : asymmetry -- smaller p pulls the baseline further below the peaks
        (0.001-0.01 is typical).

    IMPORTANT: run this on a 2theta range that excludes the huge forward-
    scattering divergence at very small angles (q->0). That feature is not
    diffraction information -- it swamps the dynamic range and wrecks the
    baseline fit everywhere else if left in. A `theta_min` of order 5-10 deg
    (well past the forward peak) is usually enough; check your own data.
    """
    y = np.asarray(y, dtype=float)
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2), dtype=float)
    D = lam * D.dot(D.transpose())
    w = np.ones(L)
    z = y.copy()
    for _ in range(niter):
        W = sparse.diags(w, 0)
        Z = W + D
        z = spsolve(Z.tocsc(), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def subtract_baseline(two_theta, intensity, theta_min=None, theta_max=None,
                       lam=3e4, p=0.001, niter=15):
    """Restrict to [theta_min, theta_max], subtract an ALS baseline, and
    return (two_theta_restricted, residual_intensity, baseline)."""
    mask = np.ones_like(two_theta, dtype=bool)
    if theta_min is not None:
        mask &= two_theta >= theta_min
    if theta_max is not None:
        mask &= two_theta <= theta_max
    x = two_theta[mask]
    y = intensity[mask]
    base = als_baseline(y, lam=lam, p=p, niter=niter)
    return x, y - base, base


# --------------------------------------------------------------------------
# 2. Peak-shape model
# --------------------------------------------------------------------------

_SQRT_PI_OVER_4LN2 = np.sqrt(np.pi / (4.0 * np.log(2.0)))


def pseudo_voigt(x, amp, center, fwhm, eta):
    """Pseudo-Voigt with shared FWHM for the Gaussian and Lorentzian parts.
    amp is the peak height (value at x=center)."""
    fwhm = max(fwhm, 1e-6)
    g = np.exp(-4.0 * np.log(2.0) * (x - center) ** 2 / fwhm ** 2)
    l = 1.0 / (1.0 + 4.0 * (x - center) ** 2 / fwhm ** 2)
    return amp * (eta * l + (1.0 - eta) * g)


def pseudo_voigt_area(amp, fwhm, eta):
    """Analytic integral of pseudo_voigt(x) dx from -inf to +inf."""
    area_g = amp * (1.0 - eta) * fwhm * _SQRT_PI_OVER_4LN2
    area_l = amp * eta * fwhm * (np.pi / 2.0)
    return area_g + area_l


# --------------------------------------------------------------------------
# 3. Grouping overlapping peaks + fitting
# --------------------------------------------------------------------------

def _merge_windows(centers: Sequence[float], half_width: float):
    """Merge nominal peak positions into groups whose fit windows
    [c-half_width, c+half_width] overlap. Returns list of (lo, hi, idx_list)."""
    order = np.argsort(centers)
    groups = []
    cur_lo, cur_hi, cur_idx = None, None, []
    for i in order:
        c = centers[i]
        lo, hi = c - half_width, c + half_width
        if cur_lo is None:
            cur_lo, cur_hi, cur_idx = lo, hi, [i]
        elif lo <= cur_hi:  # overlaps current group -> merge
            cur_hi = max(cur_hi, hi)
            cur_idx.append(i)
        else:
            groups.append((cur_lo, cur_hi, cur_idx))
            cur_lo, cur_hi, cur_idx = lo, hi, [i]
    if cur_idx:
        groups.append((cur_lo, cur_hi, cur_idx))
    return groups


@dataclass
class PeakFitResult:
    hkl_labels: list[str]      # one or more hkl sharing this position
    center: float
    center_nominal: float
    fwhm: float
    eta: float
    amp: float
    area: float
    window: tuple
    success: bool
    n_overlapping: int = 1


def _multi_pv_model(x, xref, n_peaks, n_bg, *params):
    bg_params = params[:n_bg]
    y = np.zeros_like(x)
    for p, c in zip(bg_params, range(n_bg)):
        y = y + p * (x - xref) ** c
    for k in range(n_peaks):
        amp, center, fwhm, eta = params[n_bg + 4 * k: n_bg + 4 * k + 4]
        y = y + pseudo_voigt(x, amp, center, fwhm, eta)
    return y


def _snap_to_local_min(x, y, edge, direction, search_deg):
    """Nudge a window edge outward/inward to the nearest local minimum of y
    within `search_deg` of `edge`, so the fit window boundary sits in a
    trough between reflections rather than slicing through a neighboring
    peak's rising flank. direction=+1 searches to the right (for a 'hi'
    edge), -1 to the left (for a 'lo' edge)."""
    if direction > 0:
        sel = (x >= edge - 0.2) & (x <= edge + search_deg)
    else:
        sel = (x >= edge - search_deg) & (x <= edge + 0.2)
    if sel.sum() < 3:
        return edge
    xs, ys = x[sel], y[sel]
    i_min = np.argmin(ys)
    return xs[i_min]


def fit_peak_group(two_theta, intensity, centers_nominal, hkl_groups,
                    half_width_deg=2.5, init_fwhm_deg=1.5,
                    center_tol_deg=1.0, bg_order=1, snap_to_minima=True,
                    snap_search_deg=2.0):
    """Fit one merged window that may contain several overlapping reflections.

    Parameters
    ----------
    two_theta, intensity : full-pattern arrays for this time step
    centers_nominal : nominal 2theta of each peak in this group
    hkl_groups : list of hkl-label-lists, one per peak in this group
    half_width_deg : half-window used to select data around the *group* span
    init_fwhm_deg : initial guess for FWHM of each peak
    center_tol_deg : how far the fitted center is allowed to drift from the
        nominal position (keeps the fit from sliding onto a neighboring peak)
    bg_order : degree of the local polynomial background (0=constant,
        1=linear, 2=quadratic). Small AIMD boxes (broad peaks, curved diffuse
        background) usually need 2 once the window gets wider than ~2 deg;
        a small box with genuinely sharp peaks can use 1.
    """
    lo = min(centers_nominal) - half_width_deg
    hi = max(centers_nominal) + half_width_deg
    if snap_to_minima:
        lo = _snap_to_local_min(two_theta, intensity, lo, -1, snap_search_deg)
        hi = _snap_to_local_min(two_theta, intensity, hi, +1, snap_search_deg)
    mask = (two_theta >= lo) & (two_theta <= hi)
    x = two_theta[mask]
    y = intensity[mask]
    if x.size < 8:
        return []

    xref = np.mean(centers_nominal)
    n_peaks = len(centers_nominal)
    n_bg = bg_order + 1

    # initial guesses: constant term from a low percentile, higher orders at 0
    p0 = [np.percentile(y, 10)] + [0.0] * bg_order
    lower = [-np.inf] * n_bg
    upper = [np.inf] * n_bg
    for c in centers_nominal:
        amp_init = max(np.interp(c, x, y) - p0[0], 1e-6)
        p0 += [amp_init, c, init_fwhm_deg, 0.5]
        lower += [0.0, c - center_tol_deg, 0.02, 0.0]
        upper += [np.inf, c + center_tol_deg, 2 * half_width_deg, 1.0]

    def model(xx, *params):
        return _multi_pv_model(xx, xref, n_peaks, n_bg, *params)

    try:
        popt, _ = curve_fit(model, x, y, p0=p0, bounds=(lower, upper), maxfev=30000)
        success = True
    except RuntimeError:
        popt, success = np.array(p0), False

    results = []
    for k in range(n_peaks):
        amp, center, fwhm, eta = popt[n_bg + 4 * k: n_bg + 4 * k + 4]
        area = pseudo_voigt_area(amp, fwhm, eta)
        results.append(PeakFitResult(
            hkl_labels=hkl_groups[k],
            center=center,
            center_nominal=centers_nominal[k],
            fwhm=fwhm,
            eta=eta,
            amp=amp,
            area=area,
            window=(lo, hi),
            success=success,
            n_overlapping=n_peaks,
        ))
    # stash the fit window/params on each result so plot_fit_group() can
    # redraw the curve + fit without recomputing anything
    for r in results:
        r._fit_x, r._fit_y, r._fit_popt, r._fit_xref, r._fit_nbg = x, y, popt, xref, n_bg
    return results


def plot_fit_group(results, ax=None):
    """Diagnostic plot of one fitted window: raw data, total fit, and the
    individual peak + background components. Given very broad peaks and
    curved backgrounds (typical for small AIMD boxes), ALWAYS eyeball a few
    of these before trusting the integrated areas."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    r0 = results[0]
    x, y, popt, xref, n_bg = r0._fit_x, r0._fit_y, r0._fit_popt, r0._fit_xref, r0._fit_nbg
    xx = np.linspace(x.min(), x.max(), 400)
    ax.plot(x, y, "o", ms=3, color="0.4", label="data")
    total = _multi_pv_model(xx, xref, len(results), n_bg, *popt)
    ax.plot(xx, total, "-", color="crimson", lw=2, label="total fit")
    bg = np.zeros_like(xx)
    for p, c in zip(popt[:n_bg], range(n_bg)):
        bg = bg + p * (xx - xref) ** c
    ax.plot(xx, bg, "--", color="steelblue", lw=1.5, label="background")
    for k, r in enumerate(results):
        amp, center, fwhm, eta = popt[n_bg + 4 * k: n_bg + 4 * k + 4]
        ax.plot(xx, bg + pseudo_voigt(xx, amp, center, fwhm, eta), ":",
                lw=1.5, label=f"{'+'.join(r.hkl_labels)}")
    ax.set_xlabel(r"2$\theta$ (deg)")
    ax.set_ylabel("Intensity (arb. units)")
    ax.legend(fontsize=8)
    return ax


def fit_all_peaks(two_theta, intensity, peak_info: PeakInfo,
                   half_width_deg=2.5, dedupe_tol_deg=1e-3, **fit_kwargs):
    """Fit every reflection in peak_info against one powder curve.
    Reflections that fall at (numerically) the same 2theta -- symmetry-
    equivalent hkls -- are merged into a single fitted peak and reported
    together (their intensities are indistinguishable in the powder pattern).
    """
    # group nominal positions that are essentially identical (equivalent hkls)
    order = np.argsort(peak_info.twotheta)
    unique_centers = []
    unique_hkl_groups = []
    for i in order:
        c = peak_info.twotheta[i]
        if unique_centers and abs(c - unique_centers[-1]) < dedupe_tol_deg:
            unique_hkl_groups[-1].append(peak_info.hkl[i])
        else:
            unique_centers.append(c)
            unique_hkl_groups.append([peak_info.hkl[i]])
    unique_centers = np.array(unique_centers)

    # merge overlapping fit windows (peaks closer than 2*half_width)
    groups = _merge_windows(unique_centers, half_width_deg)

    all_results = []
    for lo, hi, idx in groups:
        centers = [unique_centers[i] for i in idx]
        hkl_groups = [unique_hkl_groups[i] for i in idx]
        res = fit_peak_group(two_theta, intensity, centers, hkl_groups,
                              half_width_deg=half_width_deg, **fit_kwargs)
        all_results.extend(res)
    return all_results


# --------------------------------------------------------------------------
# 4. Driving the analysis over all time blocks of one run
# --------------------------------------------------------------------------

def analyze_run(powder_path: str, peak_info: PeakInfo, intensity_col="Total",
                 half_width_deg=2.5, theta_min=8.0, theta_max=None,
                 baseline_lam=3e4, baseline_p=0.001, **fit_kwargs) -> pd.DataFrame:
    """Fit every time block of one run's powder file. Returns a tidy
    DataFrame: one row per (time, reflection).

    theta_min/theta_max restrict the analysis range (passed to the ALS
    baseline step) -- theta_min should sit past the forward-scattering
    divergence at q->0; see subtract_baseline(). Set baseline_lam=None to
    skip baseline subtraction entirely (e.g. if you already background-
    corrected the data yourself).
    """

    col_names, blocks = parse_powder_file(powder_path)
    icol = col_names.index(intensity_col)
    acol = col_names.index("Angle")

    rows = []
    for step_idx, (t, data) in enumerate(blocks, start=1):
        start_time = datetime.now()
        print(f"[Step {step_idx}/{len(blocks)}] Starting analysis at simulation time {t:.2f} fs "
              f"(real time {start_time.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]})")

        two_theta_full = data[:, acol]
        intensity_full = data[:, icol]

        if baseline_lam is not None:
            two_theta, intensity, _ = subtract_baseline(
                two_theta_full, intensity_full, theta_min=theta_min,
                theta_max=theta_max, lam=baseline_lam, p=baseline_p)
        else:
            two_theta, intensity = two_theta_full, intensity_full

        results = fit_all_peaks(two_theta, intensity, peak_info,
                                 half_width_deg=half_width_deg, **fit_kwargs)

        end_time = datetime.now()
        print(f"               Finished time {t:.2f} fs, extracted {len(results)} peaks "
              f"(real time {end_time.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]})")

        for r in results:
            rows.append({
                "time_fs": t,
                "hkl": "+".join(r.hkl_labels),
                "n_equivalent_hkl": len(r.hkl_labels),
                "n_overlapping_peaks": r.n_overlapping,
                "center_deg": r.center,
                "center_nominal_deg": r.center_nominal,
                "fwhm_deg": r.fwhm,
                "eta": r.eta,
                "amplitude": r.amp,
                "integrated_area": r.area,
                "fit_ok": r.success,
            })

    print("✅ Analysis complete.")
    return pd.DataFrame(rows)



# --------------------------------------------------------------------------
# 5. Multi-run averaging
# --------------------------------------------------------------------------

def average_powder_curves(run_paths: Sequence[str]):
    """Average raw powder curves (Total + all element columns) point-by-point
    across several runs, block (time) by block. Assumes every run has the
    same time stamps and the same 2theta grid (true if generated with the
    same simulation settings).

    Returns (col_names, averaged_blocks, std_blocks) where averaged_blocks
    and std_blocks have the same shape/format as parse_powder_file's `blocks`
    (the std is across runs, useful as an error bar on the curve itself).
    """
    col_names_ref = None
    per_run_blocks = []
    for p in run_paths:
        col_names, blocks = parse_powder_file(p)
        if col_names_ref is None:
            col_names_ref = col_names
        elif col_names != col_names_ref:
            raise ValueError(f"{p} has different columns than the first run")
        per_run_blocks.append(blocks)

    n_blocks = len(per_run_blocks[0])
    for blocks in per_run_blocks:
        if len(blocks) != n_blocks:
            raise ValueError("Runs have different numbers of time blocks")

    avg_blocks, std_blocks = [], []
    for b in range(n_blocks):
        t_ref = per_run_blocks[0][b][0]
        stacked = np.stack([per_run_blocks[r][b][1] for r in range(len(run_paths))], axis=0)
        # keep the (identical) 2theta column, average the rest
        avg_data = stacked[0].copy()
        avg_data[:, 1:] = stacked[:, :, 1:].mean(axis=0)
        std_data = np.zeros_like(avg_data)
        std_data[:, 1:] = stacked[:, :, 1:].std(axis=0, ddof=1) if len(run_paths) > 1 else 0.0
        avg_blocks.append((t_ref, avg_data))
        std_blocks.append((t_ref, std_data))

    return col_names_ref, avg_blocks, std_blocks


def analyze_run_from_blocks(col_names, blocks, peak_info: PeakInfo,
                             intensity_col="Total", half_width_deg=2.5,
                             theta_min=8.0, theta_max=None,
                             baseline_lam=3e4, baseline_p=0.001,
                             **fit_kwargs) -> pd.DataFrame:
    """Same as analyze_run, but starting from already-parsed (e.g. averaged)
    blocks instead of re-reading a file."""
    icol = col_names.index(intensity_col)
    acol = col_names.index("Angle")
    rows = []
    for t, data in blocks:
        two_theta_full = data[:, acol]
        intensity_full = data[:, icol]
        if baseline_lam is not None:
            two_theta, intensity, _ = subtract_baseline(
                two_theta_full, intensity_full, theta_min=theta_min,
                theta_max=theta_max, lam=baseline_lam, p=baseline_p)
        else:
            two_theta, intensity = two_theta_full, intensity_full
        results = fit_all_peaks(two_theta, intensity, peak_info,
                                 half_width_deg=half_width_deg, **fit_kwargs)
        for r in results:
            rows.append({
                "time_fs": t, "hkl": "+".join(r.hkl_labels),
                "n_equivalent_hkl": len(r.hkl_labels),
                "center_deg": r.center, "center_nominal_deg": r.center_nominal,
                "fwhm_deg": r.fwhm, "eta": r.eta, "amplitude": r.amp,
                "integrated_area": r.area, "fit_ok": r.success,
            })
    return pd.DataFrame(rows)


def analyze_many_runs(run_paths: Sequence[str], peak_info: PeakInfo,
                       intensity_col="Total", half_width_deg=2.5,
                       **fit_kwargs) -> pd.DataFrame:
    """Fit each run individually, then average the fitted integrated areas
    (and FWHM) across runs and report the run-to-run std as an error bar.
    Use this in addition to fitting the pre-averaged curve (analyze_run_from_
    blocks on average_powder_curves' output) -- the pre-averaged-curve fit is
    the primary, lower-noise result; this one gives you the uncertainty.
    """
    per_run = []
    for i, p in enumerate(run_paths):
        df = analyze_run(p, peak_info, intensity_col=intensity_col,
                          half_width_deg=half_width_deg, **fit_kwargs)
        df["run"] = i
        per_run.append(df)
    all_runs = pd.concat(per_run, ignore_index=True)

    summary = (all_runs
               .groupby(["time_fs", "hkl"])
               .agg(integrated_area_mean=("integrated_area", "mean"),
                    integrated_area_std=("integrated_area", "std"),
                    fwhm_mean=("fwhm_deg", "mean"),
                    fwhm_std=("fwhm_deg", "std"),
                    center_mean=("center_deg", "mean"),
                    n_runs=("integrated_area", "count"))
               .reset_index())
    return summary, all_runs





def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    powder_files = glob.glob("OUTPUT_diffraction_powder*.dat")
    if not powder_files:
        print("No powder files found in the current directory.")
        input("Press Enter to exit...")
        return

    # Analyze each powder file
    for pf in powder_files:
        print(f"Analyzing {pf}...")
        # automatically detect peaks from the initial block
        peak_info = auto_peak_info(pf, theta_min=10.0, prominence=0.05)
        df = analyze_run(pf, peak_info, theta_min=10.0)
        out_name = os.path.splitext(pf)[0] + "_fit_automatic.csv"
        df.to_csv(out_name, index=False)
        print(f"Results saved to {out_name}")

    input("Press Enter to exit...")


if __name__ == "__main__":
    main()
