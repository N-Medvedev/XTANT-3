"""
phonon_spectrum_vacf.py

Computes the phonon (vibrational density of states) spectrum from atomic
velocities via the velocity autocorrelation function (VACF), evaluated in
sliding time windows across an MD trajectory. A damping envelope is applied
to the VACF before Fourier transforming, which both suppresses truncation
noise and gives the spectral lines a finite width (the same role a damping
constant plays in a driven-damped-oscillator picture of phonon dephasing).

Inputs
------
OUTPUT_temperatures.dat
    Header: 2 comment lines. Column 0 = time [fs].
    Only the time grid is used from this file (to get dt and total range).

OUTPUT_atomic_velocities.xyz
    Standard extended-XYZ trajectory of velocities [A/fs]:
        <natoms>
        Lattice="..." Properties=species:S:1:vel:R:3
        <species> vx vy vz
        ... (natoms lines)
    One block per time step, in the same order as the time grid above.

Output
------
A single PNG with one VACF-derived spectrum per analysis window, each
window starting every WINDOW_STEP_FS fs and spanning WINDOW_LENGTH_FS fs.

Usage
-----
Just run the script (double-click or `python phonon_spectrum_vacf.py`).
All parameters are in the CONFIG block below.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ============================== CONFIG ==============================
SCRIPT_DIR = Path(__file__).resolve().parent

TEMP_FILE = SCRIPT_DIR / "OUTPUT_temperatures.dat"
VEL_FILE = SCRIPT_DIR / "OUTPUT_atomic_velocities.xyz"
OUTPUT_PNG = SCRIPT_DIR / "phonon_spectra.png"

# Windowing: a new VACF window starts every WINDOW_STEP_FS fs; each window
# is WINDOW_LENGTH_FS fs long. Set them equal for non-overlapping windows,
# or WINDOW_LENGTH_FS > WINDOW_STEP_FS for overlapping ones.
WINDOW_STEP_FS = 200.0
WINDOW_LENGTH_FS = 400.0

# Damping applied to the VACF before the Fourier transform. This controls
# the spectral linewidth: shorter DAMPING_TIME_FS -> broader, smoother
# peaks; longer -> sharper peaks but more truncation ringing.
# 'exponential': w(t) = exp(-t/DAMPING_TIME_FS)          (Lorentzian-like broadening)
# 'gaussian':    w(t) = exp(-0.5*(t/DAMPING_TIME_FS)**2)  (Gaussian broadening)
#DAMPING_TYPE = "exponential"
DAMPING_TYPE = "gaussian"
DAMPING_TIME_FS = 100.0

# Plot / frequency axis
FREQ_UNIT = "THz"          # 'THz' or 'cm-1'
MAX_FREQ = 30.0            # upper limit of the frequency axis, in FREQ_UNIT
NORMALIZE_CURVES = False    # scale each window's spectrum to its own max
# ======================================================================


def read_time_grid(path):
    """Read the time column [fs] from OUTPUT_temperatures.dat (2 header lines)."""
    data = np.loadtxt(path, skiprows=2)
    times = data[:, 0]
    return times


def read_velocity_trajectory(path):
    """
    Parse an extended-XYZ velocity trajectory.

    Returns
    -------
    species : list[str]      length natoms, species labels of the first block
    vel     : ndarray        shape (nsteps, natoms, 3), velocities [A/fs]
    """
    species = None
    frames = []

    with open(path, "r") as f:
        lines = f.readlines()

    i = 0
    nlines = len(lines)
    while i < nlines:
        line = lines[i].strip()
        if line == "":
            i += 1
            continue
        natoms = int(line)
        i += 1  # skip natoms line
        i += 1  # skip the "Lattice=... Properties=..." comment line

        block_species = []
        block_vel = np.empty((natoms, 3), dtype=float)
        for a in range(natoms):
            parts = lines[i].split()
            block_species.append(parts[0])
            block_vel[a, :] = [float(parts[1]), float(parts[2]), float(parts[3])]
            i += 1

        if species is None:
            species = block_species
        frames.append(block_vel)

    vel = np.stack(frames, axis=0)  # (nsteps, natoms, 3)
    return species, vel


def vacf_normalized(vel_window):
    """
    Velocity autocorrelation function for a single time window, averaged
    over all atoms and Cartesian components, normalized so C(0) = 1.

    vel_window : ndarray, shape (M, natoms, 3)

    Returns
    -------
    C : ndarray, shape (M,)   C[k] = <v(t). v(t+k)> / <v(t).v(t)>
    """
    M, natoms, _ = vel_window.shape
    # Flatten atom and Cartesian axes into independent "channels" to average over.
    V = vel_window.reshape(M, natoms * 3)  # (M, channels)

    # FFT-based linear autocorrelation (avoids O(M^2) direct loop).
    nfft = 1
    while nfft < 2 * M:
        nfft *= 2
    Vf = np.fft.rfft(V, n=nfft, axis=0)
    power = Vf * np.conj(Vf)
    ac_full = np.fft.irfft(power, n=nfft, axis=0).real  # (nfft, channels)
    ac = ac_full[:M, :]  # lags 0..M-1

    # Unbiased normalization: divide lag k by the number of overlapping pairs (M-k),
    # then average over channels.
    counts = (M - np.arange(M)).astype(float)
    ac = ac / counts[:, None]
    C = ac.mean(axis=1)

    if C[0] != 0:
        C = C / C[0]
    return C


def apply_damping(C, dt_fs, tau_fs, kind="exponential"):
    """Multiply the VACF by a damping envelope w(t)."""
    t = np.arange(len(C)) * dt_fs
    if kind == "exponential":
        w = np.exp(-t / tau_fs)
    elif kind == "gaussian":
        w = np.exp(-0.5 * (t / tau_fs) ** 2)
    else:
        raise ValueError(f"Unknown DAMPING_TYPE '{kind}'")
    return C * w


def spectrum_from_vacf(C_damped, dt_fs):
    """
    Cosine-transform the (damped) VACF to get the phonon spectrum.
    C is assumed even in time (C(-t) = C(t)), so we mirror it and take a
    real FFT; the result is equivalent to a direct cosine transform.

    Returns
    -------
    freq  : ndarray, frequencies in 1/fs
    power : ndarray, spectral intensity (non-negative frequencies only)
    """
    M = len(C_damped)
    # Build the symmetric (even) extension: C[M-1], ..., C[1], C[0], C[1], ..., C[M-1]
    full = np.concatenate([C_damped[::-1], C_damped[1:]])
    spec = np.fft.rfft(full)
    power = np.abs(spec)
    freq = np.fft.rfftfreq(len(full), d=dt_fs)  # 1/fs
    return freq, power


def convert_freq(freq_per_fs, unit):
    """Convert frequency from 1/fs to THz or cm^-1."""
    freq_THz = freq_per_fs * 1000.0  # 1/fs = 1e15 Hz = 1000 THz
    if unit == "THz":
        return freq_THz
    elif unit == "cm-1":
        return freq_THz * 33.35641  # THz -> cm^-1
    else:
        raise ValueError(f"Unknown FREQ_UNIT '{unit}'")


def main():
    times = read_time_grid(TEMP_FILE)
    species, vel = read_velocity_trajectory(VEL_FILE)

    nsteps_vel = vel.shape[0]
    if nsteps_vel != len(times):
        print(
            f"Warning: {len(times)} time points in {TEMP_FILE.name} but "
            f"{nsteps_vel} frames in {VEL_FILE.name}. Using the shorter of the two."
        )
    nsteps = min(len(times), nsteps_vel)
    times = times[:nsteps]
    vel = vel[:nsteps]

    dt_fs = float(np.median(np.diff(times)))
    natoms = vel.shape[1]
    print(f"Loaded {nsteps} steps, dt = {dt_fs:.4f} fs, {natoms} atoms.")

    t0_grid = np.arange(times[0], times[-1] - WINDOW_LENGTH_FS + 1e-9, WINDOW_STEP_FS)
    if len(t0_grid) == 0:
        t0_grid = np.array([times[0]])

    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("viridis")

    any_curve = False
    for i, t0 in enumerate(t0_grid):
        t1 = t0 + WINDOW_LENGTH_FS
        idx = np.where((times >= t0) & (times < t1))[0]
        if len(idx) < 4:
            continue  # window too short to be meaningful

        vel_window = vel[idx]
        C = vacf_normalized(vel_window)
        C_damped = apply_damping(C, dt_fs, DAMPING_TIME_FS, DAMPING_TYPE)
        freq_per_fs, power = spectrum_from_vacf(C_damped, dt_fs)
        freq = convert_freq(freq_per_fs, FREQ_UNIT)

        if NORMALIZE_CURVES and power.max() > 0:
            power = power / power.max()

        mask = freq <= MAX_FREQ
        color = cmap(i / max(1, len(t0_grid) - 1))
        ax.plot(freq[mask], power[mask], color=color,
                label=f"{t0:.0f}\u2013{t1:.0f} fs")
        any_curve = True

    if not any_curve:
        raise RuntimeError("No windows produced a spectrum - check WINDOW_LENGTH_FS vs. trajectory length.")

    xlabel = "Frequency (THz)" if FREQ_UNIT == "THz" else "Wavenumber (cm$^{-1}$)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("VDOS (normalized)" if NORMALIZE_CURVES else "VDOS (arb. units)")
    ax.set_xlim(0, MAX_FREQ)
    ax.set_title("Phonon spectrum from velocity autocorrelation (damped)")
    ax.legend(title="Window", fontsize=8, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=200)
    print(f"Saved: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()