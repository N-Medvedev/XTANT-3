#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np
import re
import sys
import os
import csv
from pathlib import Path
from typing import Dict, Optional


# Module variables:
starline="*"*70



def main():


    # Define miller indices for diffraction peaks:
    hkl_list = [
        (2,0,0),
        (2,2,0),
        (3,1,1),
        (1,1,1)
    ]

    # Define the size of the supercell in the units of the unit-cell:
    n_cells=(3, 3, 3)


    filename = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "OUTPUT_atomic_coordinates.xyz"
    )

    print(f"Reading {filename}...")

    # Read the XYZ file:
    (
        cells,
        coords,
        species,
        kinetic_energies
    ) = read_extxyz(filename)

    print(
        f"Found {len(coords)} frames "
        f"with {coords[0].shape[0]} atoms."
    )


    # --------------------
    # Get MSDs:

    # Unwrapped trajectory
    traj = unwrap_trajectory(cells, coords)

    # Get the atom-by-atom displacements:
    dr0 = traj - traj[0]
    disp2_initial = np.sum(
        dr0 * dr0,
        axis=2
    )

    # Existing MSD from frame 0:
    msd = np.mean(
        disp2_initial,
        axis=1
    )


    # New time-origin averaged MSD
    msd_avg = compute_time_origin_averaged_msd(traj)

    # Analyse the MSD and get the average thermal MSD:
    (
        thermal_msd,
        thermal_msd_avg,
        mean_positions,
        disp2_thermal
    ) = compute_thermal_msd(traj)


    out = np.column_stack(
        (
            np.arange(len(msd)),
            msd,
        )
    )

    # Save the calcualted data:
    OUT_MSD_file = "OUT_MSD_analysis.dat"

    # Prepare the output file:
    header = (
    f"Average_thermal_MSD = "
    f"{thermal_msd_avg:.10e} Angstrom^2\n"
    f"RMS_thermal_displacement = "
    f"{np.sqrt(thermal_msd_avg):.10e} Angstrom\n"
    f"Frame "
    f"MSD_from_frame0_Angstrom2 "
    f"MSD_time_origin_averaged_Angstrom2 "
    f"MSD_thermal_Angstrom2"
    )


    # Save output file:
    np.savetxt(
        OUT_MSD_file,
        np.column_stack(
        (
            np.arange(len(msd)),
            msd,
            msd_avg,
            thermal_msd,
        )
        ),
        fmt=[
        "%8d",
        "%.10e",
        "%.10e",
        "%.10e",
        ],
        header=header
    )

    print(f"Saved MSD in the file: {OUT_MSD_file}")
    # And plot the MSDs:
    plot_msd(
        msd,
        msd_avg,
        thermal_msd
    )


    #------------------
    # Get Debye temperature:

    # Debye temperature evaluation:
    # Check if there is kinetic energy:
    if kinetic_energies is not None:
        print(
            "Found kinetic_energy property, evaluating Debye temperature."
        )
        # Next, get the atomic mass:
        masses_amu = get_atomic_masses(species)

        # Get the temperature:
        temperature = compute_temperature(kinetic_energies)

        # Get the Debye temperatures:
        theta_initial = compute_debye_temperature(
            temperature,
            disp2_initial,
            masses_amu
        )

        theta_thermal = compute_debye_temperature(
            temperature,
            disp2_thermal,
            masses_amu
        )

        # Get the Debye temperatures:
        OUT_Debye_file = ("OUT_Debye_temperature.dat")
        # Save them into file:
        save_debye_temperature(
            OUT_Debye_file,
            temperature,
            theta_initial,
            theta_thermal
        )
        # also plot them:
        plot_debye_temperature(
            temperature,
            theta_initial,
            theta_thermal,
            output_file="OUT_Debye_temperature.png")

    else:
        print(
            "No kinetic_energy property found, cannot obtain Debye temperature."
        )


    #---------------------
    # Get the DW diffraction intensities:

    intensities = compute_dw_peak_intensities(thermal_msd, cells, hkl_list, n_cells)

    save_dw_peak_intensities(
        "OUT_DW_thermal_peak_intensities.dat",
        intensities
    )

    plot_dw_peak_intensities(
        intensities, output_file="OUT_DW_thermal_peak_intensities.png", DW_title="DW thermal"
    )
    # Compare with the DW from MSD with respect to the zero-time:
    intensities = compute_dw_peak_intensities(msd, cells, hkl_list, n_cells)

    save_dw_peak_intensities(
        "OUT_DW_origin_peak_intensities.dat",
        intensities
    )

    plot_dw_peak_intensities(
        intensities, output_file="OUT_DW_origin_peak_intensities.png", DW_title="DW origin"
    )




def save_dw_peak_intensities(
        filename,
        intensities):
    """
    Save DW peak intensities.
    """

    hkls = list(
        intensities.keys()
    )

    nframes = len(
        intensities[hkls[0]]
    )

    data = [
        np.arange(nframes)
    ]

    labels = ["Frame"]

    for hkl in hkls:

        data.append(
            intensities[hkl]
        )

        labels.append(
            f"I_{hkl[0]}{hkl[1]}{hkl[2]}"
        )

    header = " ".join(labels)

    np.savetxt(
        filename,
        np.column_stack(data),
        fmt="%.10e",
        header=header
    )

    print(
        f"Saved diffraction intensities: "
        f"{filename}"
    )



def save_debye_temperature(
        filename,
        temperature,
        theta_initial,
        theta_thermal):

    avg_temperature = np.nanmean(temperature)
    avg_theta_initial = np.nanmean(theta_initial)
    avg_theta_thermal = np.nanmean(theta_thermal)

    header = (
        f"Average_temperature = "
        f"{avg_temperature:.10e} K\n"
        f"Average_Debye_temperature_initial = "
        f"{avg_theta_initial:.10e} K\n"
        f"Average_Debye_temperature_thermal = "
        f"{avg_theta_thermal:.10e} K\n"
        f"Frame "
        f"Temperature_K "
        f"ThetaD_initial_K "
        f"ThetaD_thermal_K"
    )
    print(header)

    np.savetxt(
        filename,
        np.column_stack(
            (
                np.arange(len(temperature)),
                temperature,
                theta_initial,
                theta_thermal
            )
        ),
        fmt=[
            "%8d",
            "%.10e",
            "%.10e",
            "%.10e"
        ],
        header=header
    )

    print(
        f"Saved Debye temperatures in file: "
        f"{filename}"
    )





def read_extxyz(filename):
    """
    Read extended XYZ trajectory.

    Returns
    -------
    cells : list of (3,3) arrays

    coords : list of (N,3) arrays

    species : ndarray (N,)
        Atomic symbols.

    kinetic_energies : list of arrays or None
        Per-atom kinetic energies [eV].
    """

    cells = []
    coords = []

    species = None

    kinetic_energies = []

    lattice_pattern = re.compile(
        r'Lattice="([^"]+)"'
    )

    with open(filename, "r") as f:

        while True:

            line = f.readline()

            if not line:
                break

            natoms = int(line.strip())

            comment = f.readline()

            # ----------------------------------
            # Cell
            # ----------------------------------

            m = lattice_pattern.search(comment)

            if not m:
                raise ValueError(
                    "Could not find Lattice field."
                )

            lattice = np.fromstring(
                m.group(1),
                sep=" "
            )

            if lattice.size != 9:
                raise ValueError(
                    f"Expected 9 lattice values, "
                    f"got {lattice.size}"
                )

            cell = lattice.reshape(3, 3)

            # ----------------------------------
            # Properties
            # ----------------------------------

            prop_map = parse_properties(comment)

            if "species" not in prop_map:
                raise ValueError(
                    "'species' property not found."
                )

            if "pos" not in prop_map:
                raise ValueError(
                    "'pos' property not found."
                )

            pos_start, pos_end = prop_map["pos"]

            has_ke = (
                "kinetic_energy" in prop_map
            )

            if has_ke:
                ke_col = prop_map[
                    "kinetic_energy"
                ][0]

                ke_frame = np.empty(
                    natoms,
                    dtype=float
                )

            xyz = np.empty(
                (natoms, 3),
                dtype=float
            )

            species_frame = []

            # ----------------------------------
            # Read atoms
            # ----------------------------------

            for i in range(natoms):

                cols = f.readline().split()

                species_frame.append(
                    cols[
                        prop_map["species"][0]
                    ]
                )

                xyz[i, :] = np.array(
                    cols[pos_start:pos_end],
                    dtype=float
                )

                if has_ke:
                    ke_frame[i] = float(
                        cols[ke_col]
                    )

            if species is None:
                species = np.array(
                    species_frame
                )

            cells.append(cell)
            coords.append(xyz)

            if has_ke:
                kinetic_energies.append(
                    ke_frame
                )

    if len(kinetic_energies) == 0:
        kinetic_energies = None
    else:
        kinetic_energies = np.array(
            kinetic_energies
        )

    return (
        cells,
        coords,
        species,
        kinetic_energies,
    )


def parse_properties(comment):

    m = re.search(r'Properties=([^ ]+)', comment)
    if not m:
        raise ValueError("Properties specification not found")

    fields = m.group(1).split(":")

    prop_map = {}

    col = 0
    i = 0

    while i < len(fields):

        name = fields[i]
        ptype = fields[i + 1]
        ncols = int(fields[i + 2])

        prop_map[name] = (col, col + ncols)

        col += ncols
        i += 3

    return prop_map



def cart_to_frac(r, cell):
    """
    Cartesian -> fractional.
    Works for arbitrary cell geometry.
    """
    return r @ np.linalg.inv(cell)


def frac_to_cart(s, cell):
    """
    Fractional -> Cartesian.
    """
    return s @ cell


def compute_msd(cells, coords):
    """
    Unwrap coordinates and calculate MSD
    relative to frame 0.
    """

    nframes = len(coords)
    natoms = coords[0].shape[0]

    # Initial frame
    frac_prev = cart_to_frac(coords[0], cells[0])

    # Continuous fractional coordinates
    frac_unwrapped = frac_prev.copy()
    frac_reference = frac_prev.copy()

    msd = np.zeros(nframes)

    # MSD at t = 0
    msd[0] = 0.0

    for iframe in range(1, nframes):

        frac_curr = cart_to_frac(coords[iframe], cells[iframe])

        ds = frac_curr - frac_prev

        # Minimum-image convention in fractional coordinates
        ds -= np.round(ds)

        frac_unwrapped += ds

        frac_prev = frac_curr

        r_unwrapped = frac_to_cart(
            frac_unwrapped,
            cells[iframe]
        )

        r_ref = frac_to_cart(
            frac_reference,
            cells[0]
        )

        dr = r_unwrapped - r_ref

        msd[iframe] = np.mean(
            np.sum(dr * dr, axis=1)
        )

    return msd


def compute_time_origin_averaged_msd(traj):
    """
    Compute time-origin averaged MSD.

    Parameters
    ----------
    traj : ndarray, shape (nframes, natoms, 3)
        Unwrapped Cartesian coordinates.

    Returns
    -------
    msd_avg : ndarray, shape (nframes,)
        Time-origin averaged MSD for each lag.
    """

    nframes, natoms, _ = traj.shape

    msd_avg = np.zeros(nframes)

    for lag in range(nframes):

        n_origins = nframes - lag

        msd_sum = 0.0

        for t0 in range(n_origins):

            dr = traj[t0 + lag] - traj[t0]

            msd_sum += np.mean(
                np.sum(dr * dr, axis=1)
            )

        msd_avg[lag] = msd_sum / n_origins

    return msd_avg


def unwrap_trajectory(cells, coords):
    """
    Construct unwrapped Cartesian coordinates.

    Returns
    -------
    traj : ndarray (nframes, natoms, 3)
    """

    nframes = len(coords)
    natoms = coords[0].shape[0]

    traj = np.empty((nframes, natoms, 3))

    frac_prev = cart_to_frac(coords[0], cells[0])
    frac_unwrapped = frac_prev.copy()

    traj[0] = coords[0]

    for iframe in range(1, nframes):

        frac_curr = cart_to_frac(coords[iframe], cells[iframe])

        ds = frac_curr - frac_prev

        # minimum-image convention
        ds -= np.round(ds)

        frac_unwrapped += ds

        traj[iframe] = frac_to_cart(
            frac_unwrapped,
            cells[iframe]
        )

        frac_prev = frac_curr

    return traj


def compute_mean_displacement_statistics(traj):
    """
    Compute time-averaged displacement statistics relative
    to the initial frame.

    Returns
    -------
    mean_disp : ndarray
        Mean displacement magnitude for each frame.

    msd_mean_disp : ndarray
        Mean squared displacement relative to the
        time-averaged mean displacement.

    avg_mean_disp : float
        Time-average of mean displacement.

    avg_msd_mean_disp : float
        Time-average of MSD relative to average displacement.
    """

    dr = traj - traj[0]

    # displacement magnitude for every atom and frame
    disp_mag = np.linalg.norm(dr, axis=2)

    # <|r(t)-r(0)|>
    mean_disp = np.mean(disp_mag, axis=1)

    avg_mean_disp = np.mean(mean_disp)

    # <(|r(t)-r(0)| - <|r|>_time)^2>
    msd_mean_disp = np.mean(
        (disp_mag - avg_mean_disp) ** 2,
        axis=1
    )

    avg_msd_mean_disp = np.mean(msd_mean_disp)

    return (
        mean_disp,
        msd_mean_disp,
        avg_mean_disp,
        avg_msd_mean_disp,
    )


def compute_thermal_msd(traj):
    """
    Compute thermal MSD:

        MSD_thermal(t)
        = (1/N) sum_i |r_i(t)-<r_i>|^2

    where <r_i> is the time-averaged position
    of atom i.

    Parameters
    ----------
    traj : ndarray (nframes, natoms, 3)
        Unwrapped Cartesian coordinates.

    Returns
    -------
    thermal_msd : ndarray (nframes,)
        Instantaneous thermal MSD.

    thermal_msd_avg : float
        Time-averaged thermal MSD.

    mean_positions : ndarray (natoms,3)
        Time-averaged atomic positions.
    """

    # <r_i>
    mean_positions = np.mean(traj, axis=0)

    dr = traj - mean_positions[None, :, :]

    disp2_thermal = np.sum(
        dr * dr,
        axis=2
    )

    thermal_msd = np.mean(
        np.sum(dr * dr, axis=2),
        axis=1
    )


    thermal_msd_avg = np.mean(thermal_msd)

    return (
        thermal_msd,
        thermal_msd_avg,
        mean_positions,
        disp2_thermal
    )



def compute_temperature(kinetic_energies):
    """
    Temperature from atomic kinetic energies.

    Parameters
    ----------
    kinetic_energies : ndarray
        Shape (nframes, natoms), eV

    Returns
    -------
    temperature : ndarray
        Temperature [K]
    """

    kb_ev = 8.617333262145e-5

    mean_ke = np.mean(
        kinetic_energies,
        axis=1
    )

    temperature = (
        2.0 * mean_ke
        / (3.0 * kb_ev)
    )

    return temperature


def compute_debye_temperature(
        temperature,
        disp2,
        masses_amu):
    """
    Compute classical Debye temperature from
    MSD and kinetic temperature.

    Formula:
        Theta_D = sqrt(
            9*hbar^2*T /
            (M_eff*kB*MSD)
        )

    Parameters
    ----------
    temperature : ndarray
        Temperature [K]

    disp2 : ndarray
        Squared atomic displacements [A^2]
        Shape (nframes, natoms)

    masses_amu : ndarray
        Atomic masses [amu]

    Returns
    -------
    theta_D : ndarray
        Debye temperature [K]
    """

    HBAR = 1.054571817e-34      # J*s
    KB = 1.380649e-23          # J/K
    AMU = 1.66053906660e-27    # kg

    masses = masses_amu * AMU

    # Effective atomic mass
    M_eff = np.mean(masses)

    theta_D = np.full(
        len(temperature),
        np.nan,
        dtype=float
    )

    for iframe in range(len(temperature)):

        T = temperature[iframe]

        if T <= 0.0:
            continue

        # MSD [A^2]
        msd = np.mean(
            disp2[iframe]
        )

        # A^2 -> m^2
        msd *= 1.0e-20

        if msd <= 0.0:
            continue

        theta_D[iframe] = np.sqrt(
            9.0
            * HBAR**2
            * T
            /
            (
                M_eff
                * KB
                * msd
            )
        )

    return theta_D



def compute_dw_peak_intensities(
        thermal_msd,
        cells,
        hkl_list,
        n_cells=(1, 1, 1),
        normalize=False):
    """
    Debye-Waller attenuation of diffraction peaks for a supercell
    built from n_cells unit cells along each axis.

    Parameters
    ----------
    thermal_msd : float or ndarray, shape (n_frames,)
        Total (3D) thermal MSD <u^2> [A^2], per frame or a single value.
    cells : list
        Simulation (super)cells, one per frame.
    hkl_list : list of (h, k, l)
        Miller indices with respect to the UNIT cell.
    n_cells : tuple of 3 ints
        Number of unit cells along each supercell axis, e.g. (2, 3, 4).
        Set once; the same numbers are used for all frames.
    normalize : bool
        If True, divide by the frame-0 value.

    Returns
    -------
    intensities : dict
        {(h, k, l): ndarray (n_frames,)}
    """
    n = np.asarray(n_cells, dtype=float)
    if n.shape != (3,) or np.any(n < 1) or not np.allclose(n, np.rint(n)):
        raise ValueError("n_cells must be three positive integers, e.g. (2, 3, 4)")

    n_frames = len(cells)
    msd = np.atleast_1d(np.asarray(thermal_msd, dtype=float))
    if msd.size == 1:
        msd = np.full(n_frames, msd[0])
    if msd.size != n_frames:
        raise ValueError("thermal_msd and cells must have the same number of frames")

    hkl = np.asarray(hkl_list, dtype=float)            # (n_hkl, 3)

    G2 = np.empty((n_frames, len(hkl)))
    for i, cell in enumerate(cells):
        B_super = np.asarray(reciprocal_lattice(cell), dtype=float)  # rows b1, b2, b3
        B_unit = B_super * n[:, None]                  # b_i(unit) = n_i * b_i(super)
        Gvec = hkl @ B_unit                            # (n_hkl, 3)
        G2[i] = np.einsum("ij,ij->i", Gvec, Gvec)

    I = np.exp(-G2 * msd[:, None] / 3.0)

    if normalize:
        I = I / I[0]

    return {tuple(int(x) for x in h): I[:, j] for j, h in enumerate(hkl)}


def reciprocal_lattice(cell):
    """
    Reciprocal lattice vectors [1/A].

    Parameters
    ----------
    cell : ndarray (3,3)

    Returns
    -------
    b1, b2, b3 : ndarray (3,)
    """

    a1 = cell[0]
    a2 = cell[1]
    a3 = cell[2]

    V = np.dot(
        a1,
        np.cross(a2, a3)
    )

    b1 = (
        2.0 * np.pi
        * np.cross(a2, a3)
        / V
    )

    b2 = (
        2.0 * np.pi
        * np.cross(a3, a1)
        / V
    )

    b3 = (
        2.0 * np.pi
        * np.cross(a1, a2)
        / V
    )

    return b1, b2, b3





##########################
# Plotting subroutines:


def plot_dw_peak_intensities(
        intensities,
        output_file="OUT_DW_peak_intensities.png",
        DW_title="Debye-Waller intensities"):
    """
    Plot Debye-Waller peak intensities.
    """

    hkls = list(
        intensities.keys()
    )

    nframes = len(
        intensities[hkls[0]]
    )

    frames = np.arange(
        nframes
    )

    fig, ax = plt.subplots(
        figsize=(5,4)
    )

    styles = [
        "-",
        "--",
        ":",
        "-."
    ]

    for i, hkl in enumerate(hkls):

        ax.plot(
            frames,
            intensities[hkl],
            lw=2.0,
            ls=styles[
                i % len(styles)
            ],
            label=f"({hkl[0]}{hkl[1]}{hkl[2]})"
        )

    ax.set_xlabel(
        "Frame",
        fontsize=14
    )

    ax.set_ylabel(
        r"$I/I_0$",
        fontsize=14
    )

    ax.set_ylim(
        0.0,
        1.05
    )

    ax.grid(
        False
    )

    ax.legend(
        frameon=False
    )

    ax.set_title(DW_title)

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(
        f"Saved figure: {output_file}"
    )




def plot_debye_temperature(
        temperature,
        theta_initial,
        theta_thermal,
        output_file="OUT_Debye_temperature.png"):
    """
    Plot Debye temperatures.
    """

    frames = np.arange(
        len(theta_initial)
    )

    avg_theta_initial = np.nanmean(
        theta_initial
    )

    avg_theta_thermal = np.nanmean(
        theta_thermal
    )

    avg_theta = max(
        avg_theta_initial,
        avg_theta_thermal
    )

    fig, ax = plt.subplots(
        figsize=(5, 4)
    )

    ax.plot(
        frames,
        theta_initial,
        lw=2.0,
        color="tab:blue",
        label=r"$\Theta_D$ from initial MSD"
    )

    ax.plot(
        frames,
        theta_thermal,
        lw=2.0,
        ls="--",
        color="tab:red",
        label=r"$\Theta_D$ from thermal MSD"
    )

    ax.set_xlabel(
        "Frame",
        fontsize=14
    )

    ax.set_ylabel(
        r"Debye temperature (K)",
        fontsize=14
    )

    # Set y-axis range
    ax.set_ylim(
        0.0,
        2.0 * avg_theta
    )

    ax.tick_params(
        axis="both",
        labelsize=12
    )

    ax.grid(False)

    ax.legend(
        fontsize=12,
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(
        f"Average Theta_D (initial): "
        f"{avg_theta_initial:.2f} K"
    )

    print(
        f"Average Theta_D (thermal): "
        f"{avg_theta_thermal:.2f} K"
    )

    print(
        f"Saved figure: {output_file}"
    )


def plot_msd(msd, msd_avg, thermal_msd,
             output_file="OUT_MSDs.png"):
    """
    Plot MSD curves and save to PNG.

    Parameters
    ----------
    msd : ndarray
        MSD relative to frame 0.

    msd_avg : ndarray
        Time-origin averaged MSD.

    thermal_msd : ndarray
        Thermal MSD with respect to the time-average position.

    output_file : str
        Output PNG filename.
    """

    frames = np.arange(len(msd))

    fig, ax = plt.subplots(figsize=(5, 4))

    ax.plot(
        frames,
        msd,
        lw=2.0,
        label="MSD from frame 0",
    )

    ax.plot(
        frames,
        msd_avg,
        lw=2.0,
        ls="--",
        label="Time-origin averaged MSD",
    )

    ax.plot(
        frames,
        thermal_msd,
        lw=2.0,
        ls=":",
        label="Thermal MSD",
    )

    ax.set_xlabel("Frame")
    ax.set_ylabel(r"MSD ($\AA^2$)")
    ax.grid(False)
    ax.legend()
    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)






##########################
# Dealing for the Periodic Table:

def get_atomic_masses(species):
    """
    Convert element symbols into masses [amu].

    Parameters
    ----------
    species : ndarray (natoms,)

    Returns
    -------
    masses : ndarray (natoms,)
        Atomic masses in amu.
    """

    atomic_data = load_atomic_data()

    masses = np.empty(len(species))

    for i, symbol in enumerate(species):

        elem = get_element(
            atomic_data,
            name=symbol
        )

        if elem is None:
            raise ValueError(
                f"Element '{symbol}' not found "
                f"in INPUT_atomic_data.txt"
            )

        masses[i] = elem["mass"]

    return masses




def load_atomic_data() -> Dict[int, dict]:
    """
    Load atomic data from INPUT/Atomic_data/INPUT_atomic_data.dat
    into a dictionary keyed by atomic number.

    Works reliably across Windows, Linux, and macOS by resolving
    the path relative to the project root (script or executable).

    Returns
    -------
    dict[int, dict]
        Dictionary of atomic properties keyed by atomic number.
    """
    # Resolve project root
    if getattr(sys, 'frozen', False):
        project_root = Path(sys.executable).parent
    else:
        project_root = Path(__file__).resolve().parent.parent

    # Build full path to atomic data file
    #full_path = project_root / "INPUT_atomic_data.dat"
    full_path = Path("INPUT_atomic_data.dat")

    if not full_path.is_file():
        print(starline)
        print(f"Error: Atomic data file not found: {full_path}", file=sys.stderr)
        print(starline)
        sys.exit(1)   # Exit with non‑zero status code

    atomic_data = {}
    with full_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            num = int(row["Atomic number"])
            atomic_data[num] = {
                "number": num,
                "name": row["Name"].strip(),
                "abbreviation": row["Abbreviation"].strip(),
                "mass": float(row["Mass"]),
                "valence_electrons": int(row["Number of valence electrons"]),
                "valent_radius": float(row["Valent radius"]),
                "electronegativity": float(row["Electronegativity"]),
            }
    return atomic_data


#def get_element(data: dict[int, dict], atomic_number: int = None, name: str = None) -> dict | None:
def get_element(data: Dict[int, Dict], atomic_number: int = None, name: str = None) -> Optional[Dict]:

    """
    Retrieve element data by atomic number or name/abbreviation.
    """
    if atomic_number is not None:
        return data.get(atomic_number)

    if name is not None:
        name = name.lower()
        for props in data.values():
            if props["name"].lower() == name or props["abbreviation"].lower() == name:
                return props

    return None



##########################
if __name__ == "__main__":
    main()
