#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import re
import sys
import os
import csv
from pathlib import Path
from typing import Dict, Optional


#----------------------
# Constants:
# Module variables:
starline="*"*70
m_two_third = 2.0 / 3.0
m_five_sixteenth = 5.0 / 16.0
m_Thom_factor = 20.6074     # Thomson factor (Eq. 2.5 [1])
g_me = 9.1093821545e-31     # [kg] Electron mass
g_alpha = 0.00729735256     # fine structure constant
g_cvel = 299792458.0        # [m/s] speed of light

#----------------------
# Output filenames:
OUT_average_positions   = "OUT_average_positions.xyz"
OUT_MSD_file            = "OUT_MSD_analysis.dat"
OUT_Debye_file          = "OUT_Debye_temperature.dat"
OUT_DW_thermal_peak_intensities         = "OUT_DW_thermal_peak_intensities.dat"
OUT_DW_origin_peak_intensities          = "OUT_DW_origin_peak_intensities.dat"
OUT_full_peak_intensities               = "OUT_full_peak_intensities.dat"
OUT_full_peak_intensities_normalized    = "OUT_full_peak_intensities_normalized.dat"


#----------------------
# Parameters of the analysis:
# Define miller indices for diffraction peaks:
hkl_list = [
        (2,0,0),
        (1,0,1),
        (2,1,1),
        (2,2,0),
        (1,1,2),
        (0,0,4)
        ]

# Define the size of the supercell in the units of the unit-cell:
n_cells=(4, 4, 7)


# File with reference coordinates (ideal structure):
#IN_reference_atomic_positions = "OUT_average_positions.xyz"
IN_reference_atomic_positions = "None"  # must be string
# Option to use the displacements from the reference points instead of the absolute coordinates:
use_displacements_set=False         # option to use : R(t) =  r_ref + (r(t) - r(0))

#----------------------



def main():

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

    # Get time grid, if possible:
    time_file = "OUTPUT_temperatures.dat"
    time = read_time_grid(time_file)   # None if the file is absent


    if time is not None:
        print(f"Time grid read from file {time_file}")
    else:
        print(f"Time grid not found in file {time_file}, using frames instead")

    # --------------------
    # Get MSDs:

    # Unwrapped trajectory
    traj = unwrap_trajectory(cells, coords)

    # Check if reference atomic positions are provided:
    ref_positions = None
    if os.path.isfile(IN_reference_atomic_positions):
        ref_positions, ref_cell = load_reference_positions(IN_reference_atomic_positions, traj[0], cells[0], species=species)
        ref_cell = np.asarray(ref_cell, dtype=float).reshape(-1, 3, 3)[0]   # reshape to ensure required shape

        use_displacements = use_displacements_set
        print(f"Using reference positions from {IN_reference_atomic_positions}")
    else:
        ref_positions = traj[0]
        ref_cell = cells[0]
        use_displacements = False
        print("Reference file not found; using frame 0 as reference.")

    # Atom-by-atom displacements from the reference
    dr0 = traj - ref_positions[None, :, :]
    disp2_initial = np.sum(dr0 * dr0, axis=2)

    # MSD from the reference
    msd = np.mean(disp2_initial, axis=1)


    # New time-origin averaged MSD
    msd_avg = compute_time_origin_averaged_msd(traj)

    # Analyse the MSD and get the average thermal MSD:
    ( thermal_msd, thermal_msd_avg, mean_positions, disp2_thermal ) = compute_thermal_msd(traj)

    # Save the average atomic positions:
    write_extxyz(
        OUT_average_positions,
        species,
        mean_positions,
        cell=np.mean(cells, axis=0),                   # average cell
        extra={"thermal_msd": disp2_thermal.mean(axis=0)},   # per-atom <|dr|^2> [A^2]
        comment_info={"thermal_msd_avg": f"{thermal_msd_avg:.8e}"},
        wrap=True,
    )



    out = np.column_stack( ( np.arange(len(msd)), msd, ) )

    if time is not None:
        x_label = "Time (fs)"
    else:
        x_label = "Frame"

    # Prepare the output file:
    header = (f"Average_thermal_MSD = "
    f"{thermal_msd_avg:.10e} Angstrom^2\n"
    f"RMS_thermal_displacement = "
    f"{np.sqrt(thermal_msd_avg):.10e} Angstrom\n"
    f"{x_label} "
    f"MSD_from_frame0_Angstrom2 "
    f"MSD_time_origin_averaged_Angstrom2 "
    f"MSD_thermal_Angstrom2"
    )

    # Save output file:
    if time is not None and len(time) == len(msd):   # use defined time:
        first_col = np.asarray(time, dtype=float)
        first_fmt = "%.6e"
    else:                                            # use frame numbers instead of time:
        first_col = np.arange(len(msd))
        first_fmt = "%8d"

    np.savetxt(
        OUT_MSD_file,
        np.column_stack((first_col, msd, msd_avg, thermal_msd)),
        fmt=[first_fmt, "%.10e", "%.10e", "%.10e"],
        header=header
    )

    print(f"Saved MSD in the file: {OUT_MSD_file}")
    # And plot the MSDs:
    plot_msd(msd, msd_avg, thermal_msd, x_label=x_label)


    #------------------
    # Get Debye temperature:

    # Debye temperature evaluation:
    # Check if there is kinetic energy:
    if kinetic_energies is not None:
        print("Found kinetic_energy property, evaluating Debye temperature.")
        # Next, get the atomic mass:
        masses_amu = get_atomic_masses(species)

        # Get the temperature:
        temperature = compute_temperature(kinetic_energies)

        # Get the Debye temperatures:
        # For displacements with respect to initial:
        theta_initial = compute_debye_temperature(temperature, disp2_initial, masses_amu)
        # For displacements with respect to thermal:
        theta_thermal = compute_debye_temperature(temperature, disp2_thermal, masses_amu)

        # Get the Debye temperatures:
        # Save them into file:
        save_debye_temperature(OUT_Debye_file, temperature, theta_initial, theta_thermal, time=time)
        # also plot them:
        plot_debye_temperature(temperature, theta_initial, theta_thermal, output_file="OUT_Debye_temperature.png", x_label=x_label)

    else:
        print(
            "No kinetic_energy property found, cannot obtain Debye temperature."
        )


    #---------------------
    # Get the DW diffraction intensities:
    # DW from MSD with respect to reference positions (or the frame from zero-time):
    intensities = compute_dw_peak_intensities(msd, cells, hkl_list, n_cells)
    save_dw_peak_intensities(OUT_DW_origin_peak_intensities, intensities)
    plot_dw_peak_intensities(intensities, time=time, output_file="OUT_DW_origin_peak_intensities.png", DW_title="DW origin")

    # DW from thermal MSD:
    intensities = compute_dw_peak_intensities(thermal_msd, cells, hkl_list, n_cells)
    save_dw_peak_intensities(OUT_DW_thermal_peak_intensities, intensities)
    plot_dw_peak_intensities(intensities, time=time, output_file="OUT_DW_thermal_peak_intensities.png", DW_title="DW thermal")


    # ----------------------------------------
    # Full diffraction calculation
    # ----------------------------------------
    intensities_full = compute_full_peak_intensities(
        cells=cells, coords=coords, species=species, hkl_list=hkl_list, n_cells=n_cells, normalize=False, )

    # Save raw intensities
    save_dw_peak_intensities(OUT_full_peak_intensities, intensities_full)

    # Plots full peaks:
    plot_dw_peak_intensities(intensities_full, time=time, output_file="OUT_full_peak_intensities_raw.png",
                             DW_title="Diffraction (raw peaks)")

    # Also calculate intensities normalized to the the given reference points (or first frame):
    #ref_positions = load_reference_positions("reference_positions.xyz", traj[0], cells[0], species=species)
    if ref_positions is not None:
        intensities_full_norm = compute_full_peak_intensities(
            cells, coords, species, hkl_list, n_cells=n_cells,
            normalize=True,
            reference_positions=ref_positions,
            reference_cell = ref_cell,
            normalize_to_reference=True,            # option 1
            use_displacements=use_displacements,    # option 2: set True to use r_ref + (r(t) - r(0))
            )
    else:   # use first frame as reference for normalization:
        intensities_full_norm = compute_full_peak_intensities(
            cells=cells, coords=coords, species=species, hkl_list=hkl_list, n_cells=n_cells, normalize=True, )

    # No need to recalculate, use precalculated full intensities:
    # for hkl in intensities_full:
    #     I0 = intensities_full[hkl][0]
    #     if I0 > 0.0:
    #         intensities_full_norm[hkl] = intensities_full[hkl] / I0[hkl]
    #     else:
    #         intensities_full_norm[hkl][:] = np.nan

    # Save normalized intensities
    save_dw_peak_intensities(OUT_full_peak_intensities_normalized, intensities_full_norm)

    # Plot normalized intensities
    plot_dw_peak_intensities(intensities_full_norm, time=time, output_file="OUT_full_peak_intensities_normalized.png",
                             DW_title="Diffraction (normalized)")




##########################
# Reading files:

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




def read_time_grid(filename="OUTPUT_temperatures.dat", skip_lines=2):
    """
    Read the time grid [fs] from the first column of the file.

    Returns
    -------
    time : ndarray or None
        Time in fs, or None if the file is absent or unreadable.
    """
    if not os.path.isfile(filename):
        return None
    try:
        time = np.loadtxt(filename, skiprows=skip_lines, usecols=0, ndmin=1)
    except (ValueError, OSError) as err:
        print(f"Warning: could not read time from {filename} ({err}); using frames.")
        return None
    return time



def get_x_axis(nframes, time=None):
    """
    Return (x, xlabel). Uses time [fs] if its length matches the number
    of frames, otherwise frame indices.
    """
    if time is not None:
        if len(time) == nframes:
            return time, "Time (fs)"
        print(f"Warning: time grid has {len(time)} points but there are "
              f"{nframes} frames; using frames instead.")
    return np.arange(nframes), "Frame"




def load_reference_positions(filename, traj0, cell0, species=None, tol_frac=0.5):
    """
    Read reference atomic positions and align them with the unwrapped trajectory.

    Parameters
    ----------
    filename : str
        Extended XYZ file with one (or more; the first is used) frame.
    traj0 : ndarray (N, 3)
        Unwrapped positions at frame 0 [A].
    cell0 : ndarray (3, 3)
        Cell at frame 0, lattice vectors as rows.
    species : ndarray (N,) or None
        Species of the trajectory, to verify atom order.

    Returns
    -------
    ref : ndarray (N, 3)
        Reference positions, each shifted by lattice vectors to the image
        closest to traj0.
    """
    ref_cells, ref_coords, ref_species, _ = read_extxyz(filename)
    ref = np.asarray(ref_coords[0], dtype=float)

    if ref.shape != traj0.shape:
        raise ValueError(
            f"Reference has {ref.shape[0]} atoms, trajectory has {traj0.shape[0]}."
        )
    if species is not None and not np.array_equal(ref_species, species):
        raise ValueError("Reference species/order differ from the trajectory.")

    # Minimum-image alignment in fractional coordinates of the frame-0 cell
    d_frac = (ref - traj0) @ np.linalg.inv(cell0)
    if np.any(np.abs(d_frac - np.rint(d_frac)) > tol_frac * 0.999):
        print("Warning: some atoms are close to half a cell from the reference; "
              "periodic alignment may be ambiguous.")
    d_frac -= np.rint(d_frac)
    return traj0 + d_frac @ cell0, ref_cells







#########################
# Calculations:

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

    V = np.dot(a1, np.cross(a2, a3))

    b1 = (2.0 * np.pi * np.cross(a2, a3) / V )
    b2 = (2.0 * np.pi * np.cross(a3, a1) / V )
    b3 = (2.0 * np.pi * np.cross(a1, a2) / V )

    return b1, b2, b3



def compute_full_peak_intensities(
        cells,
        coords,
        species,
        hkl_list,
        n_cells=(1, 1, 1),
        normalize=False,
        form_factor_file="Atomic_form_factors.dat",
        reference_positions=None,
        reference_cell=None,
        normalize_to_reference=False,
        use_displacements=False,
        coords_unwrapped=False):
    """
    Kinematic diffraction intensities from atomic positions and form factors.

    Parameters
    ----------
    cells, coords, species, hkl_list, n_cells, normalize, form_factor_file :
        As before.

    reference_positions : ndarray (natoms, 3) or None
        Reference structure (Cartesian, Angstrom), same atom order as coords.
        Align it to the trajectory first (e.g. load_reference_positions).

    reference_cell : ndarray (3, 3) or None
        Cell of the reference structure; defaults to cells[0].

    normalize_to_reference : bool
        If True (and normalize=True), I0 is the intensity of the reference
        structure instead of frame 0.

    use_displacements : bool
        If True, each frame uses  r_eff(t) = r_ref(t) + [r(t) - r(0)],
        i.e. the displacements from frame 0 are added to the reference
        positions. r_ref(t) is the reference scaled with the cell of frame t.

    coords_unwrapped : bool
        Set True if coords are already unwrapped. Otherwise they are
        unwrapped internally when use_displacements=True.
    """
    # ---------------------------------------------
    # Form-factor coefficients per atom
    # ---------------------------------------------
    ff_data = np.loadtxt(form_factor_file, skiprows=1)
    atomic_data = load_atomic_data()

    atomic_numbers, coefficients = [], []
    for symbol in species:
        element = get_element(atomic_data, name=symbol)
        if element is None:
            raise ValueError(f"Element '{symbol}' not found in INPUT_atomic_data.dat")
        Z = element["number"]
        if Z > ff_data.shape[0]:
            raise ValueError(f"No form-factor coefficients available for {symbol} (Z={Z})")
        atomic_numbers.append(Z)
        coefficients.append(ff_data[Z - 1, :])
    atomic_numbers = np.asarray(atomic_numbers)
    coefficients = np.asarray(coefficients)

    nframes = len(coords)
    natoms = len(species)
    if len(cells) != nframes:
        raise ValueError("cells and coords must contain the same number of frames")

    n = np.asarray(n_cells, dtype=float)
    if n.shape != (3,) or np.any(n < 1) or not np.allclose(n, np.rint(n)):
        raise ValueError("n_cells must contain three positive integers")

    # ---------------------------------------------
    # Reference structure
    # ---------------------------------------------
    ref = None
    if reference_positions is not None:
        ref = np.asarray(reference_positions, dtype=float)
        if ref.shape != (natoms, 3):
            raise ValueError(f"reference_positions must have shape ({natoms}, 3), got {ref.shape}")
        ref_cell = np.asarray(cells[0] if reference_cell is None else reference_cell, dtype=float)
        ref_frac = ref @ np.linalg.inv(ref_cell)       # fractional coordinates

    if (use_displacements or normalize_to_reference) and ref is None:
        raise ValueError("use_displacements / normalize_to_reference require reference_positions")

    # ---------------------------------------------
    # Displacements from frame 0 (option 2)
    # ---------------------------------------------
    if use_displacements:
        traj = np.asarray(coords, dtype=float) if coords_unwrapped \
            else np.asarray(unwrap_trajectory(cells, coords), dtype=float)
        disp = traj - traj[0][None, :, :]
        # Reference expressed in the frame-0 cell; identical to ref when
        # reference_cell is None or equals cells[0]
        r_ref0 = ref_frac @ np.asarray(cells[0], dtype=float)

    # ---------------------------------------------
    # Intensities frame by frame
    # ---------------------------------------------
    intensities = {tuple(hkl): np.zeros(nframes) for hkl in hkl_list}

    for iframe in range(nframes):
        cell = np.asarray(cells[iframe], dtype=float)

        if use_displacements:
            #r = ref_frac @ cell + disp[iframe]
            r = r_ref0 + disp[iframe]
        else:
            r = np.asarray(coords[iframe], dtype=float)
            if r.shape != (natoms, 3):
                raise ValueError(f"Unexpected coordinate shape in frame {iframe}: {r.shape}")

        I_frame = _frame_intensities(r, cell, hkl_list, n, coefficients, atomic_numbers)
        for hkl, val in I_frame.items():
            intensities[hkl][iframe] = val

    # ---------------------------------------------
    # Normalization (option 1)
    # ---------------------------------------------
    if normalize:
        if normalize_to_reference:
            I_ref = _frame_intensities(ref, ref_cell, hkl_list, n, coefficients, atomic_numbers)
        for hkl in intensities:
            I0 = I_ref[hkl] if normalize_to_reference else intensities[hkl][0]
            if I0 > 0.0:
                intensities[hkl] /= I0
            else:
                intensities[hkl][:] = np.nan

    return intensities


def _frame_intensities(r, cell, hkl_list, n, coefficients, atomic_numbers):
    """Kinematic intensity |F|^2 of every hkl for one set of positions."""
    HBAR = 1.054571817e-34      # J s
    ANGSTROM_TO_M = 1.0e-10

    B_super = np.asarray(reciprocal_lattice(cell), dtype=float)
    B_unit = B_super * n[:, None]
    natoms = len(r)

    out = {}
    for hkl in hkl_list:
        G = np.asarray(hkl, dtype=float) @ B_unit
        q = HBAR * np.linalg.norm(G) / ANGSTROM_TO_M
        f = np.array([
            form_factor(q, coefficients[j], atomic_numbers[j])
            for j in range(natoms)
        ])
        F = np.sum(f * np.exp(1j * (r @ G)))
        out[tuple(hkl)] = np.abs(F) ** 2
    return out







def compute_full_peak_intensities_OLD(
        cells,
        coords,
        species,
        hkl_list,
        n_cells=(1, 1, 1),
        normalize=False,
        form_factor_file="Atomic_form_factors.dat"):
    """
    Calculate kinematic diffraction intensities from atomic positions
    and atomic form factors.

    Parameters
    ----------
    cells : list of ndarray, each (3,3)
        Supercell vectors in Angstrom; vectors are rows.

    coords : list of ndarray, each (natoms,3)
        Atomic Cartesian coordinates in Angstrom.

    species : ndarray (natoms,)
        Atomic symbols in the same order as coords.

    hkl_list : list of (h,k,l)
        Miller indices with respect to the unit cell.

    n_cells : tuple of 3 ints
        Number of unit cells along each supercell direction.

    normalize : bool
        If True, normalize each intensity to its frame-0 value.

    form_factor_file : str
        File containing the five fitted coefficients per element.

    Returns
    -------
    intensities : dict
        {(h,k,l): ndarray (nframes,)}
    """

    # Physical constants (SI)
    HBAR = 1.054571817e-34  # J s
    ANGSTROM_TO_M = 1.0e-10

    # ---------------------------------------------
    # Read atomic form-factor coefficients
    # ---------------------------------------------
    ff_data = np.loadtxt(form_factor_file, skiprows=1)

    # Map each atom to its atomic number and form-factor coefficients
    atomic_data = load_atomic_data()

    atomic_numbers = []
    coefficients = []

    for symbol in species:
        element = get_element(atomic_data, name=symbol)

        if element is None:
            raise ValueError(
                f"Element '{symbol}' not found in INPUT_atomic_data.dat"
            )

        Z = element["number"]

        if Z > ff_data.shape[0]:
            raise ValueError(
                f"No form-factor coefficients available for {symbol} "
                f"(Z={Z})"
            )

        atomic_numbers.append(Z)
        coefficients.append(ff_data[Z - 1, :])

    atomic_numbers = np.asarray(atomic_numbers)
    coefficients = np.asarray(coefficients)

    nframes = len(coords)
    natoms = len(species)

    if len(cells) != nframes:
        raise ValueError("cells and coords must contain the same number of frames")

    n = np.asarray(n_cells, dtype=float)

    if n.shape != (3,) or np.any(n < 1) or not np.allclose(n, np.rint(n)):
        raise ValueError("n_cells must contain three positive integers")

    intensities = {
        tuple(hkl): np.zeros(nframes)
        for hkl in hkl_list
    }

    # ---------------------------------------------
    # Calculate intensity for each frame and peak
    # ---------------------------------------------
    for iframe in range(nframes):

        cell = np.asarray(cells[iframe], dtype=float)

        # Reciprocal lattice vectors of the supercell [1/Angstrom]
        B_super = np.asarray(reciprocal_lattice(cell))

        # Reciprocal lattice vectors of the unit cell
        B_unit = B_super * n[:, None]

        # Coordinates of all atoms in this frame [Angstrom]
        r = np.asarray(coords[iframe], dtype=float)

        if r.shape != (natoms, 3):
            raise ValueError(
                f"Unexpected coordinate shape in frame {iframe}: {r.shape}"
            )

        for hkl in hkl_list:

            hkl_array = np.asarray(hkl, dtype=float)

            # Scattering vector Q [1/Angstrom]
            G = hkl_array @ B_unit
            Q = np.linalg.norm(G)

            # Momentum transfer q = hbar*Q [kg m/s]
            q = HBAR * Q / ANGSTROM_TO_M

            # Atomic form factors f_j(Q)
            f = np.array([
                form_factor(q, coefficients[j], atomic_numbers[j])
                for j in range(natoms)
            ])

            # Complex scattering amplitude
            phase = r @ G
            F = np.sum(f * np.exp(1j * phase))

            # Full coherent intensity
            intensities[tuple(hkl)][iframe] = np.abs(F)**2

    if normalize:
        for hkl in intensities:
            I0 = intensities[hkl][0]

            if I0 > 0.0:
                intensities[hkl] /= I0
            else:
                intensities[hkl][:] = np.nan

    return intensities



##########################
# Calculations of the form factors:
#FF = form_factor(q, a, Z)


def form_factor(q: float,
                a: list[float],
                Z: float) -> float:
    """
    Atomic form factor (Eq. 2.7, 2.8, 2.9).

    Parameters
    ----------
    q : float
        Momentum transfer.
    a : list of float
        Fitted coefficients (length 5).
    Z : float
        Atomic number.

    Returns
    -------
    FF : float
        Form factor value.
    """
    mc = g_me * g_cvel
    x = q / mc * m_Thom_factor   # Eq. (2.5)
    x2 = x * x
    x4 = x2 * x2
    demf = 1.0 + a[3] * x2 + a[4] * x4

    numr = (1.0 + a[0] * x2 + a[1] * x * x2 + a[2] * x4)

    if ( demf > 1e8 * numr ):
        fxz = 0.0
    else :
        #fxz = Z * (1.0 + a[0] * x2 + a[1] * x * x2 + a[2] * x4) / (demf * demf)  # Eq. (2.7)
        fxz = Z * numr / (demf * demf)  # Eq. (2.7)

    FF = fxz

    if Z > 10.0 and fxz < 2.0:
        al = g_alpha * (Z - m_five_sixteenth)   # Eq. (2.9)
        b = np.sqrt(1.0 - al * al)              # Eq. (2.9)
        CapQ = q / (2.0 * mc * al)              # Eq. (2.9)

        denom = b * CapQ * (1.0 + CapQ * CapQ) ** b

        #Fk = np.sin(2.0 * b * np.arctan(CapQ)) / ( b * CapQ * (1.0 + CapQ * CapQ) ** b )    # Eq. (2.8)
        if (denom > 1.0e10):
            Fk = 0.0
        else:
            Fk = np.sin(2.0 * b * np.arctan(CapQ)) / denom    # Eq. (2.8)

        if FF < Fk:
            FF = Fk
    return FF




##########################
# Saving output files:


def write_extxyz(
        filename,
        symbols,
        positions,
        cell=None,
        pbc=(True, True, True),
        comment_info=None,
        extra=None,
        wrap=False):
    """
    Write one frame in extended XYZ format.

    extra : dict {name: ndarray (N,)} of additional per-atom scalar properties.
    wrap  : if True, wrap positions into the cell (needs `cell`).
    """
    positions = np.asarray(positions, dtype=float)
    n = len(symbols)
    if positions.shape != (n, 3):
        raise ValueError(f"positions must have shape ({n}, 3), got {positions.shape}")

    if wrap:
        if cell is None:
            raise ValueError("wrap=True requires a cell")
        frac = positions @ np.linalg.inv(cell)
        frac -= np.floor(frac)
        positions = frac @ cell

    extra = extra or {}
    props = "species:S:1:pos:R:3" + "".join(f":{k}:R:1" for k in extra)

    parts = []
    if cell is not None:
        lat = " ".join(f"{x:.10f}" for x in np.asarray(cell, float).reshape(9))
        parts.append(f'Lattice="{lat}"')
    parts.append(f"Properties={props}")
    parts.append('pbc="' + " ".join("T" if p else "F" for p in pbc) + '"')
    for key, val in (comment_info or {}).items():
        parts.append(f"{key}={val}")

    extra_cols = [np.asarray(v, float) for v in extra.values()]

    filename = unique_filename(filename) # make sure it won't overwrite existng files
    with open(filename, "x") as f:        # "x" = fail if the file appeared in the meantime
        f.write(f"{n}\n{' '.join(parts)}\n")
        for i in range(n):
            line = f"{symbols[i]:<3s} " + " ".join(f"{x:16.10f}" for x in positions[i])
            line += "".join(f" {c[i]:16.10e}" for c in extra_cols)
            f.write(line + "\n")
    print(f"Saved thermal positions: {filename}")





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
        theta_thermal,
        time=None):

    avg_temperature = np.nanmean(temperature)
    avg_theta_initial = np.nanmean(theta_initial)
    avg_theta_thermal = np.nanmean(theta_thermal)


    if time is not None and len(time) == len(temperature):   # use defined time:
        x_label = "Time (fs)"
        first_col = np.asarray(time, dtype=float)
        first_fmt = "%.6e"
    else:
        x_label = "Frame"
        first_col = np.arange(len(temperature))
        first_fmt = "%8d"

    header = (
        f"Average_temperature = "
        f"{avg_temperature:.10e} K\n"
        f"Average_Debye_temperature_initial = "
        f"{avg_theta_initial:.10e} K\n"
        f"Average_Debye_temperature_thermal = "
        f"{avg_theta_thermal:.10e} K\n"
        f"{x_label} "
        f"Temperature_K "
        f"ThetaD_initial_K "
        f"ThetaD_thermal_K"
    )
    print(header)

    np.savetxt(
        filename,
        np.column_stack((first_col, temperature, theta_initial, theta_thermal)),
        fmt=[first_fmt, "%.10e", "%.10e", "%.10e"],
        header=header
    )

    print(
        f"Saved Debye temperatures in file: "
        f"{filename}"
    )


def unique_filename(filename):
    """
    Return `filename` if it does not exist yet; otherwise insert an index
    before the extension: name.ext -> name_001.ext, name_002.ext, ...
    """
    if not os.path.exists(filename):
        return filename
    base, ext = os.path.splitext(filename)
    i = 1
    while True:
        candidate = f"{base}_{i:03d}{ext}"
        if not os.path.exists(candidate):
            return candidate
        i += 1


##########################
# Plotting subroutines:

def plot_dw_peak_intensities(
        intensities,
        time=None,
        output_file="OUT_DW_peak_intensities.png",
        DW_title="Debye-Waller intensities",
        ylabel=r"$I/I_0$"):
    """
    Plot Debye-Waller peak intensities.
    """

    hkls = list(intensities.keys())
    nframes = len(intensities[hkls[0]])
    x, xlabel = get_x_axis(nframes, time)

    frames = np.arange(nframes)

    fig, ax = plt.subplots(figsize=(5, 4))

    styles = ["-", "--", ":", "-."]

    for i, hkl in enumerate(hkls):
        ax.plot(
            x, intensities[hkl],
            lw=1.0,
            ls=styles[i % len(styles)],
            label=f"({hkl[0]}{hkl[1]}{hkl[2]})"
        )

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_ylim(0.0, 1.05)
    ax.grid(False)
    ax.legend(frameon=True)
    ax.set_title(DW_title)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {output_file}")





def plot_debye_temperature(
        temperature,
        theta_initial,
        theta_thermal,
        output_file="OUT_Debye_temperature.png",
        x_label = "Frame"):
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
        lw=1.0,
        color="tab:blue",
        label=r"$\Theta_D$ from initial MSD"
    )

    ax.plot(
        frames,
        theta_thermal,
        lw=1.0,
        ls="--",
        color="tab:red",
        label=r"$\Theta_D$ from thermal MSD"
    )

    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(r"Debye temperature (K)", fontsize=12)

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
        frameon=True
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
             output_file="OUT_MSDs.png", x_label="Frame"):
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
        lw=1.0,
        label="MSD from frame 0",
    )

    ax.plot(
        frames,
        msd_avg,
        lw=1.0,
        ls="--",
        label="Time-origin averaged MSD",
    )

    ax.plot(
        frames,
        thermal_msd,
        lw=1.0,
        ls=":",
        label="Thermal MSD",
    )

    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(r"MSD ($\AA^2$)", fontsize=12)
    ax.grid(False)
    ax.legend(fontsize=10)
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





def get_atomic_form_factors(targetpar, filename: str = "Atomic_form_factors.dat") -> None:
    """
    Read atomic form factor parameters (a1..a5) from database and attach them
    to each element in targetpar.element_data.

    Parameters
    ----------
    targetpar : TargetPar
        Target parameters with attributes:
        - elements : list of element symbols
        - element_data : dict mapping symbol -> element dict
    filename : str
        Path to the atomic form factors database file.
    """

    # Resolve path depending on frozen/executable vs source
    #if getattr(sys, 'frozen', False):
    #    project_root = Path(sys.executable).parent
    #else:
    #    project_root = Path(__file__).resolve().parent.parent
    #file_path = project_root / "INPUT" / "Atomic_data" / filename
    file_path = Path(filename)

    # --- Check file existence ---
    if not file_path.exists():
        print(f"Atomic form factors file not found: {file_path}")
        #sys.exit(1)

    # --- Load numeric data (rows = elements, cols = 5 coefficients) ---
    # Skip the first line (header with column names)
    data = np.loadtxt(file_path, skiprows=1)  # shape (N_elements, 5)

    # --- Distribute to all elements ---
    for elem in targetpar.elements:
        edata = targetpar.element_data[elem]
        Z = edata["number"]

        if Z > data.shape[0]:
            print(f"Atomic form factors: Element #{Z} not found in {filename}")
            sys.exit(1)

        # Store coefficients as a list or NumPy array
        edata["Form_factor"] = data[Z-1, :].tolist()   # Fortran arrays are 1-based



##########################
if __name__ == "__main__":
    main()
