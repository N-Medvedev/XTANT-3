#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse


def moving_average(y, window):
    """
    Compute moving average using a centered convolution.
    """
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="valid")


def main():
    parser = argparse.ArgumentParser(
        description="Plot force-velocity correlator with moving average."
    )
    parser.add_argument(
        "-w", "--window",
        type=int,
        default=20,
        help="Moving average window size (default: 20)"
    )
    parser.add_argument(
        "-i", "--input",
        default="OUTPUT_force_velocity_correlator.dat",
        help="Input file"
    )

    args = parser.parse_args()

    # Load file, skipping comments beginning with '#'
    data = np.loadtxt(args.input, comments="#")

    x = data[:, 0]
    y = data[:, 1]

    # Moving average
    y_ma = moving_average(y, args.window)

    # Corresponding x values for the averaged data
    offset = (args.window - 1) // 2
    x_ma = x[offset:offset + len(y_ma)]

    plt.figure(figsize=(5, 4))

    plt.plot(x, y, lw=1, alpha=0.5,
             label="Raw data")

    plt.plot(x_ma, y_ma,
             lw=1.5, color="red",
             label=f"{args.window}-point moving average")
    
    
    plt.axhline(
    y=0.0,
    color='black',
    linestyle='-',
    linewidth=0.5,
    alpha=0.7,
    zorder=0)

    plt.xlabel("Time (fs)", fontsize=12)
    plt.ylabel("Force-velocity correlator (eV/fs)", fontsize=12)
    plt.title(r"$\langle \mathrm{Force} \times \mathrm{velocity} \rangle$")
    plt.legend()
    plt.grid(False)
    plt.tight_layout()

    plt.savefig("OUT_Fv_moving_average.png", dpi=300, bbox_inches="tight")
    #plt.show()


if __name__ == "__main__":
    main()