import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- Parse powder file (multi-block format) ---
def parse_powder_file(path):
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
                    col_names = content.split()
                    continue
                if current_time is not None and current_rows:
                    blocks.append((current_time, np.array(current_rows)))
                current_time = t
                current_rows = []
            else:
                current_rows.append([float(x) for x in line.split()])
    if current_time is not None and current_rows:
        blocks.append((current_time, np.array(current_rows)))
    return col_names, blocks

# --- Load fitted results (automatic peaks) ---
df = pd.read_csv("OUTPUT_diffraction_powder_fit_automatic.csv")
df = df.sort_values("time_fs")

# --- Load raw powder data ---
col_names, blocks = parse_powder_file("OUTPUT_diffraction_powder.dat")

# --- Video of evolving peaks with raw data ---
def animate_peaks_with_raw(df, blocks, out_file="peak_evolution_with_raw_auto.mp4"):
    fig, ax = plt.subplots(figsize=(8, 6))
    times = sorted(df["time_fs"].unique())
    hkl_list = df["hkl"].unique()
    fit_lines = {hkl: ax.plot([], [], label=f"fit {hkl}")[0] for hkl in hkl_list}
    raw_line, = ax.plot([], [], "o", ms=2, color="0.5", label="raw data")
    ax.set_xlim(0, 180)
    ax.set_ylim(0, df["amplitude"].max() * 1.5)
    ax.set_xlabel("2θ (deg)")
    ax.set_ylabel("Intensity (arb. units)")
    ax.legend(fontsize=8)

    def init():
        raw_line.set_data([], [])
        for line in fit_lines.values():
            line.set_data([], [])
        return [raw_line] + list(fit_lines.values())

    def update(frame):
        t = times[frame]
        ax.set_title(f"Simulation time {t:.2f} fs")
        block = next(b for b in blocks if abs(b[0] - t) < 1e-6)
        two_theta = block[1][:, col_names.index("Angle")]
        intensity = block[1][:, col_names.index("Total")]
        raw_line.set_data(two_theta, intensity)
        for line in fit_lines.values():
            line.set_data([], [])
        subset = df[df["time_fs"] == t]
        xx = np.linspace(0, 180, 800)
        for _, row in subset.iterrows():
            yy = row["amplitude"] * (
                row["eta"] / (1.0 + 4.0 * (xx - row["center_deg"])**2 / row["fwhm_deg"]**2)
                + (1.0 - row["eta"]) * np.exp(-4*np.log(2)*(xx - row["center_deg"])**2 / row["fwhm_deg"]**2)
            )
            fit_lines[row["hkl"]].set_data(xx, yy)
        return [raw_line] + list(fit_lines.values())

    ani = animation.FuncAnimation(fig, update, frames=len(times),
                                  init_func=init, blit=True, repeat=False)
    ani.save(out_file, writer="ffmpeg", dpi=150)
    plt.close(fig)

# --- Integrated areas vs time ---
def plot_integrated_areas(df, out_file="peak_integrated_areas_vs_time_auto.png"):
    fig, ax = plt.subplots(figsize=(5, 4))
    for hkl in df["hkl"].unique():
        subset = df[df["hkl"] == hkl]
        ax.plot(subset["time_fs"], subset["integrated_area"], label=hkl)
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Integrated Area (arb. units)")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_file, dpi=200)
    plt.close(fig)

# --- Normalized integrated areas vs time ---
def plot_normalized_areas(df, out_file="peak_normalized_areas_vs_time_auto.png",
                          min_area_threshold=1e-3):
    fig, ax = plt.subplots(figsize=(5, 4))
    for hkl in df["hkl"].unique():
        subset = df[df["hkl"] == hkl].copy().sort_values("time_fs")
        initial_area = subset["integrated_area"].iloc[0]
        subset["normalized_area"] = subset["integrated_area"] / initial_area if initial_area != 0 else 0.0
        lw = 0.2 if initial_area < min_area_threshold else 1.5
        alpha = 0.6 if initial_area < min_area_threshold else 1.0
        ax.plot(subset["time_fs"], subset["normalized_area"], label=hkl, lw=lw, alpha=alpha)
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Normalized Integrated Area (initial=1.0)")
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_file, dpi=200)
    plt.close(fig)

# --- Run all ---
if __name__ == "__main__":
    plot_integrated_areas(df)
    plot_normalized_areas(df)
    animate_peaks_with_raw(df, blocks)
    print("✅ Video with raw data, integrated area plot, and normalized area plot saved.")
