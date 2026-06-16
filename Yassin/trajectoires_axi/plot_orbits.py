#!/home/yassin/agama_gpu_env/bin/python
"""Plotter for the fine-step axisymmetric integration.

Reads the dedicated trajectory file
(``trajectories_axi_25Myr_step0p5.npz``) and draws the cluster/maser
orbits in the three planes XY, XZ and R-Z.  

Usage
-----
    plot_reintegrated_orbits.py
    plot_reintegrated_orbits.py --tmax 10
    plot_reintegrated_orbits.py --npz <path>
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

BASE_DIR = ""
DEFAULT_NPZ = os.path.join(BASE_DIR, "outputs", "trajectories_axi_25Myr_step0p5.npz")
OUT_FIG_DIR = os.path.join(BASE_DIR, "figures", "axi", "orbits")

# Object groupings / colours 
GROUP_5MYR = ["Group A", "Group D", "Group F"]
GROUP_10MYR = ["Group B", "Group C", "Group E"]
MASERS = ["W75N", "DR21", "DR20", "IRAS20290+4052"]
GROUP_SPECS = [
    (GROUP_5MYR, "tab:red", "Group ADF"),
    (GROUP_10MYR, "deepskyblue", "Group BCE"),
    (MASERS, "tab:green", "Masers"),
]

PLANES = {
    "xy": ("X [kpc]", "Y [kpc]", True),
    "xz": ("X [kpc]", "Z [kpc]", False),
    "Rz": ("R [kpc]", "Z [kpc]", False),
}


def plane_coords(traj_j, plane):
    """Return (a, b) coordinate arrays for one object in the requested plane."""
    x, y, z = traj_j[:, 0], traj_j[:, 1], traj_j[:, 2]
    if plane == "xy":
        return x, y
    if plane == "xz":
        return x, z
    if plane == "Rz":
        return np.hypot(x, y), z
    raise ValueError(f"unknown plane {plane!r}")


def legend_handles(t_max_myr):
    handles = [
        mlines.Line2D([], [], color=c, marker="o", lw=1.2, mec="k", mew=0.5, label=lab)
        for _, c, lab in GROUP_SPECS
    ]
    handles.append(mlines.Line2D([], [], color="gray", marker="o", lw=0, mec="k",
                                 mew=0.5, label="t = 0 (today)"))
    handles.append(mlines.Line2D([], [], color="gray", marker="X", lw=0, mec="k",
                                 mew=0.5, ms=11, label=f"{int(-t_max_myr)} Myr"))
    return handles


def plot_plane(t_arr, traj, names, plane, t_max_myr, fname):
    i_end = int(np.argmin(np.abs(t_arr - (-t_max_myr))))
    xlabel, ylabel, equal = PLANES[plane]
    fig, ax = plt.subplots(figsize=(9, 8))
    for names_list, color, _ in GROUP_SPECS:
        for nm in names_list:
            if nm not in names:
                continue
            j = names.index(nm)
            a, b = plane_coords(traj[j], plane)
            ax.plot(a[i_end:], b[i_end:], "-", color=color, lw=1.0, alpha=0.6)
            ax.plot(a[-1], b[-1], "o", color=color, ms=7, mec="k", mew=0.5)       # today
            ax.plot(a[i_end], b[i_end], "X", color=color, ms=11, mec="k", mew=0.5)  # -t_max
            ax.annotate(nm, (a[-1], b[-1]), xytext=(5, 4),
                        textcoords="offset points", fontsize=8, color=color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    if equal:
        ax.set_aspect("equal")
    ax.set_title(f"Axisymmetric integration ({int(t_max_myr)} Myr) - {plane.upper()} plane")
    ax.legend(handles=legend_handles(t_max_myr), loc="upper center",
              bbox_to_anchor=(0.5, 1.10), ncol=5, fontsize=9, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {fname}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--npz", default=DEFAULT_NPZ, help="trajectory .npz file")
    ap.add_argument("--tmax", type=float, default=None,
                    help="window in Myr to display (default: full file)")
    ap.add_argument("--outdir", default=OUT_FIG_DIR, help="output figure folder")
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    t_arr = np.asarray(d["times_Myr"], dtype=float)
    traj = np.asarray(d["traj"], dtype=float)
    names = [str(n) for n in d["names"]]
    t_max = args.tmax if args.tmax is not None else abs(float(t_arr.min()))

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Loaded {args.npz}: traj{traj.shape}, t in [{t_arr.min():.1f}, {t_arr.max():.1f}] Myr")
    for plane in ("xy", "xz", "Rz"):
        plot_plane(t_arr, traj, names, plane, t_max,
                   os.path.join(args.outdir, f"reintegrated_{int(t_max)}Myr_{plane}.png"))


if __name__ == "__main__":
    main()
