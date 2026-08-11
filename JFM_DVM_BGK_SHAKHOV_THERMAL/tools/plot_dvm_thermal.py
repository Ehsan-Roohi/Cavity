#!/usr/bin/env python3
"""Unfiltered plots and parity diagnostics for reconstructed thermal DVM data."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parity_errors(a, u, v):
    return {
        "scalar_x_even_max_abs": float(np.max(np.abs(a - a[:, ::-1]))),
        "scalar_y_even_max_abs": float(np.max(np.abs(a - a[::-1, :]))),
        "u_x_odd_max_abs": float(np.max(np.abs(u + u[:, ::-1]))),
        "u_y_even_max_abs": float(np.max(np.abs(u - u[::-1, :]))),
        "v_x_even_max_abs": float(np.max(np.abs(v - v[:, ::-1]))),
        "v_y_odd_max_abs": float(np.max(np.abs(v + v[::-1, :]))),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    raw = np.genfromtxt(args.csv, delimiter=",", names=True)
    xvals = np.unique(raw["x"])
    yvals = np.unique(raw["y"])
    nx, ny = xvals.size, yvals.size
    shape = (ny, nx)
    rho = raw["rho"].reshape(shape)
    u = raw["u"].reshape(shape)
    v = raw["v"].reshape(shape)
    temp = raw["T"].reshape(shape)
    speed = np.hypot(u, v)

    with open(args.metadata, encoding="utf-8") as handle:
        metadata = json.load(handle)
    diagnostics = parity_errors(temp, u, v)
    diagnostics.update({
        "max_speed_nondimensional_from_csv": float(np.max(speed)),
        "all_fields_unfiltered": True,
        "metadata_converged": bool(metadata["converged"]),
        "metadata_mass_drift_relative": float(metadata["mass_drift_relative"]),
        "metadata_max_wall_mass_flux": float(metadata["max_physical_wall_mass_flux"]),
        "metadata_max_target_conservation_error": float(
            metadata["max_discrete_target_conservation_error"]
        ),
        "metadata_positive_projection_failures": int(
            metadata["positive_projection_failures"]
        ),
        "metadata_max_projection_iterations_used": int(
            metadata["max_projection_iterations_used"]
        ),
        "metadata_negative_target_values": int(metadata["negative_target_values"]),
        "metadata_clipped_updated_values": int(metadata["clipped_updated_values"]),
        "metadata_shakhov_limited_target_values": metadata.get(
            "shakhov_limited_target_values"
        ),
        "metadata_shakhov_limited_H_base_mass_fraction": metadata.get(
            "shakhov_limited_H_base_mass_fraction"
        ),
        "metadata_shakhov_limited_B_base_content_fraction": metadata.get(
            "shakhov_limited_B_base_content_fraction"
        ),
    })

    with open(args.output_prefix + "_postcheck.json", "w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2)

    xx, yy = np.meshgrid(xvals, yvals)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    fields = ((temp, r"$T/T_h$"), (rho, r"$\rho/\rho_0$"), (speed, r"$|\mathbf{u}|/U_{mp,h}$"))
    for ax, (field, title) in zip(axes, fields):
        cf = ax.contourf(xx, yy, field, 40, cmap="viridis")
        fig.colorbar(cf, ax=ax)
        ax.set_title(title)
        ax.set_xlabel(r"$x/L$")
        ax.set_ylabel(r"$y/L$")
        ax.set_aspect("equal")
        ax.plot([-0.5, 0.5, 0.5, -0.5, -0.5],
                [-0.5, -0.5, 0.5, 0.5, -0.5], "k-", lw=1.2)
    if np.max(speed) > 0.0:
        axes[2].streamplot(xvals, yvals, u, v, color="white", density=1.6,
                           linewidth=0.7, arrowsize=0.8)
    fig.suptitle(f"DVM {metadata['model']}, Kn={metadata['Kn_paper']:g}, RT={metadata['RT']:g}")
    fig.savefig(args.output_prefix + "_fields.png", dpi=220)
    plt.close(fig)

    ic = int(np.argmin(np.abs(xvals)))
    jc = int(np.argmin(np.abs(yvals)))
    with open(args.output_prefix + "_centerlines.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("coordinate", "u_at_x_near_0", "v_at_y_near_0"))
        n = max(nx, ny)
        for k in range(n):
            coord = yvals[k] if k < ny else ""
            uu = u[k, ic] if k < ny else ""
            vv = v[jc, k] if k < nx else ""
            writer.writerow((coord, uu, vv))

    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
