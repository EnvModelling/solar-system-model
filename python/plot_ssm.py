#!/usr/bin/env python3
"""Generate static figures from one SSM NetCDF output file."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ssm_io import AU, BODY_NAMES, choose_stride, load_output, shift_origin


def plot_orbits(data, outpath, origin, inner, max_points):
    pos = shift_origin(data["pos"], data["gm"], origin) / AU
    stride = choose_stride(len(pos), max_points)

    if inner:
        bodies = range(0, 5)
        title = f"Inner Solar System ({origin} frame)"
    else:
        bodies = range(0, 10)
        title = f"Solar System orbits ({origin} frame)"

    fig, ax = plt.subplots(figsize=(8, 8))

    for body in bodies:
        xy = pos[::stride, body, :2]
        ax.plot(xy[:, 0], xy[:, 1], "k-", lw=0.45, alpha=0.65)
        ax.plot(
            pos[-1, body, 0],
            pos[-1, body, 1],
            marker="o",
            linestyle="none",
            ms=5 if body else 8,
            label=BODY_NAMES[body],
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (AU)")
    ax.set_ylabel("y (AU)")
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_distances(data, outpath, max_points):
    years = data["years"]
    pos = data["pos"] - data["pos"][:, 0:1, :]
    radius = np.linalg.norm(pos[:, 1:, :], axis=2)
    stride = choose_stride(len(years), max_points)

    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True)
    axes = axes.ravel()

    for i, ax in enumerate(axes, start=1):
        ax.plot(
            years[::stride],
            radius[::stride, i - 1],
            "k-",
            lw=0.35,
            alpha=0.8,
            rasterized=True,
        )
        ax.set_xlabel("Earth years")
        ax.set_ylabel(f"{BODY_NAMES[i]} - Sun distance")

    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_sun_wobble(data, outpath, max_points):
    # Recentring each time on the barycentre removes any bulk centre-of-mass drift
    # from heliocentric initial state data.
    bary = shift_origin(data["pos"], data["gm"], "barycentre")
    sun = bary[:, 0, :] / 1000.0  # km
    stride = choose_stride(len(sun), max_points)

    amplitude = np.max(np.linalg.norm(sun[:, :2], axis=1))
    if amplitude < 1.0e-9:
        print("Sun is fixed in this output; skipping barycentric wobble figure.")
        return

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(sun[::stride, 0], sun[::stride, 1], "k-", lw=0.6)
    ax.plot(sun[-1, 0], sun[-1, 1], marker="o", linestyle="none", ms=6)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x relative to barycentre (km)")
    ax.set_ylabel("y relative to barycentre (km)")
    ax.set_title("Solar wobble about the system barycentre")
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def main():
    parser = ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--outdir", default="docs/figures")
    parser.add_argument("--prefix", default="")
    parser.add_argument(
        "--origin",
        choices=("sun", "barycentre", "inertial"),
        default="sun",
        help="Frame used for orbit figures",
    )
    parser.add_argument("--max-points", type=int, default=20000)
    args = parser.parse_args()

    data = load_output(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = f"{args.prefix}_" if args.prefix else ""

    plot_orbits(
        data,
        outdir / f"{prefix}orbits_all.png",
        args.origin,
        inner=False,
        max_points=args.max_points,
    )
    plot_orbits(
        data,
        outdir / f"{prefix}orbits_inner.png",
        args.origin,
        inner=True,
        max_points=args.max_points,
    )
    plot_distances(
        data,
        outdir / f"{prefix}distance_to_sun.png",
        max_points=args.max_points,
    )

    if args.origin == "barycentre":
        plot_sun_wobble(
            data,
            outdir / f"{prefix}sun_wobble.png",
            max_points=args.max_points,
        )


if __name__ == "__main__":
    main()
