#!/usr/bin/env python3
"""Animate SSM output with thin black orbital paths and moving bodies."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from ssm_io import AU, BODY_NAMES, load_output, shift_origin


def main():
    parser = ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--system", choices=("inner", "all"), default="all")
    parser.add_argument(
        "--origin",
        choices=("sun", "barycentre", "inertial"),
        default="sun",
    )
    parser.add_argument("--frames", type=int, default=300)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--max-years", type=float)
    parser.add_argument("--dpi", type=int, default=110)
    args = parser.parse_args()

    data = load_output(args.input)
    years = data["years"]
    pos = shift_origin(data["pos"], data["gm"], args.origin) / AU

    if args.max_years is not None:
        mask = years <= args.max_years
        if np.count_nonzero(mask) < 2:
            raise ValueError("max-years leaves fewer than two time samples")
        years = years[mask]
        pos = pos[mask]

    bodies = list(range(5)) if args.system == "inner" else list(range(10))

    nframes = min(max(args.frames, 2), len(years))
    frame_indices = np.unique(
        np.linspace(0, len(years) - 1, nframes).astype(int)
    )

    # Axes are based on the whole selected trajectory, so the animation does not
    # rescale from frame to frame.
    xy = pos[:, bodies, :2]
    lim = np.max(np.abs(xy))
    if not np.isfinite(lim) or lim <= 0:
        lim = 1.0
    lim *= 1.05

    fig, ax = plt.subplots(figsize=(8, 8))

    # Thin black lines show the orbital paths over the displayed time interval.
    for body in bodies:
        ax.plot(
            pos[:, body, 0],
            pos[:, body, 1],
            "k-",
            lw=0.45,
            alpha=0.45,
        )

    markers = []
    for body in bodies:
        marker, = ax.plot(
            [],
            [],
            marker="o",
            linestyle="none",
            ms=5 if body else 8,
            label=BODY_NAMES[body],
        )
        markers.append(marker)

    time_text = ax.text(
        0.02,
        0.98,
        "",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (AU)")
    ax.set_ylabel("y (AU)")
    ax.set_title(f"Solar System ({args.origin} frame)")
    ax.legend(loc="upper right", fontsize=8, ncol=2)

    def update(frame_number):
        idx = frame_indices[frame_number]
        for marker, body in zip(markers, bodies):
            marker.set_data([pos[idx, body, 0]], [pos[idx, body, 1]])
        time_text.set_text(f"{years[idx]:.2f} years")
        return [*markers, time_text]

    animation = FuncAnimation(
        fig,
        update,
        frames=len(frame_indices),
        interval=1000.0 / args.fps,
        blit=True,
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    animation.save(output, writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
