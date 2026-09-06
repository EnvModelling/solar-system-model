#!/usr/bin/env python3
"""Compare heliocentric Sun-only and all-interaction SSM runs."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ssm_io import BODY_NAMES, load_output


def sun_distance(data):
    rel = data["pos"] - data["pos"][:, 0:1, :]
    return np.linalg.norm(rel[:, 1:, :], axis=2)


def main():
    parser = ArgumentParser()
    parser.add_argument("keplerian")
    parser.add_argument("interacting")
    parser.add_argument("output")
    args = parser.parse_args()

    kep = load_output(args.keplerian)
    full = load_output(args.interacting)

    r_kep = sun_distance(kep)
    r_full = sun_distance(full)

    t = kep["years"]
    delta = np.empty_like(r_kep)

    for i in range(9):
        delta[:, i] = (
            np.interp(t, full["years"], r_full[:, i]) - r_kep[:, i]
        ) / 1000.0

    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True)
    axes = axes.ravel()

    for i, ax in enumerate(axes, start=1):
        ax.plot(t, delta[:, i - 1], "k-", lw=0.45)
        ax.set_xlabel("Earth years")
        ax.set_ylabel(r"$\Delta r$ (km)")
        ax.set_title(BODY_NAMES[i])

    fig.suptitle("Effect of mutual planetary interactions: all interactions - Sun only")
    fig.tight_layout()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
