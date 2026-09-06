#!/usr/bin/env python3
"""Small I/O helpers shared by the SSM plotting scripts."""

from pathlib import Path
import numpy as np
from netCDF4 import Dataset

SECONDS_PER_YEAR = 365.25 * 86400.0
AU = 1.495978707e11

BODY_NAMES = (
    "Sun", "Mercury", "Venus", "Earth", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto",
)

# Used for older SSM output files that pre-date the GM variable.
DEFAULT_GM = np.array(
    [
        1.327124400e11,
        22032.09,
        324858.63,
        398600.440,
        42828.3,
        126686511.0,
        37931207.8,
        5793966.0,
        6835107.0,
        872.4,
    ],
    dtype=float,
) * 1.0e9


def _get_variable(nc, *names):
    for name in names:
        if name in nc.variables:
            return nc.variables[name]
    raise KeyError(f"None of these variables are present: {', '.join(names)}")


def _normalise_vector_array(array, ntime):
    """Return a vector variable as [time, body, xyz]."""
    a = np.asarray(array)

    if a.ndim != 3:
        raise ValueError(f"Expected a 3-D position/velocity array, got {a.shape}")

    # Current SSM files, and the historical Python scripts, use [time, body, xyz].
    if a.shape[0] == ntime and a.shape[-1] == 3:
        return a

    # Be tolerant of files exposed in the Fortran declaration order.
    time_axes = [i for i, n in enumerate(a.shape) if n == ntime]
    xyz_axes = [i for i, n in enumerate(a.shape) if n == 3]

    if not time_axes or not xyz_axes:
        raise ValueError(f"Cannot identify time/xyz axes in array with shape {a.shape}")

    time_axis = time_axes[0]
    xyz_axis = next((i for i in xyz_axes if i != time_axis), xyz_axes[0])
    body_axis = next(i for i in range(3) if i not in (time_axis, xyz_axis))
    return np.transpose(a, (time_axis, body_axis, xyz_axis))


def load_output(filename):
    """Load SSM output and return a dictionary of NumPy arrays and metadata."""
    filename = Path(filename)

    with Dataset(filename) as nc:
        time = np.asarray(_get_variable(nc, "time", "TIME")[:], dtype=float)
        pos = _normalise_vector_array(_get_variable(nc, "pos", "POS")[:], len(time))
        vel = _normalise_vector_array(_get_variable(nc, "vel", "VEL")[:], len(time))

        if "gm" in nc.variables:
            gm = np.asarray(nc.variables["gm"][:], dtype=float)
        else:
            gm = DEFAULT_GM.copy()

        if "gravity_source" in nc.variables:
            gravity_source = np.asarray(nc.variables["gravity_source"][:], dtype=int)
        else:
            gravity_source = None

        attrs = {name: nc.getncattr(name) for name in nc.ncattrs()}

    return {
        "time": time,
        "years": time / SECONDS_PER_YEAR,
        "pos": pos.astype(float, copy=False),
        "vel": vel.astype(float, copy=False),
        "gm": gm,
        "gravity_source": gravity_source,
        "attrs": attrs,
    }


def shift_origin(pos, gm, origin):
    """Return positions in an inertial, Sun-centred or barycentric frame."""
    origin = origin.lower()
    p = np.asarray(pos, dtype=float)

    if origin == "inertial":
        return p.copy()

    if origin == "sun":
        return p - p[:, 0:1, :]

    if origin == "barycentre":
        weights = np.asarray(gm, dtype=float)
        bary = np.sum(p * weights[None, :, None], axis=1) / np.sum(weights)
        return p - bary[:, None, :]

    raise ValueError(f"Unknown origin {origin!r}")


def choose_stride(npoints, max_points):
    if max_points <= 0:
        return 1
    return max(1, int(np.ceil(npoints / max_points)))
