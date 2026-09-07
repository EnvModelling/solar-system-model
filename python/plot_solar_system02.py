import os
import getpass
import matplotlib
matplotlib.use("agg")

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

username = getpass.getuser()
nc = NetCDFFile("/tmp/" + username + "/output.nc")

objs = (
    "Mercury", "Venus", "Earth", "Mars",
    "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"
)

# Approximate orbital periods in Earth years
periods = np.array([
    0.241, 0.615, 1.000, 1.881,
    11.86, 29.46, 84.01, 164.8, 248.0
])

# Convert time to Earth years
time_years = nc.variables["time"][:] / (365.25 * 86400.0)

# Heliocentric distances of planets from Sun
AU = 1.495978707e11
distance = np.sqrt(np.sum(nc.variables["pos"][:, 1:10, :]**2, axis=2)) / AU

def format_years(val):
    if val >= 100:
        return f"{val:.0f} y"
    elif val >= 10:
        return f"{val:.1f} y"
    else:
        return f"{val:.2f} y"

fig, axes = plt.subplots(3, 3, figsize=(12, 10))

for i, ax in enumerate(axes.flat):
    # Main full-timescale plot
    ax.plot(time_years, distance[:, i], "k-", lw=0.25)
    ax.set_title(objs[i], fontsize=10)
    ax.set_xlabel("Time (years)", fontsize=8)
    ax.set_ylabel("Distance from Sun (AU)", fontsize=8)
    ax.tick_params(labelsize=8)

    # Inset showing first 5 orbital periods
    inset = inset_axes(ax, width="42%", height="42%", loc="upper right")

    tmax = min(5.0 * periods[i], time_years[-1])
    use = time_years <= tmax

    inset.plot(time_years[use], distance[use, i], "k-", lw=0.7)
    inset.set_xlim(0.0, tmax)

    # Red inset axes and text
    for spine in inset.spines.values():
        spine.set_color("red")
        spine.set_linewidth(1.0)

    inset.tick_params(axis="both", colors="red", labelsize=6)
    inset.set_title("5 periods", fontsize=7, color="red")

    # Red text on the main panel below the inset
    ax.text(
        0.56, 0.43,
        "5 periods = " + format_years(tmax),
        transform=ax.transAxes,
        color="red",
        fontsize=7,
        ha="left",
        va="top"
    )

fig.suptitle(
    "Planet–Sun distance over 10,000 years\n(with inset showing first 5 orbital periods)",
    fontsize=14
)

fig.subplots_adjust(left=0.07, right=0.98, bottom=0.07, top=0.90,
                    wspace=0.35, hspace=0.45)

outfile = "/tmp/" + username + "/milankovitch.png"
fig.savefig(outfile, dpi=300)
plt.close(fig)

nc.close()

print("Saved:", outfile)
