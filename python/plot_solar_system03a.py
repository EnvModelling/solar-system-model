import os
import getpass

import matplotlib
matplotlib.use('agg')
from matplotlib import rc

rc('font', family='serif')
rc('text', usetex=True)

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt

username = getpass.getuser()

nc = NetCDFFile('/tmp/' + username + '/output.nc')

# ------------------------------------------------------------
# Fourier analysis
# ------------------------------------------------------------

n_use = min(50000, len(nc.variables['time']))

dt = (
    nc.variables['time'][1] - nc.variables['time'][0]
) / (86400.0 * 365.25)   # years

NFFT = int(2.0 ** np.ceil(np.log2(n_use)))

freq = np.fft.rfftfreq(NFFT, d=dt)
freq = freq[1:]                 # remove zero frequency
period = 1.0 / freq

names = (
    'Mercury',
    'Venus',
    'Earth',
    'Mars',
    'Jupiter',
    'Saturn',
    'Uranus',
    'Neptune',
    'Pluto'
)

window = np.hanning(n_use)

AU = 1.495978707e11

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, axes = plt.subplots(
    3, 3,
    figsize=(12, 9),
    sharex=True
)

for i, ax in enumerate(axes.flat, start=1):

    # True Sun-planet distance using x, y, z
    radius = np.sqrt(
        np.sum(
            nc.variables['pos'][:n_use, i, 0:3].astype(float)**2,
            axis=1
        )
    )

    # Convert to AU
    radius = radius / AU

    # Remove mean distance so we analyse the oscillatory part
    signal = radius - np.mean(radius)

    # Window and FFT
    fft = np.fft.rfft(signal * window, n=NFFT)

    # One-sided Fourier amplitude in AU
    amplitude = 2.0 * np.abs(fft) / np.sum(window)
    amplitude = amplitude[1:]

    ax.plot(period, amplitude, 'k-', linewidth=0.8)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(0.1, 1000.0)
    ax.set_title(names[i - 1], fontsize=11)

    ax.grid(
        which='major',
        linestyle=':',
        linewidth=0.5,
        alpha=0.45
    )
    ax.grid(
        which='minor',
        linestyle=':',
        linewidth=0.25,
        alpha=0.15
    )

    ax.tick_params(labelsize=8)

for ax in axes[-1, :]:
    ax.set_xlabel(r'time period (years)', fontsize=10)

for ax in axes[:, 0]:
    ax.set_ylabel(r'Fourier amplitude (AU)', fontsize=10)

fig.suptitle(
    r'Fourier analysis of planet--Sun distance',
    fontsize=14
)

fig.subplots_adjust(
    left=0.09,
    right=0.98,
    bottom=0.08,
    top=0.93,
    wspace=0.18,
    hspace=0.25
)

nc.close()

outfile = '/tmp/' + username + '/fourier.png'
fig.savefig(outfile, dpi=300)
plt.close(fig)

print('Saved:', outfile)
