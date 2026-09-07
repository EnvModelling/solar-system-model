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
) / (86400.0 * 365.25)       # years

NFFT = int(2.0 ** np.ceil(np.log2(n_use)))

# Gravitational parameters (GM): proportional to mass
mu = np.array([
    1.327124400e11,   # Sun
    22032.09,         # Mercury
    324858.63,        # Venus
    398600.440,       # Earth
    42828.3,          # Mars
    126686511,        # Jupiter
    37931207.8,       # Saturn
    5793966,          # Uranus
    6835107,          # Neptune
    872.4             # Pluto
])

total_mu = np.sum(mu)
n_bodies = 10

# ------------------------------------------------------------
# Solar-System barycentre relative to the Sun
# ------------------------------------------------------------

x = np.zeros(n_use)
y = np.zeros(n_use)
z = np.zeros(n_use)

for i in range(n_bodies):

    weight = mu[i] / total_mu

    x += weight * nc.variables['pos'][:n_use, i, 0]
    y += weight * nc.variables['pos'][:n_use, i, 1]
    z += weight * nc.variables['pos'][:n_use, i, 2]

# Convert to AU
AU = 1.495978707e11

x = x / AU
y = y / AU
z = z / AU

# ------------------------------------------------------------
# Analyse the barycentric motion in the ecliptic plane
# ------------------------------------------------------------

# x + iy preserves the rotating vector rather than taking
# its scalar magnitude before the Fourier transform.
signal = x + 1j * y

# Remove any constant offset
signal = signal - np.mean(signal)

# Hann window to reduce spectral leakage
window = np.hanning(n_use)
signal = signal * window

# Full complex FFT
fft = np.fft.fft(signal, n=NFFT)

freq = np.fft.fftfreq(NFFT, d=dt)

# Keep positive frequencies only
use = freq > 0.0

freq = freq[use]
fft = fft[use]

period = 1.0 / freq

# Fourier amplitude in AU
#
# Do not multiply by 2 here: x + iy is already a complex
# representation of the rotating displacement.
amplitude = np.abs(fft) / np.sum(window)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(9, 6))

ax.plot(
    period,
    amplitude,
    'k-',
    linewidth=0.8
)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(0.1, 1000.0)

ax.set_xlabel(r'time period (years)')
ax.set_ylabel(r'Fourier amplitude (AU)')
ax.set_title(r"Sun's barycentric motion")

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

fig.subplots_adjust(
    left=0.13,
    right=0.97,
    bottom=0.13,
    top=0.92
)

nc.close()

outfile = '/tmp/' + username + '/fourier_sun.png'

fig.savefig(
    outfile,
    dpi=300
)

plt.close(fig)

print('Saved:', outfile)
