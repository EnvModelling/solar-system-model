#!/usr/bin/env bash
set -euo pipefail

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${root}"

PYTHON="${PYTHON:-python3}"
TFINAL_YEARS="${TFINAL_YEARS:-200}"
DT_DAYS="${DT_DAYS:-5}"
FRAMES="${FRAMES:-300}"

mkdir -p runs/keplerian runs/interacting runs/barycentric docs/figures

echo "Creating case namelists..."
"${PYTHON}" scripts/make_case.py namelist.in runs/keplerian/namelist.in \
    --mode keplerian \
    --output runs/keplerian/output.nc \
    --dt "${DT_DAYS}" \
    --tfinal "${TFINAL_YEARS}" \
    --interval-io 10 \
    --general-relativity false

"${PYTHON}" scripts/make_case.py namelist.in runs/interacting/namelist.in \
    --mode interacting \
    --output runs/interacting/output.nc \
    --dt "${DT_DAYS}" \
    --tfinal "${TFINAL_YEARS}" \
    --interval-io 10 \
    --general-relativity false

"${PYTHON}" scripts/make_case.py namelist.in runs/barycentric/namelist.in \
    --mode interacting \
    --output runs/barycentric/output.nc \
    --dt "${DT_DAYS}" \
    --tfinal "${TFINAL_YEARS}" \
    --interval-io 10 \
    --general-relativity false

echo
echo "Building the normal heliocentric executable..."
make clean
make HELIOCENTRIC=1

echo
echo "Running Sun-only / Keplerian case..."
./main.exe runs/keplerian/namelist.in

echo
echo "Running all-interaction heliocentric case..."
./main.exe runs/interacting/namelist.in

echo
echo "Making static figures..."
"${PYTHON}" scripts/plot_ssm.py runs/keplerian/output.nc \
    --outdir docs/figures --prefix keplerian --origin sun
"${PYTHON}" scripts/plot_ssm.py runs/interacting/output.nc \
    --outdir docs/figures --prefix interacting --origin sun
"${PYTHON}" scripts/compare_runs.py \
    runs/keplerian/output.nc \
    runs/interacting/output.nc \
    docs/figures/interaction_difference.png

echo
echo "Making heliocentric GIFs..."
"${PYTHON}" scripts/animate_ssm.py \
    runs/keplerian/output.nc docs/figures/keplerian_all.gif \
    --system all --origin sun --frames "${FRAMES}"

"${PYTHON}" scripts/animate_ssm.py \
    runs/interacting/output.nc docs/figures/interacting_all.gif \
    --system all --origin sun --frames "${FRAMES}"

# The first five years make a clearer inner-planet animation.
"${PYTHON}" scripts/animate_ssm.py \
    runs/interacting/output.nc docs/figures/interacting_inner.gif \
    --system inner --origin sun --max-years 5 --frames "${FRAMES}"

echo
echo "Building without -Dheliocentric so the Sun is allowed to wobble..."
make clean
make HELIOCENTRIC=0

echo
echo "Running all-interaction non-heliocentric case..."
./main.exe runs/barycentric/namelist.in

echo
echo "Plotting/animating in the instantaneous barycentric frame..."
"${PYTHON}" scripts/plot_ssm.py runs/barycentric/output.nc \
    --outdir docs/figures --prefix barycentric --origin barycentre

"${PYTHON}" scripts/animate_ssm.py \
    runs/barycentric/output.nc docs/figures/barycentric_wobble.gif \
    --system all --origin barycentre --frames "${FRAMES}"

echo
echo "Restoring the historical/default heliocentric build..."
make clean
make HELIOCENTRIC=1

echo
echo "Done. README figures are in:"
echo "  ${root}/docs/figures"
