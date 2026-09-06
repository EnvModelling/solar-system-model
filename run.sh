#!/usr/bin/env bash
#set -euo pipefail

output_dir="/tmp/${USER}"
tmp_namelist="$(mktemp "${TMPDIR:-/tmp}/solar_system_namelist.XXXXXX")"

cleanup() {
    rm -f "${tmp_namelist}"
}
trap cleanup EXIT

mkdir -p "${output_dir}"

sed "s|/tmp/output.nc|${output_dir}/output.nc|g" namelist.in > "${tmp_namelist}"

./main.exe "${tmp_namelist}"
