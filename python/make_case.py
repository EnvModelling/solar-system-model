#!/usr/bin/env python3
"""Create an SSM case namelist from the reference namelist.

Examples
--------
Sun-only / Keplerian:
    python3 scripts/make_case.py namelist.in runs/keplerian/namelist.in \
        --mode keplerian --output runs/keplerian/output.nc

All Newtonian interactions:
    python3 scripts/make_case.py namelist.in runs/interacting/namelist.in \
        --mode interacting --output runs/interacting/output.nc
"""

from argparse import ArgumentParser
from pathlib import Path
import re


def bool_fortran(value):
    return ".true." if value else ".false."


def replace_scalar(lines, name, value):
    pattern = re.compile(rf"^(\s*){re.escape(name)}\s*=", re.IGNORECASE)
    out = []
    found = False

    for line in lines:
        if not line.lstrip().startswith("!") and pattern.match(line):
            indent = pattern.match(line).group(1)
            comment = ""
            if "!" in line:
                comment = " !" + line.split("!", 1)[1].rstrip("\n")
            out.append(f"{indent}{name}={value},{comment}\n")
            found = True
        else:
            out.append(line)

    if not found:
        raise RuntimeError(f"Could not find active namelist assignment for {name}")
    return out


def replace_interact(lines, mode):
    if mode == "keplerian":
        vals = [True] + [False] * 9
    elif mode == "interacting":
        vals = [True] * 10
    else:
        raise ValueError(mode)

    first = "  interact(1:10)=" + ",".join(bool_fortran(v) for v in vals[:5]) + ", \n"
    second = "                 " + ",".join(bool_fortran(v) for v in vals[5:]) + "/ \n"

    out = []
    i = 0
    found = False

    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()

        if (not stripped.startswith("!")) and re.match(r"interact\s*\(\s*1\s*:\s*10\s*\)\s*=", stripped, re.I):
            out.extend([first, second])
            found = True
            i += 1

            # The existing active interact assignment ends at the namelist '/'.
            while i < len(lines):
                if "/" in lines[i]:
                    i += 1
                    break
                i += 1
            continue

        out.append(line)
        i += 1

    if not found:
        raise RuntimeError("Could not find active interact(1:10) assignment")
    return out


def main():
    parser = ArgumentParser()
    parser.add_argument("template")
    parser.add_argument("destination")
    parser.add_argument("--mode", choices=("keplerian", "interacting"), required=True)
    parser.add_argument("--output", required=True, help="NetCDF output filename")
    parser.add_argument("--dt", type=float, help="Output interval / DVODE HMAX in days")
    parser.add_argument("--tfinal", type=float, help="Run length in years")
    parser.add_argument("--interval-io", type=float, help="Progress print interval in years")
    parser.add_argument("--general-relativity", choices=("true", "false"))
    args = parser.parse_args()

    lines = Path(args.template).read_text().splitlines(keepends=True)

    lines = replace_scalar(lines, "outputfile01", repr(args.output))

    if args.dt is not None:
        lines = replace_scalar(lines, "dt", f"{args.dt:.16g}d0")
    if args.tfinal is not None:
        lines = replace_scalar(lines, "tfinal", f"{args.tfinal:.16g}d0")
    if args.interval_io is not None:
        lines = replace_scalar(lines, "interval_io", f"{args.interval_io:.16g}d0")
    if args.general_relativity is not None:
        lines = replace_scalar(
            lines,
            "general_relativity",
            bool_fortran(args.general_relativity == "true"),
        )

    lines = replace_interact(lines, args.mode)

    destination = Path(args.destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(lines))
    print(destination)


if __name__ == "__main__":
    main()
