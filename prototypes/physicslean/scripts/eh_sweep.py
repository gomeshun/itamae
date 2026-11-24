"""Enclose the continuous finite-band EH variance for many radii on a shared mesh."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import time

if __package__:
    from .eh_fixed_workers import mass_totals, power_cells
else:
    from eh_fixed_workers import mass_totals, power_cells

PROJECT = Path(__file__).resolve().parents[1]


def evaluate(
    grid, numerical, *, subdivisions, bits, degree, transport=None, derivative_output=None, jobs=1
):
    fixed = derivative_output is not None
    if fixed and transport is None:
        raise ValueError("fixed-cell evaluation requires the power and integral transport")

    if min(subdivisions, bits, degree, jobs) <= 0:
        raise ValueError("positive subdivisions, precision, degree and worker count required")
    masses = [Fraction(x) for x in grid["masses_msun"]]
    host = Fraction(numerical["host_progenitor"]["mass_msun"])
    input_text = "".join(f"{m}\n" for m in [*masses, host])
    start = time.monotonic()
    if fixed:
        cells = power_cells(subdivisions, bits, degree, jobs)
        output = mass_totals(cells, [*masses, host], bits, degree, jobs)
    else:
        result = subprocess.run(
            [
                "lake",
                "env",
                "lean",
                "--run",
                "scripts/EhSweep.lean",
                str(subdivisions),
                str(bits),
                str(degree),
            ],
            cwd=PROJECT,
            input=input_text,
            text=True,
            capture_output=True,
        )
        if result.returncode:
            raise RuntimeError(result.stderr + result.stdout[-1200:])
        output = result.stdout

    def box(lo, hi):
        if Fraction(lo) > Fraction(hi):
            raise RuntimeError("inverted rational enclosure")
        return {"lo": lo, "hi": hi}

    bins = []
    rows = []
    anchor = None
    slopes = []
    for line in output.splitlines():
        match line.split():
            case ["bin", lo, hi, p, q]:
                bins.append({"left": lo, "right": hi, "weight": box(p, q)})
            case ["cell", lo, hi, p, q, al, au] if fixed:
                bins.append({"left": lo, "right": hi, "weight": box(p, q), "integral": box(al, au)})
            case ["slope", mass, lo, hi] if fixed:
                slope = box(lo, hi)
                resolved = Fraction(hi) < 0 or Fraction(lo) > 0
                slopes.append(
                    {
                        "mass_msun": mass,
                        "status": "enclosed" if resolved else "unresolved",
                        "dS_dM": slope if resolved else None,
                    }
                )
            case ["anchor", lo, hi]:
                anchor = {"radius": box("8", "8"), "raw": box(lo, hi)}
            case ["mass", mass, r, s, l, u, low, high]:
                rows.append(
                    {"mass_msun": mass, "radius": box(r, s), "raw": box(l, u), "S": box(low, high)}
                )
            case _:
                raise RuntimeError(f"unexpected Lean output {line[:120]}")
    if (
        [Fraction(r["mass_msun"]) for r in rows] != [*masses, host]
        or not bins
        or anchor is None
        or Fraction(bins[0]["left"]) != Fraction(1, 10**8)
        or Fraction(bins[-1]["right"]) != 10**10
        or any(left["right"] != right["left"] for left, right in zip(bins, bins[1:]))
    ):
        raise RuntimeError("incomplete or unordered mass or wavenumber grid")
    if transport is not None:
        transport.write_text(
            json.dumps(
                {
                    "a": bins[0]["left"],
                    "b": bins[-1]["right"],
                    "bits": bits,
                    "degree": degree,
                    "coefficient_bits": bits,
                    "coefficient_degree": degree,
                    **({"integrals": [x["integral"] for x in bins]} if fixed else {}),
                    "bins": [[x["right"], x["weight"]["lo"], x["weight"]["hi"]] for x in bins],
                }
            )
            + "\n"
        )
    if fixed:
        if [Fraction(s["mass_msun"]) for s in slopes] != [*masses, host]:
            raise RuntimeError("incomplete signed derivative sweep")
        derivative_output.write_text(
            json.dumps(
                {
                    "target": "finite-band dS/dM per physical solar mass",
                    "mesh_sha256": hashlib.sha256(transport.read_bytes()).hexdigest(),
                    "k_band": [bins[0]["left"], bins[-1]["right"]],
                    "cells": len(bins),
                    "bits": bits,
                    "degree": degree,
                    "absolute_width_budget": str(10**20),
                    "anchor": anchor["raw"],
                    "runtime_seconds": time.monotonic() - start,
                    "verification": "Fixed-cell slope_sound; shared finite inputs evaluated by Lean runtime.",
                    "rows": slopes[:-1],
                },
                indent=2,
            )
            + "\n"
        )
    native = [*numerical["mass_rows"], numerical["host_progenitor"]]
    outside = [
        r["mass_msun"]
        for r, n in zip(rows, native)
        if not Fraction(r["S"]["lo"]) <= Fraction.from_float(n["S"]) <= Fraction(r["S"]["hi"])
    ]
    return {
        "status": "proved evaluator; concrete runtime (kernel replay tracked separately)",
        "target": "continuous, sigma8-normalized finite-band EH variance; raw excludes 1/(2*pi^2)",
        "k_band": [bins[0]["left"], bins[-1]["right"]],
        "mesh": {
            "bits": bits,
            "degree": degree,
            "subdivisions_per_octave": subdivisions,
            "cells": len(bins),
            **({"enclosure_method": "fixed_cell_enclosure"} if fixed else {}),
        },
        "runtime_seconds": time.monotonic() - start,
        "shared_bins": bins,
        "anchor": anchor,
        "rows": rows[:-1],
        "host_progenitor": rows[-1],
        "native_quadrature_outliers": outside,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("grid", "numerical", "output"):
        parser.add_argument("--" + arg, required=True, type=Path)
    parser.add_argument(
        "--jobs", type=int, default=2, help="number of independent fixed-cell workers"
    )
    parser.add_argument("--transport", type=Path)
    parser.add_argument(
        "--fixed-cell-slopes",
        type=Path,
        help="use tighter enclosures on the identical cells and write signed derivatives",
    )
    parser.add_argument("--subdivisions", type=int, default=16)
    parser.add_argument("--bits", type=int, default=112)
    parser.add_argument("--degree", type=int, default=16)
    args = parser.parse_args()
    r = evaluate(
        json.loads(args.grid.read_text()),
        json.loads(args.numerical.read_text()),
        subdivisions=args.subdivisions,
        bits=args.bits,
        degree=args.degree,
        transport=args.transport,
        derivative_output=args.fixed_cell_slopes,
        jobs=args.jobs,
    )
    args.output.write_text(json.dumps(r, indent=2) + "\n")


if __name__ == "__main__":
    main()
