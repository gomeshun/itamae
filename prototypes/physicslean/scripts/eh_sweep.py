"""Enclose the continuous finite-band EH variance for many radii on a shared mesh."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import time

PROJECT = Path(__file__).resolve().parents[1]


def evaluate(grid, numerical, *, subdivisions, bits, degree, transport=None):
    if min(subdivisions, bits, degree) <= 0:
        raise ValueError("positive subdivisions, precision, and degree required")
    masses = [Fraction(x) for x in grid["masses_msun"]]
    host = Fraction(numerical["host_progenitor"]["mass_msun"])
    input_text = "".join(f"{m}\n" for m in [*masses, host])
    start = time.monotonic()
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

    def box(lo, hi):
        if Fraction(lo) > Fraction(hi):
            raise RuntimeError("inverted rational enclosure")
        return {"lo": lo, "hi": hi}

    bins = []
    rows = []
    anchor = None
    for line in result.stdout.splitlines():
        match line.split():
            case ["bin", lo, hi, p, q]:
                bins.append({"left": lo, "right": hi, "weight": box(p, q)})
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
                    "bins": [[x["right"], x["weight"]["lo"], x["weight"]["hi"]] for x in bins],
                }
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
    parser.add_argument("--transport", type=Path)
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
    )
    args.output.write_text(json.dumps(r, indent=2) + "\n")


if __name__ == "__main__":
    main()
