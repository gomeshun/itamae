"""Evaluate the signed window derivative using the analytic EH sweep's shared power mesh."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import time


def evaluate(mesh_file, grid_file, tolerance):
    mesh_data = mesh_file.read_bytes()
    mesh = json.loads(mesh_data)
    grid = json.loads(grid_file.read_text())
    masses = [Fraction(x) for x in grid["masses_msun"]]
    cells = mesh["bins"]
    if not cells or Fraction(mesh["a"]) <= 0 or Fraction(mesh["b"]) != Fraction(cells[-1][0]):
        raise ValueError("invalid mesh bounds")
    source = [f"cell {r} {lo} {hi}" for r, lo, hi in cells]
    source.extend(f"mass {m}" for m in masses)
    cmd = [
        "lake",
        "env",
        "lean",
        "--run",
        "scripts/EhDerivativeBounds.lean",
        mesh["a"],
        str(mesh["bits"]),
        str(mesh["degree"]),
        str(tolerance),
    ]
    started = time.monotonic()
    run = subprocess.run(
        cmd,
        input="\n".join(source) + "\n",
        text=True,
        capture_output=True,
        cwd=Path(__file__).resolve().parents[1],
    )
    if run.returncode:
        raise RuntimeError((run.stdout + run.stderr)[-2000:])
    lines = [line.split() for line in run.stdout.splitlines()]
    if not lines or lines[0][0] != "anchor" or len(lines) != 1 + len(masses):
        raise RuntimeError("incomplete derivative response")
    rows = []
    for m, label, line in zip(masses, grid["masses_msun"], lines[1:]):
        status, value, *bounds = line
        if Fraction(value) != m or status not in ("enclosed", "unresolved"):
            raise RuntimeError("unexpected derivative response")
        row = {"mass_msun": label, "status": status, "dS_dM": None}
        if status == "enclosed":
            lo, hi = map(Fraction, bounds)
            if lo > hi or lo <= 0 <= hi or hi - lo > tolerance:
                raise RuntimeError("invalid derivative enclosure")
            row["dS_dM"] = {"lo": str(lo), "hi": str(hi)}
        rows.append(row)
    return {
        "target": "finite-band dS/dM per physical solar mass",
        "mesh_sha256": hashlib.sha256(mesh_data).hexdigest(),
        "k_band": [mesh["a"], mesh["b"]],
        "cells": len(cells),
        "bits": mesh["bits"],
        "degree": mesh["degree"],
        "absolute_width_budget": str(tolerance),
        "anchor": {"lo": lines[0][1], "hi": lines[0][2]},
        "runtime_seconds": time.monotonic() - started,
        "verification": "Kernel-proved derivative_sound requires valid source power intervals. This pass reuses the recorded mesh; it does not replay each supplied cell in the kernel.",
        "rows": rows,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mesh", type=Path)
    p.add_argument("grid", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--tolerance", type=Fraction, default=Fraction(10**20))
    args = p.parse_args()
    args.output.write_text(
        json.dumps(evaluate(args.mesh, args.grid, args.tolerance), indent=2) + "\n"
    )
