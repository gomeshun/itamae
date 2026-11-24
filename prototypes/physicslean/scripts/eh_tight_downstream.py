"""Propagate fixed-cell EH enclosures to the full spectrum, EPS and host at the published grid."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from time import monotonic

SCRIPTS = Path(__file__).resolve().parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(script, *args):
    start = monotonic()
    subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], check=True)
    print(f"{script} completed in {monotonic() - start:.3f} s", flush=True)


def record_inputs(path, mesh, **inputs):
    report = json.loads(path.read_bytes())
    report["mesh"] = mesh
    report["inputs_sha256"] = {k: digest(v) for k, v in inputs.items()}
    path.write_text(json.dumps(report, indent=2) + "\n")


def verify_mesh(candidate, reference):
    control = json.loads(reference.read_bytes())

    def edges(doc):
        return [(Fraction(c["left"]), Fraction(c["right"])) for c in doc["shared_bins"]]

    def radii(doc):
        return [
            (Fraction(r["mass_msun"]), Fraction(r["radius"]["lo"]), Fraction(r["radius"]["hi"]))
            for r in [*doc["rows"], doc["host_progenitor"]]
        ]

    grid, masses = edges(control), radii(control)
    if (
        grid != edges(candidate)
        or masses != radii(candidate)
        or candidate["anchor"]["radius"] != control["anchor"]["radius"]
    ):
        raise ValueError(
            "Fixed-cell and reference sweeps differ in their physical integration grid"
        )
    for k in ("bits", "degree", "cells", "subdivisions_per_octave"):
        if candidate["mesh"][k] != control["mesh"][k]:
            raise ValueError(f"Fixed-cell and reference sweeps differ in {k}")
    if len(grid) != candidate["mesh"]["cells"]:
        raise ValueError("Reported cell count disagrees with the integration grid")
    return {
        "reference_sweep_sha256": digest(reference),
        "matched_cells": len(grid),
        "matched_radii_including_progenitor": len(masses),
    }


def propagate(sweep, growth, out, reference_sweep, finite_slope=None):
    candidate = json.loads(sweep.read_bytes())
    mesh = candidate["mesh"]
    if mesh.get("enclosure_method") != "fixed_cell_enclosure":
        raise ValueError("Expected the fixed-cell sweep, with its original mesh metadata")
    same_grid = verify_mesh(candidate, reference_sweep)
    out.mkdir(parents=True, exist_ok=True)
    full = out / "full-tight.json"
    run("eh_full_grid.py", sweep, full)
    report = json.loads(full.read_bytes())
    report["same_grid"] = same_grid
    full.write_text(json.dumps(report, indent=2) + "\n")
    run("eps_eh_kernel_bounds.py", full, growth, out / "kernel-bounds-tight.json")
    host = out / "host-bounds-tight.json"
    run("eh_host_bounds.py", sweep, growth, host)
    record_inputs(host, mesh, sweep=sweep, growth=growth)
    if finite_slope:
        slopes = out / "derivative-full-tight.json"
        run(
            "eh_full_derivative.py",
            sweep,
            slopes,
            "--finite-evaluation",
            finite_slope,
            "--certificates",
            out / "proofs" / "full-slopes-tight",
        )
        record_inputs(slopes, mesh, sweep=sweep, finite_derivative=finite_slope)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("sweep", type=Path)
    cli.add_argument("growth", type=Path)
    cli.add_argument("output", type=Path)
    cli.add_argument("--reference-sweep", type=Path, required=True)
    cli.add_argument("--finite-slope", type=Path)
    args = cli.parse_args()
    propagate(args.sweep, args.growth, args.output, args.reference_sweep, args.finite_slope)
