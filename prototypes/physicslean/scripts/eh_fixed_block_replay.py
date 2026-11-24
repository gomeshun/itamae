"""Replay a continuous top-hat raw integral on unchanged cells of consecutive chosen octaves.

The generated proof recomputes each cell's EH spectrum, its integral and its
window from exact edges, before joining consecutive blocks. No table input.
"""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess
import tempfile
import time

from eh_plot_variance import PROJECT
from eh_sweep_certificate import add, box
from table_replay import literal
from table_window import check_modules


def declarations(octave, last_octave, per_octave, block_size, radius):
    edges = [Q(2**octave, 10**8)]
    for _ in range(octave, last_octave + 1):
        start = edges[-1]
        stop = min(2 * start, Q(10**10))
        if start >= stop:
            raise ValueError("octave is outside the matched band")
        edges.extend(start + (stop - start) * j / per_octave for j in range(1, per_octave + 1))
    groups = [
        edges[i : min(i + block_size, len(edges) - 1) + 1]
        for i in range(0, len(edges) - 1, block_size)
    ]
    src = [
        "import SashimiEhFixedReplay",
        "open FloatLib.Numerics SashimiEhFixedReplay SashimiEhFixedCells SashimiEisensteinIntegral SashimiEhFullGrid",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 60000000",
        "set_option Elab.async false",
        f"private def R : RationalInterval := .point {literal(radius)}",
    ]
    return edges, groups, src


def replay(
    octave, last_octave, per_octave, block_size, radius, bits, degree, output, timeout, jobs
):
    edges, groups, src = declarations(octave, last_octave, per_octave, block_size, radius)
    blocks = [f"({literal(g[0])}, [{', '.join(literal(k) for k in g[1:])}])" for g in groups]
    native = [
        *src,
        f"private def blocks : List (ℚ × List ℚ) := [{', '.join(blocks)}]",
        "def main : IO Unit := do",
        f'  let some C := boxes? matched {bits} {degree} | throw (IO.userError "coefficients unresolved")',
        '  for I in C do IO.println s!"c {I.lo} {I.hi}"',
        "  for (a, ends) in blocks do",
        f'    let some rows := build? matched C {bits} {degree} a ends | throw (IO.userError "block unresolved")',
        '    for x in rows do IO.println s!"cell {x.right} {x.pointwise.lo} {x.pointwise.hi} {x.integral.lo} {x.integral.hi}"',
        f"    let I := variance R {bits} {degree} a rows",
        '    IO.println s!"sum {I.lo} {I.hi}"',
    ]
    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "Evaluate.lean"
        path.write_text("\n".join(native) + "\n")
        r = subprocess.run(
            ["lake", "env", "lean", "--run", str(path)],
            cwd=PROJECT,
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout,
        )
    wall = time.perf_counter() - start
    coefficients, cells, values = [], [], []
    for line in r.stdout.splitlines():
        tag, *raw = line.split()
        x = tuple(map(Q, raw))
        if tag == "c" and len(x) == 2:
            coefficients.append(x)
        elif tag == "cell" and len(x) == 5:
            cells.append(x)
        elif tag == "sum" and len(x) == 2:
            values.append(x)
        else:
            raise ValueError("unexpected numerical output")
    if len(values) != len(groups) or not coefficients or [c[0] for c in cells] != edges[1:]:
        raise ValueError("incomplete runtime blocks or coefficients")
    nodes = []
    rr = literal(radius)
    src.extend(
        [
            f"private def exactFixedCoefficients : List RationalInterval := [{', '.join(box(v) for v in coefficients)}]",
            f"private theorem hCoefficients : boxes? matched {bits} {degree} = some exactFixedCoefficients := by decide +kernel",
            f"private theorem hR : R.Contains ({rr}:ℝ) := by norm_num [R,RationalInterval.point,RationalInterval.Contains]",
        ]
    )
    header = [line.removeprefix("private ") for line in src]
    leaves = []
    j = 0
    for i, (g, v) in enumerate(zip(groups, values)):
        leaf = []
        names = []
        valid = []
        for lo, hi in zip(g, g[1:]):
            right, pl, pu, al, au = cells[j]
            name = f"item{j}"
            leaf.extend(
                [
                    f"private def {name} : Cell := ⟨{literal(right)}, {box((pl, pu))}, {box((al, au))}⟩",
                    f"private theorem evaluated{j} : buildCell? matched exactFixedCoefficients {literal(lo)} {literal(hi)} {bits} {degree} = some {name} := by decide +kernel",
                    f"private theorem valid{j} : Valid (SashimiEhSweep.weightReal matched) {literal(lo)} {name} := (buildCell_sound hCoefficients evaluated{j}).2",
                ]
            )
            names.append(name)
            valid.append(f"valid{j}")
            j += 1
        lo, hi = literal(g[0]), literal(g[-1])
        leaf.extend(
            [
                f"private def cells{i} : List Cell := [{', '.join(names)}]",
                f"private theorem good{i} : Good (SashimiEhSweep.weightReal matched) {lo} cells{i} := ⟨{', '.join(valid)}, trivial⟩",
                f"private theorem total{i} : variance R {bits} {degree} {lo} cells{i} = {box(v)} := by decide +kernel",
                f"private theorem part{i} : {box(v)}.Contains (raw matched {lo} {hi} {rr}) := by",
                f"  have h := variance_sound (R:=R) (r:=({rr}:ℝ)) (bits:={bits}) (degree:={degree}) good{i} hR (by decide +kernel)",
                f"  have he : lastEdge {lo} cells{i} = {hi} := by decide +kernel",
                f"  rwa [total{i}, he] at h",
            ]
        )
        leaves.append(leaf)
        src.extend(leaf)
        nodes.append((i, g[0], g[-1], v))
    assembly_start = len(src)
    i = len(groups)
    while len(nodes) > 1:
        joined = []
        for pair_start in range(0, len(nodes), 2):
            if pair_start + 1 == len(nodes):
                joined.append(nodes[pair_start])
                continue
            left, lo, mid, x = nodes[pair_start]
            right, mid2, hi, y = nodes[pair_start + 1]
            if mid != mid2:
                raise ValueError("disconnected blocks")
            v = add(x, y)
            src.extend(
                [
                    f"private theorem part{i} : {box(v)}.Contains (raw matched {literal(lo)} {literal(hi)} {rr}) := by",
                    f"  have h := @join ({rr}:ℝ) {literal(lo)} {literal(mid)} {literal(hi)} {box(x)} {box(y)} (by norm_num) (by decide +kernel) (by decide +kernel) (by decide +kernel) part{left} part{right}",
                    f"  have he : {box(x)}.add {box(y)} = {box(v)} := by decide +kernel",
                    "  rwa [he] at h",
                ]
            )
            joined.append((i, lo, hi, v))
            i += 1
        nodes = joined
    j, lo, hi, v = nodes[0]
    src.extend(
        [
            f"theorem certifiedFixedBand : {box(v)}.Contains (raw matched {literal(lo)} {literal(hi)} {rr}) := part{j}",
            "#print axioms certifiedFixedBand",
        ]
    )
    output.write_text("\n".join(src) + "\n")
    start = time.perf_counter()
    if jobs == 1:
        checked = subprocess.run(
            ["lake", "env", "lean", "-j1", "-M5100", str(output.resolve())],
            cwd=PROJECT,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
        expected = "'certifiedFixedBand' depends on axioms: [propext, Classical.choice, Quot.sound]"
        if checked.returncode or expected not in checked.stdout.splitlines():
            raise RuntimeError(checked.stdout + checked.stderr)
    else:
        checked_statement = f"theorem concrete : {box(v)}.Contains (raw matched {literal(lo)} {literal(hi)} {rr}) := certifiedFixedBand"
        assembly = [*src[assembly_start:], checked_statement, "#print axioms concrete"]
        check_modules(header, leaves, assembly, 5100, timeout, jobs)
    return {
        "kernel_replayed": True,
        "theorem": "certifiedFixedBand",
        "kernel_axioms": ["propext", "Classical.choice", "Quot.sound"],
        "scope": "Unnormalized top-hat variance integral on the stated finite k band; the universal constant 1/(2π²) is omitted.",
        "octave": octave,
        "last_octave": last_octave,
        "cells_per_octave": per_octave,
        "k_band_h_per_Mpc": [str(edges[0]), str(edges[-1])],
        "cells": len(edges) - 1,
        "block_size": block_size,
        "parallel_jobs": jobs,
        "radius_mpc_h": str(radius),
        "bits": bits,
        "degree": degree,
        "raw": {"lo": str(v[0]), "hi": str(v[1])},
        "runtime_seconds": wall,
        "kernel_seconds": time.perf_counter() - start,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("proof", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--octave", type=int, default=23)
    p.add_argument(
        "--last-octave", type=int, help="inclusive last original octave; defaults to --octave"
    )
    p.add_argument("--cells-per-octave", type=int, default=128)
    p.add_argument("--block-size", type=int, default=4)
    p.add_argument("--radius", type=Q, default=Q(8))
    p.add_argument("--bits", type=int, default=112)
    p.add_argument("--degree", type=int, default=16)
    p.add_argument(
        "--jobs", type=int, default=1, help="parallel independent groups of up to eight blocks"
    )
    p.add_argument("--timeout", type=int, default=3600)
    a = p.parse_args()
    last_octave = a.octave if a.last_octave is None else a.last_octave
    if not (
        0 <= a.octave <= last_octave < 60
        and min(a.cells_per_octave, a.block_size, a.radius, a.bits, a.degree, a.timeout, a.jobs) > 0
    ):
        p.error("positive numerical settings and an octave in [0,60) are required")
    data = replay(
        a.octave,
        last_octave,
        a.cells_per_octave,
        a.block_size,
        a.radius,
        a.bits,
        a.degree,
        a.proof,
        a.timeout,
        a.jobs,
    )
    a.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
