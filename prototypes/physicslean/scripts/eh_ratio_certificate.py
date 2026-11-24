"""Kernel-replay shared-weight normalized endpoints on a certified EH mesh.

Create the upstream WindowCertificate modules with eh_sweep_certificate.py.
The output Lean file recomputes both signed sums and proves the full positive-k
physical S for one mass. Python-proposed endpoints and subtotals are untrusted.
"""

import argparse
from fractions import Fraction as Q
from pathlib import Path
import json
import os
import subprocess
import tempfile
import time

from eh_shared_ratio import PROJECT, read_cells, residual, write_request
from eh_sweep_certificate import add, box, pair
from eh_sweep_masses import _chunks
from table_replay import literal


def source(sweep, row, qbits, modules, parts):
    bits, degree = (sweep["mesh"][k] for k in ("bits", "degree"))
    i, mass = row["index"], Q(row["mass_msun"])
    R, A = box(pair(row["radius"])), box(pair(sweep["anchor"]["raw"]))
    qs = pair(row["ratio"])
    starts, joins = _chunks(modules)
    edges = [Q(x["left"]) for x in sweep["shared_bins"]]
    position = {k: j for j, k in enumerate(edges)}
    ordered = sorted(starts)
    if starts[ordered[0]] != edges[0]:
        raise ValueError("certified chunk starts differ from the cell grid")
    src = [
        "import WindowCertificate",
        "import CertifiedIntegralRatio",
        "open CertifiedIntegralRatio SashimiEhSweep SashimiEhFullGrid SashimiEisensteinIntegral FloatLib.Numerics",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 60000000",
        "set_option Elab.async false",
        f"private def R : RationalInterval := {R}",
        "private theorem low_eq : lowEdge = (1/100000000:ℚ) := by decide +kernel",
        "private theorem high_eq : highEdge = (10000000000/1:ℚ) := by decide +kernel",
    ]
    for name, q in zip(("L", "U"), qs):
        args = f"residual R {bits} {degree} {literal(q)}"
        sums = {}
        for j, ix in enumerate(ordered):
            lower = position[starts[ix]]
            upper = position[starts[ordered[j + 1]]] if j + 1 < len(ordered) else len(edges)
            if (q * 2**qbits).denominator != 1:
                raise ValueError("expected dyadic endpoint on the selected ratio grid")
            sums[ix] = residual(parts[lower:upper], int(q * 2**qbits), bits, qbits)
            src.append(
                f"private theorem {name}{ix} : {args} {literal(starts[ix])} ehBins{ix} = {box(sums[ix])} := by decide +kernel"
            )
        for ix, left, right in joins:
            starts[ix] = starts[left]
            sums[ix] = add(sums[left], sums[right])
            src.append(
                f"private theorem {name}{ix} : {args} {literal(starts[ix])} ehBins{ix} = {box(sums[ix])} := by\n"
                f"  rw [ehBins{ix}, residual_append, end{left}, {name}{left}, {name}{right}]; decide +kernel"
            )
        ix = max(sums)
        declared = pair(row["residual_at_lower" if name == "L" else "residual_at_upper"])
        if declared != sums[ix]:
            raise ValueError("declared signed sum differs from the recomputed cell sum")
        src.extend(
            [
                f"private theorem sign{name} {{x:ℝ}} (hx : R.Contains x) : {box(declared)}.Contains "
                f"(raw matched lowEdge highEdge x - {literal(q)} * raw matched lowEdge highEdge 8) := by",
                f"  have h := residual_sound (bits := {bits}) (degree := {degree}) (q := {literal(q)}) "
                f"(by rw [low_eq]; exact certifiedSpectrum) (by rw [low_eq,high_eq]; exact certifiedBandEnd) hx (by decide +kernel) (by decide +kernel)",
                f"  rw [low_eq, {name}{ix}] at h",
                "  exact h",
            ]
        )
    src.extend(
        [
            f"theorem sharedNormalized{i} : {box(pair(row['S']))}.Contains (S {literal(mass)}) := by",
            f"  have hr := radiusBox_sound (by decide +kernel : radiusBox? matched {literal(mass)} {bits} {degree} = some R)",
            f"  have ha : {A}.Contains (raw matched lowEdge highEdge 8) := by rw [low_eq, high_eq]; exact concrete.2.1",
            "  have h := full_sound hr.1 (by decide +kernel) (by decide +kernel) ha (by decide +kernel) "
            "(signL hr) (signU hr) (by decide +kernel)",
            f"  have he : ((⟨{literal(qs[0])}, {literal(qs[1])}⟩ : RationalInterval).scaleNonnegative (matched.sigma8^2)) = {box(pair(row['S']))} := by decide +kernel",
            "  rwa [he] at h",
            f"#print axioms sharedNormalized{i}",
        ]
    )
    return "\n".join(src) + "\n"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("sweep", "result", "modules", "output"):
        p.add_argument(name, type=Path)
    p.add_argument("--index", type=int, required=True)
    a = p.parse_args()
    sweep = json.loads(a.sweep.read_text())
    if sweep["mesh"].get("enclosure_method") == "fixed_cell_enclosure":
        p.error("This witness builder expects original rectangular spectral bins")
    report = json.loads(a.result.read_text())
    (row,) = [r for r in report["rows"] if r["index"] == a.index]
    ids = [a.index]
    sr = sweep["rows"] + [sweep["host_progenitor"]]
    if row["mass_msun"] != sr[a.index]["mass_msun"] or row["radius"] != sr[a.index]["radius"]:
        raise ValueError("mismatched input and ratio row")
    with tempfile.TemporaryDirectory() as temp:
        request, cells = Path(temp) / "input", Path(temp) / "cells"
        write_request(sweep, ids, request)
        with request.open() as f, cells.open("w") as out:
            subprocess.run(
                [
                    "lake",
                    "env",
                    "lean",
                    "--run",
                    "scripts/EhRatioCells.lean",
                    str(sweep["mesh"]["bits"]),
                    str(sweep["mesh"]["degree"]),
                ],
                stdin=f,
                stdout=out,
                cwd=PROJECT,
                check=True,
            )
        ((_, e, w, d),) = read_cells(cells, ids, len(sweep["shared_bins"]))
        a.output.write_text(source(sweep, row, report["qbits"], a.modules, list(zip(e, w, d))))
    env = os.environ.copy()
    env["LEAN_PATH"] = os.pathsep.join(
        [
            str(a.modules.resolve()),
            subprocess.check_output(
                ["lake", "env", "printenv", "LEAN_PATH"], cwd=PROJECT, text=True
            ).strip(),
        ]
    )
    start = time.monotonic()
    run = subprocess.run(
        ["lean", "-j1", "-M5100", str(a.output.resolve())],
        env=env,
        cwd=PROJECT,
        text=True,
        capture_output=True,
    )
    ax = f"'sharedNormalized{a.index}' depends on axioms: [propext, Classical.choice, Quot.sound]"
    if run.returncode or ax not in run.stdout.splitlines():
        raise RuntimeError(f"concrete replay failed:\n{run.stdout}{run.stderr}")
    print(ax, f"({time.monotonic() - start:.3f} seconds)")


if __name__ == "__main__":
    main()
