"""Replay each requested analytic EH mass using a shared, kernel-checked spectrum.

First generate a shared module with eh_sweep_certificate.py --module-output.
This script checks radius, cell sum and normalization for every selected mass.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

if __package__:
    from .eh_sweep_certificate import add, box, pair
    from .table_replay import literal
else:
    from eh_sweep_certificate import add, box, pair
    from table_replay import literal

PROJECT = Path(__file__).resolve().parents[1]
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}


def _chunks(modules):
    starts = {}
    for group in modules.glob("WindowGroup*.lean"):
        for i, q in re.findall(
            r"^theorem end(\d+) : last \(([^ ]+) : ℚ\) ehBins\d+ = ", group.read_text(), re.M
        ):
            starts[int(i)] = Fraction(q)
    joins = [
        tuple(map(int, m))
        for m in re.findall(
            r"^def ehBins(\d+) := ehBins(\d+) \+\+ ehBins(\d+)$",
            (modules / "WindowCertificate.lean").read_text(),
            re.M,
        )
    ]
    if set(starts) != set(range(len(starts))) or not starts:
        raise ValueError("No complete certified leaf list")
    return starts, joins


def _subtotals(output, env, ids, rows, starts, bits, degree, *, paired_parent=None, depth=0):
    chunks = ", ".join(f"({literal(starts[j])}, ehBins{j})" for j in sorted(starts))
    radii = ", ".join(f"({i}, {box(pair(rows[i]['radius']))})" for i in ids)
    evaluation = f"total r {bits} {degree} a cells"
    imports = ""
    if paired_parent is not None:
        imports = "import SashimiEhYang"
        evaluation = f"SashimiEhYang.pairedTotal r {box(pair(paired_parent))} {bits} {degree} {depth} a cells"
    src = f"""import WindowCertificate
{imports}
open FloatLib.Numerics SashimiEhSweep
set_option maxRecDepth 100000
def chunks : List (ℚ × List Bin) := [{chunks}]
def radii : List (Nat × RationalInterval) := [{radii}]
def main : IO Unit := do
  for (i, r) in radii do
    for ((a, cells), j) in chunks.zipIdx do
      let v := {evaluation}
      IO.println s!"{{i}} {{j}} {{v.lo}} {{v.hi}}"
"""
    path = output / f"MassSubtotals{ids[0]}.lean"
    path.write_text(src)
    r = subprocess.run(
        ["lean", "--run", str(path.resolve())],
        cwd=PROJECT,
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
    )
    if r.returncode:
        raise RuntimeError(r.stdout + r.stderr)
    results = {i: {} for i in ids}
    for line in r.stdout.splitlines():
        i, j, lo, hi = line.split()
        results[int(i)][int(j)] = Fraction(lo), Fraction(hi)
    return results


def _program(i, row, a, b, A, bits, degree, starts, joins, values):
    m, R, raw, S = (
        literal(Fraction(row["mass_msun"])),
        box(pair(row["radius"])),
        box(pair(row["raw"])),
        box(pair(row["S"])),
    )
    proofs = []
    starts = dict(starts)
    for j in sorted(starts):
        proofs.append(
            f"private theorem sub{j} : total rowRadius {bits} {degree} {literal(starts[j])} ehBins{j} = {box(values[j])} := by decide +kernel"
        )
    for j, x, y in joins:
        starts[j] = starts[x]
        values[j] = add(values[x], values[y])
        proofs.append(
            f"private theorem sub{j} : total rowRadius {bits} {degree} {literal(starts[j])} ehBins{j} = {box(values[j])} := by\n  rw [ehBins{j}, total_append, end{x}, sub{x}, sub{y}]; decide +kernel"
        )
    last = max(values)
    if pair(row["raw"]) != values[last]:
        raise ValueError(f"Inconsistent exact subtotals for mass index {i}")
    assembled = "\n".join(proofs)
    return f"""import WindowCertificate
import SashimiEhSweep
open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhSweep SashimiEhCertificate
set_option maxRecDepth 100000
set_option maxHeartbeats 30000000
set_option Elab.async false
private theorem radiusCalc{i} : radiusBox? matched {m} {bits} {degree} = some {R} := by decide +kernel
theorem ehRadius{i} : {R}.Contains (radius matched {m}) := radiusBox_sound radiusCalc{i}
private def rowRadius : RationalInterval := {R}
{assembled}
private theorem rawCalc{i} : total {R} {bits} {degree} {a} certifiedBins = {raw} := sub{last}
theorem ehRaw{i} : {raw}.Contains (raw matched {a} {b} (radius matched {m})) := by
  have h := total_sound (bits := {bits}) (degree := {degree}) certifiedSpectrum ehRadius{i}
    (by decide +kernel : 0 ≤ {R}.lo)
  rwa [rawCalc{i}, certifiedBandEnd] at h
theorem ehMass{i} : {S}.Contains (SashimiEhVariance.normalizedS {m} {a} {b}) := by
  have h := total_sound (bits := {bits}) (degree := {degree}) certifiedSpectrum
    (show (RationalInterval.point 8).Contains (8:ℝ) from by norm_num [RationalInterval.Contains,RationalInterval.point])
    (by decide +kernel : 0 ≤ (RationalInterval.point 8).lo)
  rw [certifiedBandEnd,certifiedAnchor] at h
  have hn : SashimiNormalization.apply? matched.sigma8 (10^20) {raw} {A} = some {S} := by decide +kernel
  have hz := (SashimiNormalization.apply_sound ehRaw{i} h hn).1
  rw [normalized_eq (by decide +kernel : 0 < {a}) (by decide +kernel : {a} ≤ {b})] at hz
  exact hz
#print axioms ehMass{i}
"""


def replay(sweep, modules, output, ids, jobs=1):
    source = json.loads(sweep.read_text())
    rows = [*source["rows"], source["host_progenitor"]]
    if not ids or len(set(ids)) != len(ids) or min(ids) < 0 or max(ids) >= len(rows):
        raise ValueError("invalid or duplicate mass index")
    if not (modules / "WindowCertificate.olean").is_file():
        raise ValueError("build WindowCertificate.olean before replaying masses")
    p = source["mesh"]
    a, b = (literal(Fraction(x)) for x in source["k_band"])
    A = box(pair(source["anchor"]["raw"]))
    output.mkdir(parents=True, exist_ok=True)
    pkg = subprocess.check_output(
        ["lake", "env", "printenv", "LEAN_PATH"], cwd=PROJECT, text=True
    ).strip()
    env = dict(
        os.environ, LEAN_PATH=os.pathsep.join([str(output.resolve()), str(modules.resolve()), pkg])
    )

    starts, joins = _chunks(modules)
    values = _subtotals(output, env, ids, rows, starts, p["bits"], p["degree"])

    def run(i):
        name = f"EhMass{i}"
        begin = time.monotonic()
        path = output / (name + ".lean")
        obj = output / (name + ".olean")
        path.write_text(
            _program(i, rows[i], a, b, A, p["bits"], p["degree"], starts, joins, values[i])
        )
        cmd = [
            "lean",
            "-j1",
            "-M4096",
            "-R",
            str(output.resolve()),
            "-o",
            str(obj.resolve()),
            str(path.resolve()),
        ]
        r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=3600, cwd=PROJECT)
        (output / (name + ".log")).write_text(r.stdout + r.stderr)
        audited = re.search(rf"'ehMass{i}' depends on axioms: \[(.*?)\]", r.stdout)
        if r.returncode or not audited or set(audited[1].split(", ")) - ALLOWED:
            raise RuntimeError((r.stdout + r.stderr)[-3600:])
        return {
            "index": i,
            "mass_msun": rows[i]["mass_msun"],
            "seconds": time.monotonic() - begin,
            "axioms": audited[1].split(", "),
        }

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        outcomes = list(pool.map(run, ids))
    return {
        "sweep_sha256": hashlib.sha256(sweep.read_bytes()).hexdigest(),
        "kernel_replay": "passed",
        "rows": outcomes,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("sweep", "modules", "output"):
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument(
        "--indices", help="comma-separated indices; default all mass rows plus progenitor"
    )
    p.add_argument("--jobs", type=int, default=1)
    args = p.parse_args()
    if args.jobs < 1:
        raise ValueError("jobs must be positive")
    ids = (
        [int(s) for s in args.indices.split(",")]
        if args.indices
        else list(range(len(json.loads(args.sweep.read_text())["rows"]) + 1))
    )
    result = replay(args.sweep, args.modules, args.output, ids, args.jobs)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
