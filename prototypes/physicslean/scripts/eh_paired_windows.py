"""Replay signed window differences on previously checked analytic EH spectrum bins.

WindowCertificate supplies the spectral inputs and normalized radius-8 anchor.
One EhPaired{index} module per child proves the full variance difference for
any radii in its input boxes. This permits an exact fitted host progenitor with
its own radius proof. Four window subpanels per spectral bin recover the nearest
three progenitor masses on the 960-cell grid.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess
import time

if __package__:
    from .eh_sweep_certificate import add, box, pair
    from .eh_sweep_masses import PROJECT, _chunks, _subtotals
    from .eh_full_grid import budget
    from .table_replay import literal
else:
    from eh_sweep_certificate import add, box, pair
    from eh_sweep_masses import PROJECT, _chunks, _subtotals
    from eh_full_grid import budget
    from table_replay import literal


def normalized_gap(raw_gap, child, parent, anchor):
    lower, upper = raw_gap
    lower -= budget(parent)
    upper += budget(child)
    al, au = anchor
    au += budget(Fraction(8))
    if al <= 0:
        raise ValueError("positive anchor required")
    candidates = [n / d for n in (lower, upper) for d in (al, au)]
    return tuple(Fraction(811, 1000) ** 2 * f(candidates) for f in (min, max))


def module(i, row, parent, anchor, values, starts, joins, bits, degree, depth):
    child_r, parent_r = pair(row["radius"]), pair(parent["radius"])
    chunks = dict(values)
    starts = dict(starts)
    lines = [
        "import WindowCertificate",
        "import SashimiEhYang",
        "open SashimiEhYang SashimiEhSweep SashimiEisensteinIntegral SashimiEhFullGrid FloatLib.Numerics",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 60000000",
        "set_option Elab.async false",
        f"private def child : RationalInterval := {box(child_r)}",
        f"private def host : RationalInterval := {box(parent_r)}",
        f"private def a : RationalInterval := {box(anchor)}",
        "private theorem low_eq : lowEdge = (1/100000000 : ℚ) := by decide +kernel",
        "private theorem high_eq : highEdge = (10000000000/1 : ℚ) := by decide +kernel",
    ]
    prefix = f"pairedTotal child host {bits} {degree} {depth}"
    for j in sorted(starts):
        lines.append(
            f"private theorem sum{j} : {prefix} {literal(starts[j])} ehBins{j} = {box(chunks[j])} := by decide +kernel"
        )
    for j, left, right in joins:
        starts[j] = starts[left]
        chunks[j] = add(chunks[left], chunks[right])
        lines.append(
            f"private theorem sum{j} : {prefix} {literal(starts[j])} ehBins{j} = {box(chunks[j])} := by\n"
            f"  rw [ehBins{j}, paired_append, end{left}, sum{left}, sum{right}]; decide +kernel"
        )
    last = max(chunks)
    J = chunks[last]
    delta = normalized_gap(J, child_r[0], parent_r[0], anchor)
    lines.extend(
        [
            f"theorem ehPaired{i} {{x y : ℝ}} (hx : {box(child_r)}.Contains x) (hy : {box(parent_r)}.Contains y) : "
            f"{box(J)}.Contains (raw matched lowEdge highEdge x - raw matched lowEdge highEdge y) := by",
            f"  have h := pairedTotal_sound (c := child) (p := host) (bits := {bits}) (degree := {degree}) (depth := {depth}) certifiedSpectrum certifiedBandEnd low_eq.symm high_eq.symm hx hy (by decide +kernel) (by decide +kernel)",
            "  rw [low_eq,high_eq]",
            f"  rw [show {prefix} (1/100000000 : ℚ) certifiedBins = {box(J)} from sum{last}] at h",
            "  exact h",
            f"theorem ehPairedDelta{i} {{x y : ℝ}} (hx : {box(child_r)}.Contains x) (hy : {box(parent_r)}.Contains y) : "
            f"{box(delta)}.Contains (SashimiSpectrum.normalize8 matched.sigma8 (SashimiEhFull.fullRaw matched 8) (SashimiEhFull.fullRaw matched x) - "
            "SashimiSpectrum.normalize8 matched.sigma8 (SashimiEhFull.fullRaw matched 8) (SashimiEhFull.fullRaw matched y)) := by",
            "  have ha : a.Contains (raw matched lowEdge highEdge 8) := by rw [low_eq,high_eq]; exact concrete.2.1",
            f"  have h := pairedDifference_sound (rc := child.lo) (rp := host.lo) hx.1 hy.1 (ehPaired{i} hx hy) ha (by decide +kernel)",
            f"  have he : pairedDifference child.lo host.lo {box(J)} a = {box(delta)} := by decide +kernel",
            "  rwa [he] at h",
            f"#print axioms ehPairedDelta{i}",
        ]
    )
    return "\n".join(lines) + "\n", J, delta


def run(source, modules, dest, ids, depth, jobs):
    sweep = json.loads(source.read_text())
    rows = [*sweep["rows"], sweep["host_progenitor"]]
    if (
        not ids
        or len(ids) != len(set(ids))
        or min(ids) < 0
        or max(ids) >= len(rows)
        or depth < 0
        or jobs < 1
    ):
        raise ValueError("invalid child indices, depth or job count")
    parent = sweep["rows"][-1]
    if Fraction(parent["mass_msun"]) != 10**12 or list(map(Fraction, sweep["k_band"])) != [
        Fraction("1e-8"),
        Fraction("1e10"),
    ]:
        raise ValueError("expected the matched low/high spectral edges and physical host mass")
    starts, joins = _chunks(modules)
    dest.mkdir(parents=True, exist_ok=True)
    lpath = subprocess.check_output(
        ["lake", "env", "printenv", "LEAN_PATH"], cwd=PROJECT, text=True
    ).strip()
    env = dict(os.environ, LEAN_PATH=os.pathsep.join((str(modules.resolve()), lpath)))
    bits, degree = sweep["mesh"]["bits"], sweep["mesh"]["degree"]
    totals = _subtotals(
        dest, env, ids, rows, starts, bits, degree, paired_parent=parent["radius"], depth=depth
    )

    def one(i):
        src, J, d = module(
            i,
            rows[i],
            parent,
            pair(sweep["anchor"]["raw"]),
            totals[i],
            starts,
            joins,
            bits,
            degree,
            depth,
        )
        name = f"EhPaired{i}"
        path = dest / (name + ".lean")
        obj = dest / (name + ".olean")
        path.write_text(src)
        start = time.monotonic()
        r = subprocess.run(
            [
                "lean",
                "-j1",
                "-M5100",
                "-R",
                str(dest.resolve()),
                "-o",
                str(obj.resolve()),
                str(path.resolve()),
            ],
            env=env,
            cwd=PROJECT,
            capture_output=True,
            text=True,
            timeout=1600,
        )
        log = dest / (name + ".log")
        log.write_text(r.stdout + r.stderr)
        if (
            r.returncode
            or f"'ehPairedDelta{i}' depends on axioms: [propext, Classical.choice, Quot.sound]"
            not in r.stdout.splitlines()
        ):
            raise RuntimeError(f"paired integral/normalization proof failed; see {log}")
        return {
            "index": i,
            "mass_msun": rows[i]["mass_msun"],
            "radius": rows[i]["radius"],
            "parent_radius": parent["radius"],
            "raw_difference": {"lo": str(J[0]), "hi": str(J[1])},
            "delta_S": {"lo": str(d[0]), "hi": str(d[1])},
            "module": name,
            "theorem": f"ehPairedDelta{i}",
            "kernel_replayed": True,
            "replay_seconds": time.monotonic() - start,
        }

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        result = list(pool.map(one, ids))
    return {
        "scope": "Full EH normalized variance difference, quantified over both radius boxes; tails are bounded with their separate radius-dependent budgets.",
        "shape_cells": sweep["mesh"]["cells"],
        "window_subpanels_per_cell": 2**depth,
        "rows": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("sweep", "modules", "output_modules", "output"):
        parser.add_argument(arg, type=Path)
    parser.add_argument("--indices", type=int, nargs="+", default=[69, 70, 71, 73])
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=1)
    a = parser.parse_args()
    data = run(a.sweep, a.modules, a.output_modules, a.indices, a.depth, a.jobs)
    a.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
