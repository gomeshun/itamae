"""Replay finite-band analytic EH98b power cells and normalized mass variance in Lean.

Both JSON inputs are untrusted. The proof recomputes every transfer range,
both windowed integrals, the mass radius and the normalized endpoints.
"""

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic

if __package__:
    from .table_replay import literal
    from .table_window import PROJECT, check_modules
else:
    from table_replay import literal
    from table_window import PROJECT, check_modules


def box(bounds):
    lo, hi = bounds
    return f"(⟨{literal(lo)}, {literal(hi)}⟩ : RationalInterval)"


def pair(obj):
    return Q(obj["lo"]), Q(obj["hi"])


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def report(bins_path, sweep_path, mass, *, chunk_size=8, jobs=1, timeout=3600, module_output=None):
    """Verify weighted transfer cells and both integrals for one physical mass in Msun."""
    if min(chunk_size, jobs, timeout) <= 0:
        raise ValueError("chunk size, jobs and timeout must be positive")
    data = json.loads(Path(bins_path).read_text())
    result = json.loads(Path(sweep_path).read_text())
    bits, degree = int(data["bits"]), int(data["degree"])
    cb, cd = int(data["coefficient_bits"]), int(data["coefficient_degree"])
    a, b = Q(data["a"]), Q(data["b"])
    records = [[Q(x) for x in row] for row in data["bins"]]
    edges = [a, *(right for right, _, _ in records)]
    if (
        not records
        or edges[-1] != b
        or any(lower >= upper for lower, upper in zip(edges, edges[1:]))
    ):
        raise ValueError("Expected an increasing grid with last endpoint equal to b")
    mass = Q(mass)
    row = next(r for r in result["rows"] if Q(r["mass_msun"]) == mass)
    R = pair(row["radius"])
    chunks = [records[i : i + chunk_size] for i in range(0, len(records), chunk_size)]
    left = [a] + [ch[-1][0] for ch in chunks[:-1]]
    header = [
        "import SashimiEhCertificate",
        "open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhSweep SashimiEhCertificate",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 30000000",
        "set_option Elab.async false",
        f"def ehRadius := {box(R)}",
    ]
    definitions = []
    for i, ch in enumerate(chunks):
        cells = ", ".join("⟨" + literal(right) + "," + box((lo, hi)) + "⟩" for right, lo, hi in ch)
        definitions.append(f"def ehBins{i} : List Bin := [{cells}]")
    runtime = [
        *header,
        *definitions,
        "def main : IO Unit := do",
        f'  let some cs := boxes? matched {cb} {cd} | throw (IO.userError "coefficient enclosure failed")',
        '  for c in cs do IO.println s!"COEFF {c.lo} {c.hi}"',
    ]
    for i, lower_k in enumerate(left):
        for field, rad in (("raw", "ehRadius"), ("anchor", "(.point 8)")):
            runtime += [
                f"  let v := total {rad} {bits} {degree} {literal(lower_k)} ehBins{i}",
                f'  IO.println s!"SUB {i} {field} {{v.lo}} {{v.hi}}"',
            ]
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "MomentBlocks.lean"
        p.write_text("\n".join(runtime), encoding="utf-8")
        run = subprocess.run(
            ["lake", "env", "lean", "--run", str(p)], cwd=PROJECT, text=True, capture_output=True
        )
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    coeffs = []
    subtotals = {}
    for line in run.stdout.splitlines():
        if line.startswith("COEFF "):
            _, lo, hi = line.split()
            coeffs.append((Q(lo), Q(hi)))
        elif line.startswith("SUB "):
            _, i, k, lo, hi = line.split()
            subtotals[(int(i), k)] = (Q(lo), Q(hi))
    header += [
        "def ehCoefficients : List RationalInterval := [" + ",".join(box(x) for x in coeffs) + "]",
        f"theorem hEhCoefficients : boxes? matched {cb} {cd} = some ehCoefficients := by decide +kernel",
        f"theorem hEhRadius : radiusBox? matched {literal(mass)} {bits} {degree} = some ehRadius := by decide +kernel",
    ]
    blocks = []
    nodes = []
    bounds = {}
    starts = {}
    ends = {}
    for i, (ch, lo) in enumerate(zip(chunks, left)):
        starts[i], ends[i] = lo, ch[-1][0]
        bins = f"ehBins{i}"
        lit = literal(lo)
        hi = literal(ends[i])
        group = [
            definitions[i],
            f"theorem good{i} : Good (weightReal matched) {lit} {bins} :=",
            f"  build_sound (ends := {bins}.map (·.right)) (bits := {bits}) (degree := {degree}) hEhCoefficients (by decide +kernel)",
            f"theorem end{i} : last {lit} {bins} = {hi} := by decide +kernel",
        ]
        for key, r in (("raw", "ehRadius"), ("anchor", "(.point 8)")):
            bounds[i, key] = subtotals[i, key]
            group.append(
                f"theorem {key}{i} : total {r} {bits} {degree} {lit} {bins} = {box(bounds[i, key])} := by decide +kernel"
            )
        blocks.append(group)
        nodes.append(i)
    assembly = []
    label = len(nodes)
    while len(nodes) > 1:
        merged = []
        for j in range(0, len(nodes), 2):
            x = nodes[j]
            if j + 1 == len(nodes):
                merged.append(x)
                continue
            y = nodes[j + 1]
            i = label
            label += 1
            starts[i], ends[i] = starts[x], ends[y]
            lo = literal(starts[i])
            hi = literal(ends[i])
            bins = f"ehBins{i}"
            assembly.extend(
                [
                    f"def {bins} := ehBins{x} ++ ehBins{y}",
                    f"theorem good{i} : Good (weightReal matched) {lo} {bins} := good_append good{x} (end{x} ▸ good{y})",
                    f"theorem end{i} : last {lo} {bins} = {hi} := by rw [{bins}, last_append, end{x}, end{y}]",
                ]
            )
            for key, r in (("raw", "ehRadius"), ("anchor", "(.point 8)")):
                bounds[i, key] = add(bounds[x, key], bounds[y, key])
                assembly.extend(
                    [
                        f"theorem {key}{i} : total {r} {bits} {degree} {lo} {bins} = {box(bounds[i, key])} := by",
                        f"  rw [{bins}, total_append, end{x}, {key}{x}, {key}{y}] <;> decide +kernel",
                    ]
                )
            merged.append(i)
        nodes = merged
    node = nodes[0]
    N, A, S = pair(row["raw"]), pair(result["anchor"]["raw"]), pair(row["S"])
    if not (bounds[node, "raw"] == N and bounds[node, "anchor"] == A):
        raise ValueError("Transport subtotals or band do not match supplied report")
    raw = f"raw matched {literal(a)} {literal(b)}"
    ns = f"{raw} (radius matched {literal(mass)})"
    aa = f"{raw} 8"
    assembly.extend(
        [
            f"abbrev certifiedBins := ehBins{node}",
            f"theorem certifiedSpectrum : Good (weightReal matched) {literal(a)} certifiedBins := good{node}",
            f"theorem certifiedBandEnd : last {literal(a)} certifiedBins = {literal(b)} := end{node}",
            f"theorem certifiedAnchor : total (.point 8) {bits} {degree} {literal(a)} certifiedBins = {box(A)} := anchor{node}",
            f"theorem concrete : {box(N)}.Contains ({ns}) ∧ {box(A)}.Contains ({aa}) ∧",
            f"  {box(S)}.Contains (SashimiSpectrum.normalize8 matched.sigma8 ({aa}) ({ns})) := by",
            f"  have outcome := SashimiEhCertificate.report good{node} end{node} hEhRadius (by decide +kernel) raw{node} anchor{node} (by decide +kernel)",
            f"  have e : ({box(N)}.scaleNonnegative (matched.sigma8^2)).divPositive {box(A)} = {box(S)} := by decide +kernel",
            "  rwa [e] at outcome",
            "#print axioms concrete",
        ]
    )
    start = monotonic()
    axioms = check_modules(
        header, blocks, assembly, 5800, timeout, jobs, module_output=module_output
    )
    return {
        "kernel_replayed": True,
        "kernel_axioms": axioms,
        "mass_msun": str(mass),
        "k_band": [str(a), str(b)],
        "shape_cells": len(records),
        "proof_block_size": chunk_size,
        "bins_sha256": hashlib.sha256(Path(bins_path).read_bytes()).hexdigest(),
        "bits": bits,
        "degree": degree,
        "raw": row["raw"],
        "anchor": result["anchor"]["raw"],
        "S": row["S"],
        "kernel_seconds": monotonic() - start,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("bins", type=Path)
    p.add_argument("sweep", type=Path)
    p.add_argument("--mass", type=Q, default=Q(10**12))
    p.add_argument("--chunk-size", type=int, default=8)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--timeout", type=int, default=3600)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--module-output",
        type=Path,
        help="save checked Lean sources and olean files for further composition",
    )
    args = p.parse_args()
    data = report(
        args.bins,
        args.sweep,
        args.mass,
        chunk_size=args.chunk_size,
        jobs=args.jobs,
        timeout=args.timeout,
        module_output=args.module_output,
    )
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
