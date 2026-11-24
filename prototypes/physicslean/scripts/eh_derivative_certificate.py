"""Replay signed EH derivatives using the checked-host modules of eh_sweep_certificate.py."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

if __package__:
    from .eh_sweep_certificate import add, box, literal, pair
    from .spectrum_table import PROJECT
else:
    from eh_sweep_certificate import add, box, literal, pair
    from spectrum_table import PROJECT

OPTIONS = [
    "open FloatLib.Numerics SashimiEhSweep SashimiEisensteinIntegral SashimiEhDerivativeBounds",
    "set_option maxRecDepth 100000",
    "set_option maxHeartbeats 16000000",
    "set_option Elab.async false",
]


def write_report(path, data):
    # Readers can consume completed rows while other kernels are still working.
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as output:
        json.dump(data, output, indent=2)
        output.write("\n")
    os.replace(output.name, path)


def prepare(mesh_file, derivative_file, directory, chunk_size):
    if chunk_size <= 0:
        raise ValueError("chunk size must be positive and equal the shared replay block size")
    mesh = json.loads(mesh_file.read_bytes())
    result = json.loads(derivative_file.read_bytes())
    if hashlib.sha256(mesh_file.read_bytes()).hexdigest() != result["mesh_sha256"]:
        raise ValueError("derivative and shared mesh digests differ")
    a, b = Q(mesh["a"]), Q(mesh["b"])
    bits, degree = mesh["bits"], mesh["degree"]
    rows = [
        (i, Q(r["mass_msun"]), pair(r["dS_dM"]))
        for i, r in enumerate(result["rows"])
        if r["status"] == "enclosed"
    ]
    records = [tuple(map(Q, x)) for x in mesh["bins"]]
    parts = [records[i : i + chunk_size] for i in range(0, len(records), chunk_size)]
    left = [a] + [p[-1][0] for p in parts[:-1]]
    stdin = [f"cell {r} {lo} {hi}" for r, lo, hi in records]
    stdin.extend(f"mass {m}" for _, m, _ in rows)
    run = subprocess.run(
        [
            "lake",
            "env",
            "lean",
            "--run",
            "scripts/EhDerivativeBounds.lean",
            str(a),
            str(bits),
            str(degree),
            result["absolute_width_budget"],
            str(chunk_size),
        ],
        cwd=PROJECT,
        input="\n".join(stdin) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    raw, radii, outputs = {}, {}, {}
    for line in run.stdout.splitlines():
        typ, *data = line.split()
        if typ == "part":
            raw[Q(data[0]), int(data[1])] = tuple(map(Q, data[2:]))
        elif typ == "radius":
            radii[Q(data[0])] = tuple(map(Q, data[1:]))
        elif typ == "enclosed":
            outputs[Q(data[0])] = tuple(map(Q, data[1:]))
    directory.mkdir(parents=True, exist_ok=True)

    def write(name, lines):
        (directory / f"{name}.lean").write_text("\n".join(lines) + "\n", encoding="utf-8")

    base = ["import SashimiEhDerivativeBounds", *OPTIONS]
    for i, part in enumerate(parts):
        cells = ", ".join(f"⟨{literal(r)}, {box((lo, hi))}⟩" for r, lo, hi in part)
        base.append(f"def inputBins{i} : List Bin := [{cells}]")
    write("SlopeInputs", base)
    tol = literal(Q(result["absolute_width_budget"]))
    for idx, m, enclosure in rows:
        if outputs[m] != enclosure:
            raise ValueError(f"standalone endpoints disagree at {m}")
        pref = f"slope{idx}_"
        R = f"{pref}radius"
        slope = f"SashimiEhDerivativeBounds.sum {R} {bits} {degree}"
        leaf = [
            "import SlopeInputs",
            *OPTIONS,
            f"def {R} := {box(radii[m])}",
            f"theorem {pref}radiusCorrect : radiusBox? matched {literal(m)} {bits} {degree} = some {R} := by decide +kernel",
        ]
        final = [f"import SlopePanels{idx}", "import WindowCertificate", *OPTIONS]
        bounds = {}
        starts = {}
        nodes = list(range(len(parts)))
        for i, lo in enumerate(left):
            bounds[i], starts[i] = raw[m, i], lo
            pred = f"{slope} {literal(lo)}"
            leaf.append(
                f"theorem {pref}part{i} : {pred} inputBins{i} = {box(bounds[i])} := by decide +kernel"
            )
            final.append(
                f"theorem {pref}sum{i} : {pred} ehBins{i} = {box(bounds[i])} := {pref}part{i}"
            )
        label = len(nodes)
        while len(nodes) > 1:
            merged = []
            for j in range(0, len(nodes), 2):
                x = nodes[j]
                if j + 1 == len(nodes):
                    merged.append(x)
                    continue
                y, i = nodes[j + 1], label
                label += 1
                starts[i], bounds[i] = starts[x], add(bounds[x], bounds[y])
                final.extend(
                    [
                        f"theorem {pref}sum{i} : {slope} {literal(starts[i])} ehBins{i} = {box(bounds[i])} := by",
                        f"  rw [ehBins{i}, sum_append, end{x}, {pref}sum{x}, {pref}sum{y}]; decide +kernel",
                    ]
                )
                merged.append(i)
            nodes = merged
        pred = f"deriv (fun M => SashimiEhVariance.normalizedS M {literal(a)} {literal(b)}) ({literal(m)} : ℝ)"
        if enclosure[1] < 0:
            sign = "negative"
            conclusion = f"{pred} < 0 := {pref}contains.2.trans_lt (by norm_num)"
        elif 0 < enclosure[0]:
            sign = "positive"
            conclusion = f"0 < {pred} := lt_of_lt_of_le (by norm_num) {pref}contains.1"
        else:
            raise ValueError(f"accepted derivative has indeterminate sign at {m}")
        final.extend(
            [
                f"theorem {pref}contains : {box(enclosure)}.Contains ({pred}) := by",
                f"  apply replayed_deriv (tol := {tol}) (by decide +kernel) (by decide +kernel) (by decide +kernel) "
                f"certifiedSpectrum certifiedBandEnd {pref}radiusCorrect (by decide +kernel) {pref}sum{nodes[0]} concrete.2.1",
                "  decide +kernel",
                f"theorem {pref}{sign} : {conclusion}",
                f"#print axioms {pref}{sign}",
            ]
        )
        write(f"SlopePanels{idx}", leaf)
        write(f"SlopeFinal{idx}", final)
    return {
        "mesh_sha256": result["mesh_sha256"],
        "k_band": [str(a), str(b)],
        "cells": len(records),
        "proof_chunk_size": chunk_size,
        "unresolved_masses": [r["mass_msun"] for r in result["rows"] if r["status"] != "enclosed"],
        "masses": [
            {"index": i, "mass_msun": str(m), "dS_dM": {"lo": str(v[0]), "hi": str(v[1])}}
            for i, m, v in rows
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["prepare", "panels", "finish", "report"])
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, help="destination of the grid replay report")
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--derivative", type=Path)
    parser.add_argument(
        "--shared",
        type=Path,
        help="module-output from eh_sweep_certificate.py for the same mesh and block size",
    )
    parser.add_argument("--chunk-size", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--memory-mb", type=int, default=4600)
    parser.add_argument(
        "--exclude",
        type=Q,
        action="append",
        default=[],
        help="masses checked separately, for distributing the numerical replay",
    )
    parser.add_argument("--mass", type=Q, action="append", help="check only these exact masses")
    args = parser.parse_args()
    directory = args.directory.resolve()
    manifest_path = directory / "derivative-manifest.json"
    if args.phase == "prepare":
        if args.mesh is None or args.derivative is None:
            parser.error("prepare requires --mesh and --derivative")
        write_report(manifest_path, prepare(args.mesh, args.derivative, directory, args.chunk_size))
        return
    data = json.loads(manifest_path.read_text())
    if args.phase == "report":
        if args.output is None:
            parser.error("report requires --output")
        rows = []
        for row in data["masses"]:
            idx = row["index"]
            proofs = [directory / f"{prefix}{idx}.json" for prefix in ("SlopePanels", "SlopeFinal")]
            item = {**row, "kernel_replayed": False}
            if all(path.exists() for path in proofs):
                parts = [json.loads(path.read_text()) for path in proofs]
                for proof, name in zip(parts, ("panels", "finish")):
                    if (
                        any(proof[key] != row[key] for key in row)
                        or proof["mesh_sha256"] != data["mesh_sha256"]
                        or proof["phase"] != name
                    ):
                        raise ValueError(f"stale or inconsistent replay at mass {row['mass_msun']}")
                item.update(
                    kernel_replayed=True,
                    panels_seconds=parts[0]["seconds"],
                    composition_seconds=parts[1]["seconds"],
                )
            rows.append(item)
        write_report(
            args.output,
            {
                "target": "signed derivative of finite-band normalizedS by physical mass",
                "mesh_sha256": data["mesh_sha256"],
                "k_band": data["k_band"],
                "cells": data["cells"],
                "proof_chunk_size": data["proof_chunk_size"],
                "unresolved_masses": data.get("unresolved_masses", []),
                "kernel_replayed": bool(rows)
                and not data.get("unresolved_masses")
                and all(r["kernel_replayed"] for r in rows),
                "replayed_count": sum(r["kernel_replayed"] for r in rows),
                "rows": rows,
            },
        )
        return
    package = subprocess.check_output(
        ["lake", "env", "printenv", "LEAN_PATH"], cwd=PROJECT, text=True
    ).strip()
    paths = [str(directory), package]
    if args.shared:
        paths.insert(1, str(args.shared.resolve()))
    env = dict(os.environ, LEAN_PATH=os.pathsep.join(paths))

    def compile_one(name):
        cmd = [
            "lean",
            "-j1",
            f"-M{args.memory_mb}",
            "-R",
            str(directory),
            "-o",
            str(directory / f"{name}.olean"),
            str(directory / f"{name}.lean"),
        ]
        start = time.monotonic()
        proc = subprocess.run(cmd, cwd=PROJECT, env=env, capture_output=True, text=True)
        output = proc.stdout + proc.stderr
        (directory / f"{name}.log").write_text(output)
        if proc.returncode:
            raise RuntimeError(f"{name}: {output[-1200:]}")
        return time.monotonic() - start, output

    if args.phase == "panels":
        compile_one("SlopeInputs")
    elif not args.shared:
        parser.error("finish requires --shared")
    rows = [r for r in data["masses"] if not args.mass or Q(r["mass_msun"]) in args.mass]
    rows = [r for r in rows if Q(r["mass_msun"]) not in args.exclude]
    if args.mass and len(rows) != len(set(args.mass)):
        parser.error("unknown mass in requested grid")

    def verify(row):
        idx = row["index"]
        prefix = "SlopePanels" if args.phase == "panels" else "SlopeFinal"
        stamp = directory / f"{prefix}{idx}.json"
        stamp.unlink(missing_ok=True)
        if args.phase == "panels":
            (directory / f"SlopeFinal{idx}.json").unlink(missing_ok=True)
        seconds, output = compile_one(f"{prefix}{idx}")
        if args.phase == "finish":
            sign = "negative" if Q(row["dS_dM"]["hi"]) < 0 else "positive"
            expected = (
                f"'slope{idx}_{sign}' depends on axioms: [propext, Classical.choice, Quot.sound]"
            )
            if expected not in output:
                raise RuntimeError("missing standard-axiom audit")
        report = dict(row, phase=args.phase, seconds=seconds, mesh_sha256=data["mesh_sha256"])
        write_report(stamp, report)
        print(f"{args.phase} {row['mass_msun']} {seconds:.3f}s", flush=True)

    with ThreadPoolExecutor(max_workers=args.jobs) as workers:
        list(workers.map(verify, reversed(rows)))


if __name__ == "__main__":
    main()
