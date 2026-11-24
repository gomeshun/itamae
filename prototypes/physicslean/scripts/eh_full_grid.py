"""Promote matched EH finite-band raw intervals to the full (0,∞) variance.

The input raw integrals omit 1/(2π²), which cancels in σ₈ normalization.
Replay the exact propagation for each row. Optional per-mass independent
WindowCertificate modules discharge both finite-band integrals and the radius.
"""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import isqrt
import os
from pathlib import Path
import subprocess

PROJECT = Path(__file__).resolve().parents[1]
LOW, HIGH = Q(1, 10**8), Q(10**10)
C0, C1 = Q(3), Q(1, 10**9)


def budget(rmin):
    if rmin <= 0 or HIGH * rmin < 1:
        raise ValueError("Tail requires positive radius and K*rmin ≥ 1")
    return C0**2 * LOW**3 + 18 * C1**2 / (rmin**4 * HIGH**2)


def pair(box):
    low, high = Q(box["lo"]), Q(box["hi"])
    if low > high:
        raise ValueError("Inverted interval")
    return low, high


def box(lo, hi):
    return {"lo": str(lo), "hi": str(hi)}


def rounded(lo, hi, grid):
    return box(Q((lo * grid).__floor__(), grid), Q((hi * grid).__ceil__(), grid))


def sqrt_grid(x, grid, upper):
    n = isqrt((x.numerator * grid * grid) // x.denominator)
    return Q(n + int(upper and n * n * x.denominator < x.numerator * grid**2), grid)


def promote(raw, anchor, radius, bits):
    lo, hi = pair(raw)
    denlo, denhi = pair(anchor)
    if denlo <= 0 or lo <= 0:
        raise ValueError("Positive numerator and anchor lower bounds required")
    b = budget(radius)
    coef = Q(811, 1000) ** 2
    scale = 2**bits
    variance = rounded(coef * lo / (denhi + budget(Q(8))), coef * (hi + b) / denlo, scale)
    sl, su = pair(variance)
    return {
        "S": variance,
        "sigma": box(sqrt_grid(sl, scale, False), sqrt_grid(su, scale, True)),
        "raw_tail_budget": str(b),
    }


def literal(x):
    value = Q(x)
    return f"({value.numerator}/{value.denominator} : ℚ)"


def interval_expr(interval):
    lo, hi = pair(interval)
    return f"(⟨{literal(lo)}, {literal(hi)}⟩ : RationalInterval)"


def statement(row):
    return (
        f"{interval_expr(row['S'])}.Contains (S {literal(row['mass_msun'])}) ∧ "
        f"{interval_expr(row['sigma'])}.Contains (Real.sqrt (S {literal(row['mass_msun'])}))"
    )


def proof(row, anchor, i, bits, degree, concrete=False, mass_module=False):
    mass = literal(row["mass_msun"])
    r = literal(row["radius"]["lo"])
    j, a, v, s = map(interval_expr, (row["raw"], anchor, row["S"], row["sigma"]))
    args = f"(m := {mass}) (rmin := {r}) (j := {j}) (a := {a}) (v := {v}) (sigma := {s})"
    if concrete:
        rad = f"ehRadius{i}.1" if mass_module else "(radiusBox_sound hEhRadius).1"
        numerator = f"ehRaw{i}" if mass_module else "concrete.1"
        return (
            f"theorem fullComplete{i} : {statement(row)} :=\n"
            f"  check_sound {args} {rad} (by rw [lowEdge_eq, highEdge_eq]; exact {numerator}) "
            f"(by rw [lowEdge_eq, highEdge_eq]; exact concrete.2.1) (by decide +kernel)\n"
            f"#print axioms fullComplete{i}\n"
        )
    return (
        f"theorem fullRadius{i} : radiusBox? matched {mass} {bits} {degree} = some "
        f"{interval_expr(row['radius'])} := by decide +kernel\n"
        f"theorem fullPropagation{i}\n"
        f"    (hn : {j}.Contains (raw matched lowEdge highEdge (radius matched {mass})))\n"
        f"    (ha : {a}.Contains (raw matched lowEdge highEdge 8)) : {statement(row)} :=\n"
        f"  check_sound {args} (radiusBox_sound fullRadius{i}).1 hn ha (by decide +kernel)\n"
    )


def check_source(path, text, modules=None):
    path.write_text(text)
    env = os.environ.copy()
    if modules is None:
        cmd = ["lake", "env", "lean", "-j1", str(path.resolve())]
    else:
        current = subprocess.check_output(
            ["lake", "env", "printenv", "LEAN_PATH"], cwd=PROJECT, text=True
        ).rstrip()
        env["LEAN_PATH"] = os.pathsep.join([*(str(path.resolve()) for path in modules), current])
        cmd = [
            subprocess.check_output(
                ["lake", "env", "which", "lean"], cwd=PROJECT, text=True
            ).strip(),
            "-j1",
            str(path.resolve()),
        ]
    run = subprocess.run(cmd, cwd=PROJECT, env=env, capture_output=True, text=True)
    path.with_suffix(".log").write_text(run.stdout + run.stderr)
    if run.returncode:
        raise RuntimeError(f"Lean failed; see {path.with_suffix('.log')}")
    if modules is not None:
        reports = [line for line in run.stdout.splitlines() if "depends on axioms" in line]
        allowed = {"propext", "Classical.choice", "Quot.sound"}
        expected = text.count("#print axioms fullComplete")
        if len(reports) != expected or any(
            set(r.partition("[")[2].rstrip("]").split(", ")) - allowed for r in reports
        ):
            raise ValueError("Missing or unexpected concrete theorem axiom dependencies")


def evaluate(source, output, bits=80, certificates=None, mass_modules=None, anchor_modules=None):
    data = json.loads(source.read_text())
    if list(map(Q, data["k_band"])) != [LOW, HIGH] or pair(data["anchor"]["radius"]) != (8, 8):
        raise ValueError("Expected matched finite band and radius-8 anchor")
    if bits < 0:
        raise ValueError("Nonnegative endpoint precision required")
    anchor = data["anchor"]["raw"]
    rb, rd = int(data["mesh"]["bits"]), int(data["mesh"]["degree"])
    rows = []
    for raw in [*data["rows"], data["host_progenitor"]]:
        rows.append(
            {key: raw[key] for key in ("mass_msun", "radius", "raw")}
            | promote(raw["raw"], anchor, pair(raw["radius"])[0], bits)
            | {"kernel_replayed": False}
        )
    shared = "open SashimiEhFullGrid SashimiEisensteinIntegral FloatLib.Numerics\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\n"
    bridge = (
        f"private theorem lowEdge_eq : lowEdge = {literal(LOW)} := by decide +kernel\n"
        f"private theorem highEdge_eq : highEdge = {literal(HIGH)} := by decide +kernel\n"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    check_source(
        output.with_suffix(".lean"),
        "import SashimiEhFullGrid\n"
        + shared
        + "\n".join(proof(row, anchor, i, rb, rd) for i, row in enumerate(rows))
        + "\n#print axioms check_sound\n",
    )
    for mass, directory in (certificates or {}).items():
        i = next(i for i, row in enumerate(rows) if Q(row["mass_msun"]) == mass)
        witness = output.with_name(output.stem + f"-mass-{i}.lean")
        check_source(
            witness,
            "import SashimiEhFullGrid\nimport WindowCertificate\n"
            + shared
            + bridge
            + proof(rows[i], anchor, i, rb, rd, concrete=True),
            [directory],
        )
        rows[i]["kernel_replayed"] = True
        rows[i]["witness"] = witness.name
    if mass_modules is not None:
        if anchor_modules is None:
            raise ValueError("Supply --anchor-modules for the shared radius-8 proof")
        available = [i for i in range(len(rows)) if (mass_modules / f"EhMass{i}.olean").is_file()]
        if not available:
            raise ValueError("No checked EhMass modules found")
        witness = output.with_name(output.stem + "-closed.lean")
        header = "import SashimiEhFullGrid\nimport WindowCertificate\n"
        header += "".join(f"import EhMass{i}\n" for i in available)
        bodies = [
            proof(rows[i], anchor, i, rb, rd, concrete=True, mass_module=True) for i in available
        ]
        check_source(
            witness, header + shared + bridge + "\n".join(bodies), [mass_modules, anchor_modules]
        )
        for i in available:
            rows[i]["kernel_replayed"] = True
            rows[i]["witness"] = witness.name
    return {
        "model": "Flat matched baryonic EH98 linear transfer, z=0, spherical top hat, independent full σ8 normalization",
        "integration_domain_k_h_per_Mpc": "(0,infinity)",
        "finite_band_source_sha256": sha256(source.read_bytes()).hexdigest(),
        "finite_band": data["k_band"],
        "mesh": data["mesh"],
        "raw_omits_factor": "1/(2*pi^2)",
        "envelopes": {
            "low": "|T(k)| <= 3",
            "high": "|T(k)| <= 1e-9/k",
            "bound": "9*L^3+18*10^(-18)/(Rmin^4*K^2)",
        },
        "anchor": data["anchor"] | {"raw_tail_budget": str(budget(Q(8)))},
        "propagation_kernel_replayed": True,
        "radius_kernel_replayed": True,
        "row_status": "kernel_replayed true requires a closed real theorem. Other rows use the proved evaluator with runtime finite-band sums; each tail, normalization and square-root implication has kernel replay.",
        "rows": rows[:-1],
        "host_progenitor": rows[-1],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sweep", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--bits", type=int, default=80)
    parser.add_argument(
        "--certificate", nargs=2, action="append", default=[], metavar=("MASS", "MODULE_DIRECTORY")
    )
    parser.add_argument(
        "--mass-modules", type=Path, help="directory of checked EhMassN numerator/radius modules"
    )
    parser.add_argument(
        "--anchor-modules",
        type=Path,
        help="directory of the shared checked WindowCertificate module",
    )
    args = parser.parse_args()
    certs = {Q(m): Path(p) for m, p in args.certificate}
    result = evaluate(
        args.sweep, args.output, args.bits, certs, args.mass_modules, args.anchor_modules
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
