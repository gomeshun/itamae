"""Generate and kernel-check the fixed two-epoch EPS conditional density.

From prototypes/physicslean after lake build, with all three files on one grid:
  python3 scripts/eps_plot_kernel.py BACKGROUND.json VARIANCE.json KERNEL.json

The adjacent .lean certificate can be replayed with lake env lean CERT.lean.
Every input interval used in a density is reproduced from the CPT barrier or
the fiducial SASHIMI-C Ludlow variance in that certificate.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic

if __package__:
    from .table_replay import check_proof, literal
else:
    from table_replay import check_proof, literal

PROJECT = Path(__file__).resolve().parents[1]
M0 = 10**12


def box(endpoints):
    lo, hi = (Fraction(endpoints[k]) for k in ("lo", "hi"))
    if lo > hi:
        raise ValueError("reversed input interval")
    return f"⟨{literal(lo)}, {literal(hi)}⟩"


def inputs(background, variance):
    parameters = background["parameters"]
    if any(
        Fraction(parameters[k]) != v
        for k, v in (
            ("omega_m", Fraction(".315")),
            ("h", Fraction(".674")),
            ("host_mass_Msun", M0),
        )
    ):
        raise ValueError("requires the fiducial SASHIMI-C cosmology and 1e12 Msun parent")
    if any(
        Fraction(variance[k]) != v
        for k, v in (
            ("omega_m", Fraction(".315")),
            ("omega_lambda", Fraction(".685")),
            ("h", Fraction(".674")),
            ("pivot_msun", 10**10 / Fraction(".674")),
        )
    ):
        raise ValueError("requires the SASHIMI-C Ludlow fit at the fiducial cosmology")
    barriers = {}
    for r in background["rows"]:
        z = Fraction(r["z"])
        if z in (0, 1):
            if z in barriers:
                raise ValueError(f"duplicate barrier at redshift {z}")
            barriers[z] = r["delta_c"]
    mass_boxes = {}
    for r in variance["rows"]:
        if Fraction(r["z"]) == 0:
            m = Fraction(r["mass_msun"])
            if not 0 < m <= M0 or m in mass_boxes:
                raise ValueError("expected unique positive reference masses ≤ parent mass")
            mass_boxes[m] = r["S"]
    if len(barriers) != 2 or M0 not in mass_boxes or len(mass_boxes) < 2:
        raise ValueError("missing z=0 or z=1 barrier, parent or child mass variance")
    return barriers, mass_boxes


def generate(background, variance, certificate, *, bits=96, degree=40):
    barriers, mass_boxes = inputs(background, variance)
    precisions = (
        background["bits"],
        background["degree"],
        variance["bits"],
        variance["degree"],
        bits,
        degree,
    )
    if any(type(x) is not int or x <= 0 for x in precisions):
        raise ValueError("precision and degree must be positive integers")
    bg_bits, bg_degree, var_bits, var_degree, _, _ = precisions
    masses = sorted(m for m in mass_boxes if m < M0)
    definitions = [
        "import SashimiEpsLudlow",
        "import SashimiEpsKernel",
        "open FloatLib.Numerics",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 8000000",
        "set_option Elab.async false",
        f"private def parentBarrier : RationalInterval := {box(barriers[0])}",
        f"private def childBarrier : RationalInterval := {box(barriers[1])}",
        f"private def parent : RationalInterval := {box(mass_boxes[M0])}",
        f"private def compute (s : RationalInterval) := EpsKernel.enclose? childBarrier parentBarrier s parent {bits} {degree}",
    ]
    for i, m in enumerate(masses):
        definitions.append(f"private def child{i} : RationalInterval := {box(mass_boxes[m])}")
    prefix = "\n".join(definitions)
    cases = ", ".join(f"child{i}" for i in range(len(masses)))
    runner = f"""
def main : IO Unit := do
  for (s, i) in ([{cases}] : List RationalInterval).zipIdx do
    let some I := compute s | throw (IO.userError s!"unresolved EPS density at index {{i}}")
    IO.println s!"EPS {{i}} {{I.lo}} {{I.hi}}"
"""
    start = monotonic()
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "Eval.lean"
        src.write_text(prefix + runner)
        run = subprocess.run(
            ["lake", "env", "lean", "--run", str(src)], cwd=PROJECT, text=True, capture_output=True
        )
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    rows = []
    for line in run.stdout.splitlines():
        if line.startswith("EPS "):
            _, i, lo, hi = line.split()
            if int(i) != len(rows) or not 0 < Fraction(lo) <= Fraction(hi):
                raise ValueError(f"invalid output {line}")
            rows.append(dict(mass_msun=str(masses[len(rows)]), F=dict(lo=lo, hi=hi)))
    if len(rows) != len(masses):
        raise ValueError("incomplete kernel output")
    runtime = monotonic() - start

    # Prove each input box describes this particular growth and variance model,
    # then compose the guarded interval computation to obtain closed real claims.
    checks = [
        prefix,
        f"""
private noncomputable def barrier (z : ℝ) := SashimiGrowth.value .collapseThreshold 0.315 0.685 0.674 z
private noncomputable def variance0 (m : ℝ) := (EpsLudlow.sigmaAt 0 m)^2
private theorem variance_sound (m : ℚ) (B : RationalInterval)
    (h : (EpsLudlow.evaluate? m 0 {var_bits} {var_degree}).map (·.variance) = some B) :
    B.Contains (variance0 m) := by
  cases hr : EpsLudlow.evaluate? m 0 {var_bits} {var_degree} with
  | none => simp [hr] at h
  | some r =>
    simp only [hr, Option.map_some, Option.some.injEq] at h
    subst B
    exact (EpsLudlow.sound hr).2.1
private theorem parent_checked : parent.Contains (variance0 {M0}) :=
  variance_sound {M0} parent (by decide +kernel)
""",
    ]
    for z, name in ((0, "parentBarrier"), (1, "childBarrier")):
        checks.append(f"""
private theorem {name}_checked : {name}.Contains (barrier {z}) := by
  apply (SashimiGrowth.enclose_sound (tol := 1) (bits := {bg_bits}) (degree := {bg_degree})
    (z := {z}) (q := .collapseThreshold) (m := 0.315) (l := 0.685) (h := 0.674)
    (I := {name}) (by decide +kernel)).1
""")
    claims = []
    for i, (m, row) in enumerate(zip(masses, rows, strict=True)):
        claim = (
            f"out{i}.Contains (EpsKernel.density (barrier 1 - barrier 0) "
            f"(variance0 {literal(m)} - variance0 {M0}))"
        )
        checks.extend(
            [
                f"private def out{i} : RationalInterval := {box(row['F'])}",
                f"private theorem c{i} : child{i}.Contains (variance0 {literal(m)}) := "
                f"variance_sound {literal(m)} child{i} (by decide +kernel)",
                f"private theorem f{i} : {claim} :=\n"
                f"  (EpsKernel.sound childBarrier_checked parentBarrier_checked c{i} parent_checked "
                f"(bits := {bits}) (degree := {degree}) (I := out{i}) (by decide +kernel)).2.2",
            ]
        )
        claims.append(f"({claim})")
    term = "f0" if len(rows) == 1 else "⟨" + ", ".join(f"f{i}" for i in range(len(rows))) + "⟩"
    checks.append(
        "theorem concrete : " + " ∧\n  ".join(claims) + " := " + term + "\n#print axioms concrete\n"
    )
    proof = "\n".join(checks)
    start = monotonic()
    check_proof(proof, PROJECT, max_memory_mb=4096, timeout_seconds=600)
    certificate.write_text(proof, encoding="utf-8")
    return {
        "model": "Yang et al. Eq. (14), fixed two-epoch EPS conditional mass-fraction density per dS",
        "parameters": {
            "parent_mass_msun": str(M0),
            "parent_z": "0",
            "child_z": "1",
            "variance_z": "0",
        },
        "kernel_replayed": True,
        "certificate": certificate.name,
        "bits": bits,
        "degree": degree,
        "evaluation_seconds": runtime,
        "replay_seconds": monotonic() - start,
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("background", "variance", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--certificate", type=Path)
    parser.add_argument("--bits", type=int, default=96)
    parser.add_argument("--degree", type=int, default=40)
    args = parser.parse_args()
    result = generate(
        json.loads(args.background.read_text()),
        json.loads(args.variance.read_text()),
        args.certificate or args.output.with_suffix(".lean"),
        bits=args.bits,
        degree=args.degree,
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
