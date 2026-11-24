"""Certify flat matter+Λ Heath growth, its z derivative, and the EPS barrier.

Generate rational endpoint proposals, replay one Lean kernel theorem per row,
and write the JSON only after every check succeeds. The .lean witness sits
next to the JSON and can also be checked with "lake env lean <file>".
"""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess

PROJECT = Path(__file__).resolve().parents[1]


def enclose(lo, hi, bits):
    grid = 2**bits
    return {"lo": str(Q((lo * grid).__floor__(), grid)), "hi": str(Q((hi * grid).__ceil__(), grid))}


def factor(z, degree, matter):
    a = Q(1) / (1 + z)
    density = matter + (1 - matter) * a**3
    y = (1 - matter) * a**3 / density
    term = total = Q(1)
    for n in range(degree):
        term *= Q(6 * n + 9, 6 * n + 11) * y
        total += term
    remainder = Q(6 * degree + 9, 5) * term * y
    if remainder >= 1:
        raise ValueError(
            "Series residual exceeds the positive upper-enclosure threshold; increase degree"
        )
    lower = Q(2, 5) * a / density * total
    return lower, lower / (1 - remainder), density, a


def candidates(matter, redshifts, degree, bits):
    ref_lo, ref_hi, _, _ = factor(Q(0), degree, matter)
    for z in redshifts:
        lo, hi, density, a = factor(z, degree, matter)
        dl, du = (lo / ref_hi, hi / ref_lo) if z else (Q(1), Q(1))
        alpha, beta = Q(3, 2) * matter * a / density, a * a / density
        yield {
            "z": str(z),
            "D": enclose(dl, du, bits),
            "dD_dz": enclose(alpha * dl - beta / ref_lo, alpha * du - beta / ref_hi, bits),
            "delta_c": enclose(Q(1686, 1000) / du, Q(1686, 1000) / dl, bits),
        }


def render(rows, degree, matter):
    lines = [
        "import SashimiHeath",
        "open SashimiHeath",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 0",
    ]
    for i, row in enumerate(rows):
        boxes = (
            f"⟨({row[name]['lo']}), ({row[name]['hi']})⟩" for name in ("D", "dD_dz", "delta_c")
        )
        lines.append(f"def heathRow{i} : Row := ⟨{', '.join(boxes)}⟩")
        lines.append(
            f"theorem heathCheck{i} : check {degree} ({matter}) ({row['z']}) heathRow{i} = true := by decide +kernel"
        )
        m, z = f"(({matter} : ℚ) : ℝ)", f"(({row['z']} : ℚ) : ℝ)"
        lines.append(
            f"theorem heathContains{i} : heathRow{i}.D.Contains (D {m} {z}) ∧ "
            f"heathRow{i}.dD_dz.Contains (deriv (D {m}) {z}) ∧ "
            f"heathRow{i}.delta_c.Contains (barrier {m} {z}) := check_sound heathCheck{i}"
        )
    lines.extend(("#print axioms check_sound", f"#print axioms heathContains{len(rows) - 1}"))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--degree", type=int, default=90)
    parser.add_argument("--bits", type=int, default=72)
    parser.add_argument("--matter", type=Q, default=Q(63, 200))
    parser.add_argument("--z", type=Q, nargs="+", default=[Q(i, 4) for i in range(29)])
    args = parser.parse_args()
    if not 0 < args.matter <= 1 or min(args.z) < 0 or min(args.degree, args.bits) < 0:
        parser.error("Require 0 < Ωm ≤ 1, z ≥ 0, nonnegative degree and bits")
    rows = list(candidates(args.matter, args.z, args.degree, args.bits))
    witness = args.output.with_suffix(".lean")
    if witness == args.output:
        parser.error("Output path must differ from the Lean witness path")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    witness.write_text(render(rows, args.degree, args.matter))
    subprocess.run(["lake", "env", "lean", "-j1", str(witness.resolve())], cwd=PROJECT, check=True)
    args.output.write_text(
        json.dumps(
            {
                "parameters": {
                    "omega_m": str(args.matter),
                    "omega_lambda": str(1 - args.matter),
                    "radiation": "0",
                    "collapse_threshold_today": "843/500",
                },
                "definition": "a=1/(1+z), R(t)=m+(1-m)t^6, G(a)=sqrt(R(sqrt(a)))/a^(3/2)*integral_0^sqrt(a) 2t^4/R(t)^(3/2) dt, D=G(a)/G(1)",
                "derivative": "dD/dz=(3/2)*m*a/[m+(1-m)a^3]*D - a^2/([m+(1-m)a^3]*G(1))",
                "barrier": "1.686/D",
                "degree": args.degree,
                "outward_rounding_bits": args.bits,
                "kernel_replayed": True,
                "theorem": "SashimiHeath.check_sound",
                "witness": witness.name,
                "rows": rows,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
