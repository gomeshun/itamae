"""Select a Saikawa–Shirai 2018 source cell and verify its exact error bands in Lean.

Supply the v2 ancillary file https://arxiv.org/src/1803.01038v2/anc/standardmodel2018.dat.
The SHA-256 pins the external table. No table download is needed to check the built-in examples.
"""

import argparse
from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

TABLE_SHA256 = "a4bf5f18dacb084a00e7202a29c689b4ae385950cfeb00be69c98aeb30e8014b"
PROJECT = Path(__file__).resolve().parents[1]


def source_rows(path):
    data = path.read_bytes()
    if sha256(data).hexdigest() != TABLE_SHA256:
        raise ValueError("unexpected SHA-256 for arXiv:1803.01038v2 standardmodel2018.dat")
    return [
        tuple(map(Fraction, line.split()))
        for line in data.decode("ascii").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def source_cell(rows, temperature):
    if not rows[0][0] <= temperature <= rows[-1][0]:
        raise ValueError("temperature outside the tabulated range in GeV")
    i = min(bisect_right([row[0] for row in rows], temperature) - 1, len(rows) - 2)
    return i, rows[i], rows[i + 1]


def decimal_or_ratio(q):
    return f"({q.numerator}/{q.denominator})" if q.denominator != 1 else str(q.numerator)


def render(rows, temperature):
    i, left, right = source_cell(rows, temperature)
    w = (temperature - left[0]) / (right[0] - left[0])
    t = decimal_or_ratio(temperature)

    def knot(row):
        return "⟨" + ", ".join(map(decimal_or_ratio, row)) + "⟩"

    lines = [
        "import SashimiThermal",
        "",
        "open Cosmology.EffectiveDegreesOfFreedom",
        "open FloatLib.Numerics FloatLib.Numerics.Interval",
        "",
        f"-- arXiv:1803.01038v2 data rows {i + 1} and {i + 2} (one-based, excluding header).",
        f"-- Table SHA-256 {TABLE_SHA256}",
        f"def inputCell : Cell := ⟨{knot(left)}, {knot(right)}⟩",
        "",
    ]
    for q, field, model in [("radiation", 1, "gρ"), ("entropy", 3, "gs")]:
        center, error = [left[k] + w * (right[k] - left[k]) for k in (field, field + 1)]
        lo, hi = map(decimal_or_ratio, (center - error, center + error))
        lines += [
            f"theorem {q}_band (gρ : GRho) (gs : GS) (h : SaikawaShiraiDoF inputCell gρ gs) :",
            f"    ({lo} : ℝ) ≤ {model} {t} ∧ {model} {t} ≤ ({hi} : ℝ) := by",
            f"  have v := tableBand_sound h (T := {t}) (q := .{q})",
            f"    (I := ⟨{lo}, {hi}⟩) (by decide +kernel)",
            "  norm_num [RationalInterval.Contains] at v ⊢",
            "  exact v",
            f"#print axioms {q}_band",
            "",
        ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path, help="standardmodel2018.dat (exact v2 file)")
    parser.add_argument("temperature", type=Fraction, help="photon temperature in GeV")
    parser.add_argument("--lean", type=Path, help="also retain the input-specific Lean source")
    args = parser.parse_args()
    try:
        proof = render(source_rows(args.table), args.temperature)
    except ValueError as exc:
        parser.error(str(exc))
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "ThermalTable.lean"
        path.write_text(proof, encoding="utf-8")
        subprocess.run(["lake", "env", "lean", str(path)], cwd=PROJECT, check=True)
        if args.lean:
            args.lean.write_text(proof, encoding="utf-8")
    print(proof)


if __name__ == "__main__":
    main()
