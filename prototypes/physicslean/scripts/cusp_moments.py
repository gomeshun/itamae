"""Check finite-band SASHIMI-C BBKS moments for a supplied z=31 CSV (k,Delta²)."""

import argparse
import hashlib
import json
import subprocess
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory

PROJECT = Path(__file__).resolve().parents[1]
QUANTITIES = {
    **{f"variance{j}": f"(.variance {j})" for j in range(3)},
    **{f"sigma{j}": f"(.sigma {j})" for j in range(3)},
    "radius": ".radius",
    "gamma": ".gamma",
}


def rat(x):
    return f"({x.numerator}/{x.denominator})"


def read_rows(path):
    data = path.read_bytes()
    rows = []
    for raw in data.decode().splitlines():
        line = raw.partition("#")[0].strip()
        if line:
            k, delta = map(Fraction, line.split(","))
            if k <= 0 or delta <= 0 or (rows and k <= rows[-1][0]):
                raise ValueError("k and Delta² must be positive and k strictly increasing")
            rows.append((k, delta))
    if len(rows) < 2:
        raise ValueError("at least two rows required")
    return rows, hashlib.sha256(data).hexdigest()


def source(args, rows):
    pairs = ",\n  ".join(f"⟨{rat(k)}, {rat(y)}⟩" for k, y in rows)
    prefix = f"""import SashimiCusp
open SashimiCusp SashimiTable FloatLib.Numerics FloatLib.Numerics.Interval
set_option maxRecDepth 32768
set_option maxHeartbeats 0
namespace CuspInput
def rows : List Knot := [
  {pairs}]
def policy := EdgePolicy.{args.edges}
def low : ℚ := {rat(args.k_min)}
def high : ℚ := {rat(args.k_max)}
def kfs : ℚ := {rat(args.kfs)}
def z : ℚ := {rat(args.redshift)}
def tol : ℚ := {rat(args.tolerance)}
def boxes := moments? policy rows low high kfs {args.depth} {args.bits} {args.degree}
"""
    probe = [
        "#eval do",
        '  let some ms := boxes | IO.println "unresolved"; return',
        '  for m in ms do IO.println s!"box {m.lo} {m.hi}"',
    ]
    for key in args.quantities:
        probe += [
            f"  match peak? {QUANTITIES[key]} z tol ms {args.bits} {args.degree} with",
            f'  | some box => IO.println s!"{key} {{box.lo}} {{box.hi}}"',
            f'  | none => IO.println "{key} unresolved"',
        ]
    return prefix, "\n".join(probe) + "\nend CuspInput\n"


def certificate(prefix, args, boxes, outputs):
    def interval(pair):
        return "⟨" + ", ".join(map(rat, pair)) + "⟩"

    ms = "[" + ", ".join(map(interval, boxes)) + "]"
    lines = [
        prefix,
        f"def stored : List RationalInterval := {ms}",
        "theorem checkedMoments : boxes = some stored := by decide +kernel",
    ]
    for key, pair in outputs.items():
        if pair is None:
            continue
        q = QUANTITIES[key]
        lines += [
            f"def {key}Value : RationalInterval := {interval(pair)}",
            f"theorem {key}_checked : {key}Value.Contains",
            f"    (statistic {q} z (moment policy rows low high kfs)) ∧ {key}Value.hi - {key}Value.lo ≤ tol :=",
            f"  peak_sound checkedMoments (bits' := {args.bits}) (degree' := {args.degree}) (by decide +kernel)",
            f"#print axioms {key}_checked",
        ]
    return "\n".join(lines) + "\nend CuspInput\n"


def simpson_comparison(args, rows):
    import numpy as np
    from scipy.integrate import simpson
    from scipy.interpolate import interp1d

    k, delta = np.array(rows, dtype=float).T
    grid = np.geomspace(float(args.k_min), float(args.k_max), num=1000)
    interp = interp1d(np.log(k), np.log(delta), fill_value="extrapolate")
    y = np.exp(interp(np.log(grid))) * (32 / (1 + float(args.redshift))) ** 1.802
    y *= np.exp(-((grid / float(args.kfs)) ** 2))
    variances = [simpson(y * grid ** (2 * j), x=np.log(grid)) for j in range(3)]
    s = np.sqrt(variances)
    vals = {
        **{f"variance{i}": v for i, v in enumerate(variances)},
        **{f"sigma{i}": v for i, v in enumerate(s)},
        "radius": np.sqrt(3) * s[1] / s[2],
        "gamma": s[1] ** 2 / (s[0] * s[2]),
    }
    return {key: float(vals[key]) for key in args.quantities}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table", type=Path)
    ap.add_argument("--k-min", type=Fraction, default=Fraction("0.1"))
    ap.add_argument("--k-max", type=Fraction, default=Fraction("1e7"))
    ap.add_argument("--kfs", type=Fraction, default=Fraction("1.06e6"))
    ap.add_argument("--redshift", type=Fraction, default=Fraction(0))
    ap.add_argument("--edges", choices=["inside", "extrapolate"], default="extrapolate")
    ap.add_argument("--depth", type=int, default=2)
    ap.add_argument("--bits", type=int, default=64)
    ap.add_argument("--degree", type=int, default=20)
    ap.add_argument(
        "--tolerance",
        type=Fraction,
        required=True,
        help="absolute output width budget in each selected quantity's unit",
    )
    ap.add_argument(
        "--quantities",
        nargs="+",
        choices=QUANTITIES,
        default=["sigma0", "sigma1", "sigma2", "radius", "gamma"],
    )
    ap.add_argument(
        "--simpson",
        action="store_true",
        help="also compute upstream-style floating 1000-point Simpson reference",
    )
    ap.add_argument(
        "--certificate", type=Path, help="save the standalone Lean file after kernel checking"
    )
    args = ap.parse_args()
    try:
        rows, checksum = read_rows(args.table)
        if not (
            0 < args.k_min < args.k_max
            and args.kfs > 0
            and args.redshift > -1
            and min(args.bits, args.depth, args.degree, args.tolerance) >= 0
        ):
            raise ValueError("invalid physical domain or precision setting")
        if args.edges == "inside" and not (rows[0][0] <= args.k_min < args.k_max <= rows[-1][0]):
            raise ValueError("inside policy requires full data coverage")
        if len(set(args.quantities)) != len(args.quantities):
            raise ValueError("duplicate quantities")
    except (ValueError, ZeroDivisionError) as err:
        ap.error(str(err))

    base, probe = source(args, rows)
    info = dict(
        table_sha256=checksum,
        rows=len(rows),
        redshift=str(args.redshift),
        kfs_Mpc_inverse=str(args.kfs),
        k_band_Mpc_inverse=[str(args.k_min), str(args.k_max)],
        table_band_Mpc_inverse=[str(rows[0][0]), str(rows[-1][0])],
        edges=args.edges,
        depth=args.depth,
        bits=args.bits,
        degree=args.degree,
        absolute_width_budget=str(args.tolerance),
    )
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "CuspInput.lean"
        path.write_text(base + probe, encoding="utf8")

        def run():
            checked = subprocess.run(
                ["lake", "env", "lean", str(path)], cwd=PROJECT, capture_output=True, text=True
            )
            if checked.returncode:
                raise SystemExit(checked.stdout + checked.stderr)
            return checked.stdout

        output = run().splitlines()
        boxes = []
        outputs = {}
        for line in output:
            words = line.split()
            if words and words[0] == "box":
                boxes.append(tuple(map(Fraction, words[1:])))
            elif words and words[0] in QUANTITIES:
                outputs[words[0]] = (
                    None if words[1] == "unresolved" else tuple(map(Fraction, words[1:]))
                )
        if len(boxes) != 3 or any(key not in outputs for key in args.quantities):
            if "unresolved" not in output:
                raise RuntimeError(f"invalid Lean output {output!r}")
            print(json.dumps(dict(info, status="unresolved")))
            return 2
        checked = certificate(base, args, boxes, outputs)
        path.write_text(checked, encoding="utf8")
        run()
        if args.certificate:
            args.certificate.write_text(checked, encoding="utf8")
    info["enclosures"] = {
        k: None if v is None else dict(lower=str(v[0]), upper=str(v[1])) for k, v in outputs.items()
    }
    info["status"] = "enclosed" if all(v is not None for v in outputs.values()) else "unresolved"
    if args.simpson:
        info["floating_simpson_1000"] = simpson_comparison(args, rows)
    print(json.dumps(info, indent=2))
    return 0 if info["status"] == "enclosed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
