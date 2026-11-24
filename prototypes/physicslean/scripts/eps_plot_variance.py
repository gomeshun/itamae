"""Batch SASHIMI-C's default CDM Ludlow fit through the verified Lean evaluator.

Input JSON is a list, or an object with a "rows" list, of {"mass_msun": "1e9",
"z": "0"} records. A numeric JSON literal is read as its exact written decimal.
For the shared plot grid {"masses_msun":[...], "redshifts":[...]}, evaluate S(M)
at redshift zero as in the native SASHIMI-C s_func and Na_calc.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import tempfile
import time

if __package__:
    from .spectrum_table import PROJECT, number
    from .table_replay import check_proof, literal
else:
    from spectrum_table import PROJECT, number
    from table_replay import check_proof, literal


QUANTITIES = ("sigma", "S", "dS_dM")
PIN = "9334d14622a0d66082375ea9b66ce21a6d8f1ad2"

# Run one Lean process for all measurements. The distinct replay process
# certifies literal inputs and output rationals, independent of line parsing.
RUNTIME = """
private def epsRun (bits degree : Nat) : IO Unit := do
  let input ← IO.getStdin
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    let [m,z] := line.trimAscii.toString.splitOn " " | throw (IO.userError "expected M z")
    let .ok m := SashimiTableRequest.rational m | throw (IO.userError "invalid mass")
    let .ok z := SashimiTableRequest.rational z | throw (IO.userError "invalid redshift")
    let some r := EpsLudlow.evaluate? m z bits degree | throw (IO.userError s!"unresolved at {m} {z}")
    IO.println s!"{m} {z} {r.sigma.lo} {r.sigma.hi} {r.variance.lo} {r.variance.hi} {r.derivative.lo} {r.derivative.hi}"
"""


def _eval(requests, bits, degree):
    with tempfile.TemporaryDirectory() as tempdir:
        path = Path(tempdir) / "Run.lean"
        path.write_text(
            "import SashimiTableRequest\nimport SashimiEpsLudlow\n"
            + RUNTIME
            + f"\ndef main : IO Unit := epsRun {bits} {degree}\n"
        )
        run = subprocess.run(
            ["lake", "env", "lean", "--run", str(path)],
            cwd=PROJECT,
            input="".join(f"{m} {z}\n" for m, z in requests),
            text=True,
            capture_output=True,
        )
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    out = {}
    for line in run.stdout.splitlines():
        m, z, *values = map(Fraction, line.split())
        if len(values) != 6 or any(lo > hi for lo, hi in zip(values[::2], values[1::2])):
            raise RuntimeError(f"invalid Lean response {line}")
        out[m, z] = dict(
            zip(
                QUANTITIES,
                [dict(lo=str(lo), hi=str(hi)) for lo, hi in zip(values[::2], values[1::2])],
            )
        )
    if list(out) != requests:
        raise RuntimeError("returned Lean inputs differ from the request")
    return out


def _replay(batch, values, bits, degree):
    source = [
        "import SashimiEpsLudlow",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 8000000",
        "open EpsLudlow",
    ]
    claims = []
    for i, (m, z) in enumerate(batch):
        boxes = values[m, z]
        box_terms = []
        for key in QUANTITIES:
            lo, hi = (literal(Fraction(boxes[key][end])) for end in ("lo", "hi"))
            box_terms.append(f"⟨{lo}, {hi}⟩")
        claim = f"Refines {literal(m)} {literal(z)} row{i}"
        claims.append(claim)
        source.extend(
            [
                f"private def row{i} : Row := ⟨{', '.join(box_terms)}⟩",
                f"private theorem check{i} : evaluate? {literal(m)} {literal(z)} {bits} {degree} = "
                f"some row{i} := by decide +kernel",
                f"private theorem ok{i} : {claim} := sound check{i}",
            ]
        )
    if len(claims) == 1:
        prop, term = claims[0], "ok0"
    else:
        prop = " ∧ ".join(f"({claim})" for claim in claims)
        term = "⟨" + ", ".join(f"ok{i}" for i in range(len(batch))) + "⟩"
    source.extend([f"theorem concrete : {prop} := {term}", "#print axioms concrete"])
    check_proof("\n".join(source) + "\n", PROJECT, max_memory_mb=4096, timeout_seconds=600)


def evaluate_variance(rows, *, bits=96, degree=40, replay=True, batch_size=12):
    """Preserve input order and exact unit annotations; evaluate duplicate inputs once."""
    if not rows or any(type(p) is not int or p < 1 for p in (bits, degree, batch_size)):
        raise ValueError(
            "nonempty grid and positive integer precision, degree and batch size required"
        )
    inputs = [(number(row["mass_msun"]), number(row.get("z", 0))) for row in rows]
    if any(m <= 0 or z < 0 for m, z in inputs):
        raise ValueError("mass must be positive and redshift nonnegative")
    unique = list(dict.fromkeys(inputs))
    start = time.monotonic()
    values = _eval(unique, bits, degree)
    seconds = time.monotonic() - start
    start = time.monotonic()
    if replay:
        for i in range(0, len(unique), batch_size):
            _replay(unique[i : i + batch_size], values, bits, degree)
    return {
        "model": "SASHIMI-C default non-cusp CDM Ludlow (2016) fit with CPT growth",
        "source_commit": PIN,
        "h": "337/500",
        "omega_m": "63/200",
        "omega_lambda": "137/200",
        "pivot_msun": str(Fraction(10**10) / Fraction("0.674")),
        "units": {"mass_msun": "Msun", "sigma": "1", "S": "1", "dS_dM": "1/Msun"},
        "bits": bits,
        "degree": degree,
        "unique_inputs": len(unique),
        "evaluation_seconds": seconds,
        "kernel_replay": "passed" if replay else "not-requested",
        "replay_seconds": time.monotonic() - start if replay else None,
        "rows": [dict(mass_msun=str(m), z=str(z), **values[m, z]) for m, z in inputs],
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "Example: python scripts/eps_plot_variance.py --grid requests.json --output variance.json. "
            'Example input: [{"mass_msun":"1e8"},{"mass_msun":"1e11"},{"mass_msun":"1e8","z":1},'
            '{"mass_msun":"1e11","z":1}].'
        ),
    )
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bits", type=int, default=96)
    parser.add_argument("--degree", type=int, default=40)
    parser.add_argument("--batch-size", type=int, default=12)
    parser.add_argument("--no-replay", action="store_true")
    args = parser.parse_args()
    data = json.loads(args.grid.read_text(), parse_float=Decimal)
    if isinstance(data, dict) and "masses_msun" in data:
        rows = [{"mass_msun": m, "z": "0"} for m in data["masses_msun"]]
    else:
        rows = data["rows"] if isinstance(data, dict) else data
    result = evaluate_variance(
        rows,
        bits=args.bits,
        degree=args.degree,
        replay=not args.no_replay,
        batch_size=args.batch_size,
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
