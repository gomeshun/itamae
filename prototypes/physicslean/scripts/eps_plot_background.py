"""Enclose the SASHIMI-C cold-EPS background on a redshift grid with Lean.

The public C model sets zi=0, disables prompt cusps and uses the Ludlow CDM
variance. dM_dz is its negative fitted accretion-rate expression, including
the Hu–Kravtsov virial mass conversion used by Na_calc.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic

if __package__:
    from .sashimi import rational
    from .table_replay import check_proof, literal
else:
    from sashimi import rational
    from table_replay import check_proof, literal

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_GRID = tuple(Fraction(i, 4) for i in range(25))
UNITS = {
    "H": "km/s/Mpc",
    "D": "1",
    "dD_dz": "1/redshift",
    "delta_c": "1",
    "M200": "Msun",
    "dM_dz": "Msun/redshift",
}


def _source(m, h, mass, bits, degree):
    return f"""import SashimiHalo
import SashimiBackground
open FloatLib.Numerics

private def m : ℚ := {literal(m)}
private def h : ℚ := {literal(h)}
private def m0 : ℚ := {literal(mass)}
private def value (key : String) (z : ℚ) : Option RationalInterval :=
  let growth := fun q => SashimiGrowth.enclose q m (1-m) h z 1 {bits} {degree}
  let host := fun q => SashimiHalo.evaluate? q [m,h,m0,z,0] {bits} m0 {degree}
  match key with
  | "H" => SashimiBackground.enclose .hubble (SashimiBackground.flat m) (100*h) 4.30091e-9 z 1 {bits} {degree}
  | "D" => growth .growth
  | "dD_dz" => growth .growthDerivative
  | "delta_c" => growth .collapseThreshold
  | "M200" => host "cold-m200"
  | "dM_dz" => host "cold-virial-rate"
  | _ => none
"""


def _run(source, zs):
    keys = ", ".join(json.dumps(key) for key in UNITS)
    grid = ", ".join(map(literal, zs))
    program = f"""
def main : IO Unit := do
  for (z, i) in ([{grid}] : List ℚ).zipIdx do
    for key in [{keys}] do
      let some v := value key z | throw (IO.userError s!"unresolved {{key}} at {{z}}")
      IO.println s!"EPS {{i}} {{key}} {{v.lo}} {{v.hi}}"
"""
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "Background.lean"
        path.write_text(source + program, encoding="utf-8")
        command = subprocess.run(
            ["lake", "env", "lean", "--run", str(path)],
            cwd=PROJECT,
            capture_output=True,
            text=True,
        )
    if command.returncode:
        raise RuntimeError(command.stdout + command.stderr)
    rows = [{"z": str(z)} for z in zs]
    for line in command.stdout.splitlines():
        if line.startswith("EPS "):
            _, index, key, lo_s, hi_s = line.split()
            lo, hi = Fraction(lo_s), Fraction(hi_s)
            if not (0 <= hi - lo <= max(abs(lo), abs(hi)) * Fraction("1e-10")):
                raise ArithmeticError(
                    f"Relative precision target missed at z={zs[int(index)]}: {key}"
                )
            row = rows[int(index)]
            if key not in UNITS or key in row:
                raise RuntimeError(f"Unexpected Lean output: {line}")
            row[key] = {"lo": str(lo), "hi": str(hi)}
    if any(row.keys() != {"z", *UNITS} for row in rows):
        raise RuntimeError("Incomplete Lean response")
    return rows


def _replay(source, rows):
    for start in range(0, len(rows), 5):
        proofs = []
        goals = []
        names = []
        for i, row in enumerate(rows[start : start + 5]):
            name = f"row{i}"
            bounds = _row_bounds(row)
            goals.append(f"{name}Goal")
            names.append(name)
            proofs.append(
                f"private def {name}Goal : Prop := {bounds}\n"
                f"private theorem {name} : {name}Goal := by unfold {name}Goal; decide +kernel"
            )
        all_proofs = "⟨" + ", ".join(names) + "⟩" if len(names) > 1 else names[0]
        check_proof(
            source
            + "\nset_option maxRecDepth 200000\nset_option maxHeartbeats 12000000\n"
            + "\n".join(proofs)
            + f"\ntheorem concrete : {' ∧ '.join(goals)} := {all_proofs}\n"
            + "#print axioms concrete\n",
            PROJECT,
            max_memory_mb=6000,
            timeout_seconds=300,
        )


def _row_bounds(row):
    bounds = []
    for key in UNITS:
        r = row[key]
        expr = f'value "{key}" {literal(Fraction(row["z"]))}'
        lo, hi = (literal(Fraction(r[k])) for k in ("lo", "hi"))
        bounds.append(f"({expr} = some ⟨{lo}, {hi}⟩)")
    return " ∧ ".join(bounds)


def evaluate(
    redshifts=DEFAULT_GRID,
    *,
    omega_m="0.315",
    h="0.674",
    host_mass=10**12,
    bits=96,
    degree=32,
    replay=True,
):
    """Return exact rational endpoint strings for each physical input and observable."""
    zs = [rational(z) for z in redshifts]
    m, h, mass = map(rational, (omega_m, h, host_mass))
    if not zs or any(z < 0 for z in zs) or not 0 < m <= 1 or min(h, mass) <= 0:
        raise ValueError(
            "Expected nonnegative redshifts, positive host mass, h and 0 < omega_m <= 1"
        )
    if any(type(v) is not int or v < 1 for v in (bits, degree)):
        raise ValueError("bits and degree must be positive integers")
    if type(replay) is not bool:
        raise ValueError("replay must be a boolean")
    source = _source(m, h, mass, bits, degree)
    start = monotonic()
    rows = _run(source, zs)
    computed = monotonic()
    if replay:
        _replay(source, rows)
    return {
        "model": "SASHIMI-C cold Ludlow CDM, normalized CPT, fitted virial accretion rate",
        "parameters": {"omega_m": str(m), "h": str(h), "host_mass_Msun": str(mass), "zi": "0"},
        "units": UNITS,
        "bits": bits,
        "degree": degree,
        "relative_width_budget": "1/10000000000",
        "kernel_replayed": replay,
        "runtime_seconds": computed - start,
        "replay_seconds": monotonic() - computed if replay else None,
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--grid", type=Path, help="JSON grid with exact redshifts, default 0:0.25:6"
    )
    parser.add_argument("--no-replay", action="store_true")
    args = parser.parse_args()
    grid = DEFAULT_GRID
    if args.grid is not None:
        grid = json.loads(args.grid.read_text(encoding="utf-8"))["redshifts"]
    result = json.dumps(evaluate(grid, replay=not args.no_replay), indent=2) + "\n"
    if args.output is None:
        print(result, end="")
    else:
        args.output.write_text(result, encoding="utf-8")


if __name__ == "__main__":
    main()
