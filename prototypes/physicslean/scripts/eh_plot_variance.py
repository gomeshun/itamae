"""Integrate the verified-model EH amplitude numerically in Lean binary64.

The finite-band integral, a widened-band rerun and nested Simpson refinements
are independent numerical evidence. The reported binary64 values are not
continuous-integral enclosures. k is in h/Mpc; R is in Mpc/h; M in Msun.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import struct
import subprocess
import tempfile
import time


PROJECT = Path(__file__).resolve().parents[1]
FIELDS = ("mass_msun", "radius_mpc_h", "sigma", "S", "dS_dM", "dlnS_dlnM")
COMPARE = ("sigma", "S", "dS_dM")


def _as_bits(x):
    return int.from_bytes(struct.pack(">d", x), "big")


def _from_bits(s):
    return struct.unpack(">d", int(s).to_bytes(8, "big"))[0]


def _source(masses, meshes, *, entrypoint="main"):
    args = ", ".join(f"Float.ofBits {_as_bits(m)}" for m in masses)
    commands = [
        "import SashimiEhVariance",
        "import SashimiEisensteinHu",
        f"def {entrypoint} : IO Unit := do",
        f"  let masses : Array Float := #[{args}]",
        "  let c := SashimiEisensteinHu.prepareFloat 0.674 0.315 0.049 2.7255",
        "  let transfer := SashimiEisensteinHu.transferFloat c",
    ]
    for a, b, n in meshes:
        commands.append(f"  SashimiEhVariance.runMesh transfer masses ({a}) ({b}) {n}")
    return "\n".join(commands) + "\n"


def _read_output(output, meshes, masses):
    runs = []
    for line in output.splitlines():
        fields = line.split()
        if fields[:1] == ["anchor"] and len(fields) == 8:
            mesh = tuple(map(int, fields[1:4]))
            if len(runs) >= len(meshes) or mesh != meshes[len(runs)]:
                raise ValueError("unexpected mesh in Lean response")
            anchor, amplitude, q, zf = map(_from_bits, fields[4:])
            runs.append(
                dict(mesh=mesh, anchor_raw=anchor, amplitude=amplitude, q=q, zf=zf, rows=[])
            )
        elif fields[:1] == ["row"] and len(fields) == 7 and runs:
            vals = list(map(_from_bits, fields[1:]))
            if not all(map(math.isfinite, vals)):
                raise RuntimeError("non-finite EH result")
            runs[-1]["rows"].append(dict(zip(FIELDS, vals)))
        else:
            raise RuntimeError(f"unexpected Lean output {line[:200]}")
    if len(runs) != len(meshes):
        raise RuntimeError("incomplete mesh evaluation")
    for r in runs:
        if [row["mass_msun"] for row in r["rows"][:-1]] != masses:
            raise RuntimeError("missing or reordered mass")
        if not math.isclose(r["rows"][-1]["mass_msun"] * r["q"], 1.0e12, rel_tol=1e-15):
            raise RuntimeError("inconsistent progenitor mass")
        if not 0 < r["anchor_raw"] or not 0 < r["amplitude"]:
            raise RuntimeError("invalid normalization")
    return runs


def _differences(other, final):
    differences = {}
    per_mass = []
    for old, new in zip(other["rows"], final["rows"]):
        per_mass.append({key: abs(old[key] - new[key]) / abs(new[key]) for key in COMPARE})
    for key in COMPARE:
        i = max(range(len(per_mass)), key=lambda j: per_mass[j][key])
        differences[key] = dict(
            max_relative_change=per_mass[i][key], mass_msun=final["rows"][i]["mass_msun"]
        )
    return dict(
        anchor_raw_relative_change=abs(other["anchor_raw"] / final["anchor_raw"] - 1),
        extrema=differences,
        per_mass_relative_change=per_mass,
    )


def evaluate(mass_labels, *, ppd=8192):
    masses = [float(Fraction(str(m))) for m in mass_labels]
    if not masses or not all(math.isfinite(x) and x > 0 for x in masses):
        raise ValueError("expected positive finite mass grid")
    if type(ppd) is not int or ppd < 16 or ppd % 16:
        raise ValueError("panels per decade must be a positive multiple of 16")
    meshes = [(-8, 10, ppd // 16), (-8, 10, ppd // 4), (-8, 10, ppd), (-10, 12, ppd)]
    start = time.monotonic()
    with tempfile.TemporaryDirectory() as tmp:
        request = Path(tmp) / "EhRun.lean"
        request.write_text(_source(masses, meshes), encoding="utf-8")
        command = ["lake", "env", "lean", "--run", str(request)]
        run = subprocess.run(command, cwd=PROJECT, capture_output=True, text=True)
        if run.returncode:
            raise RuntimeError(run.stdout + run.stderr)
    elapsed = time.monotonic() - start
    coarse, middle, fine, extended = _read_output(run.stdout, meshes, masses)
    rows = [dict(row, mass_msun=str(m)) for row, m in zip(fine["rows"][:-1], mass_labels)]
    progenitor = dict(
        fine["rows"][-1],
        mass_msun=repr(fine["rows"][-1]["mass_msun"]),
        host_mass_msun=1.0e12,
        q=fine["q"],
        zf=fine["zf"],
    )
    return dict(
        model="Eisenstein & Hu 1998 with baryons, spherical top hat, z=0",
        cosmology=dict(
            h=0.674,
            omega_m=0.315,
            omega_b=0.049,
            ns=0.965,
            sigma8=0.811,
            Tcmb_K=2.7255,
            rho_crit0_coefficient_Msun_h2_Mpc3=277536627245.708,
        ),
        normalization=dict(
            R8_Mpc_h=8.0,
            sigma8=0.811,
            raw_R8_variance=fine["anchor_raw"],
            amplitude=fine["amplitude"],
            normalized_R8_sigma=math.sqrt(fine["amplitude"] * fine["anchor_raw"]),
        ),
        mass_rows=rows,
        host_progenitor=progenitor,
        quadrature=dict(
            method="composite Simpson in ln k with compensated sums",
            units=dict(k="h/Mpc", R="Mpc/h", M="Msun", dS_dM="1/Msun"),
            status="native Lean Float; finite-band numerical convergence",
            panels_per_decade=ppd,
            total_panels=18 * ppd,
            k_band_h_per_Mpc=[1e-8, 1e10],
            extended_band_h_per_Mpc=[1e-10, 1e12],
            coarse_to_mid=_differences(coarse, middle),
            mid_to_fine=_differences(middle, fine),
            finite_to_extended=_differences(fine, extended),
            runtime_seconds=elapsed,
        ),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--grid", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--panels-per-decade", type=int, default=8192)
    a = p.parse_args()
    grid = json.loads(a.grid.read_text())
    result = evaluate(grid["masses_msun"], ppd=a.panels_per_decade)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
