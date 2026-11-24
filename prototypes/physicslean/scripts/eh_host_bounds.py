"""Certified matched host assembly from the exact Correa formation mass ratio.

Without --paired use each supplied full-band finite numerator and a shared
sigma8 anchor. With --paired import the independently checked paired-window
integral and prove unconditional mass and signed accretion-rate bounds.
"""

import argparse
from fractions import Fraction as Q
import json
import os
from pathlib import Path
import subprocess

if __package__:
    from .eh_full_grid import PROJECT, box, budget, interval_expr as I, literal, pair
else:
    from eh_full_grid import PROJECT, box, budget, interval_expr as I, literal, pair


def difference(child, parent, anchor):
    c, p, a = (pair(v["raw"]) for v in (child, parent, anchor))
    rc, rp = (pair(v["radius"])[0] for v in (child, parent))
    low = c[0] - p[1] - budget(rp)
    high = c[1] - p[0] + budget(rc)
    if low <= 0 or a[0] <= 0:
        raise ValueError(
            "Independent raw enclosures fail to separate; use a paired-window certificate"
        )
    return box(Q(".811") ** 2 * low / (a[1] + budget(Q(8))), Q(".811") ** 2 * high / a[0])


def check(path, text, modules=(), run=False):
    path.write_text(text)
    env = os.environ.copy()
    pkg = subprocess.check_output(
        ["lake", "env", "printenv", "LEAN_PATH"], text=True, cwd=PROJECT
    ).strip()
    env["LEAN_PATH"] = os.pathsep.join([*(str(p.resolve()) for p in modules), pkg])
    lean = subprocess.check_output(["lake", "env", "which", "lean"], text=True, cwd=PROJECT).strip()
    cmd = [lean, "-j1", "-M4800"] + (["--run"] if run else []) + [str(path.resolve())]
    result = subprocess.run(cmd, cwd=PROJECT, env=env, capture_output=True, text=True)
    path.with_suffix(".log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Proof or evaluation failed. See {path.with_suffix('.log')}")
    return result.stdout


def evaluate(sweep_path, growth_path, output, paired=None, modules=(), degree=32):
    sweep = json.loads(sweep_path.read_text())
    growth = json.loads(growth_path.read_text())
    child, parent, anchor = sweep["host_progenitor"], sweep["rows"][-1], sweep["anchor"]
    if Q(parent["mass_msun"]) != 10**12 or list(map(Q, sweep["k_band"])) != [
        Q(1, 10**8),
        Q(10**10),
    ]:
        raise ValueError("Expected matched 1e12-solar-mass host and positive certified band")
    at_zero = next(r for r in growth["rows"] if Q(r["z"]) == 0)
    if paired:
        pairs = json.loads(paired.read_text())
        items = pairs if isinstance(pairs, list) else pairs["rows"]
        data = next(row for row in items if row["index"] == len(sweep["rows"]))
        if pair(data["radius"]) != pair(child["radius"]) or pair(data["parent_radius"]) != pair(
            parent["radius"]
        ):
            raise ValueError("Paired and input radii differ")
        delta = data["delta_S"]
    else:
        delta = difference(child, parent, anchor)
    if pair(delta)[0] <= 0:
        raise ValueError("The host fit requires strictly positive ΔS")
    header = "import SashimiEhHost\n" + ("import EhPaired73\n" if paired else "")
    header += "open SashimiEhHost SashimiEisensteinIntegral FloatLib.Numerics\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\n"
    base = [
        f"private def cr := {I(child['radius'])}",
        f"private def pr := {I(parent['radius'])}",
        f"private def sval := {I(at_zero['dD_dz'])}",
        f"private def dsval := {I(delta)}",
        "private theorem hc : cr.Contains childRadius := radius_check_sound (bits:=128) (degree:=40) (by decide +kernel)",
        f"private def heathRecord : SashimiHeath.Row := ⟨{I(at_zero['D'])}, sval, {I(at_zero['delta_c'])}⟩",
        f"private theorem hSlope : sval.Contains (deriv (SashimiHeath.D matched.om) 0) := by\n"
        f"  have h := (SashimiHeath.check_sound (m:=matched.om) (z:=0) (n:={growth['degree']}) (out:=heathRecord) (by decide +kernel)).2.1\n"
        "  simpa only [heathRecord, Rat.cast_zero] using h",
    ]
    if paired:
        base += [
            "private theorem hp : pr.Contains (radius matched host) := by\n"
            f"  apply radiusBox_sound (bits:={sweep['mesh']['bits']}) (degree:={sweep['mesh']['degree']})\n"
            "  decide +kernel",
            "private theorem hDS : dsval.Contains deltaS := by\n"
            "  have r := ehPairedDelta73 hc hp\n"
            "  simpa only [dsval, deltaS, SashimiSpectrum.normalize8, mul_sub, sub_div] using r",
        ]
    epochs = [Q(r["z"]) for r in growth["rows"]]
    runtime = (
        "\n".join(
            base
            + [
                "def main : IO Unit := do",
                '  for (tag,e) in [("q",qExpr),("child",childExpr),("radius",childRadiusExpr)] do',
                '    let some J := e.gridEnclosure? 128 40 | throw (IO.userError "fit interval failed")',
                '    IO.println s!"{tag} {J.lo} {J.hi}"',
                "  for z in ([" + ",".join(map(literal, epochs)) + "] : List ℚ) do",
                '    for (tag,q) in [("M200",SashimiHalo.AssemblyOutput.m200),("dM_dz",.virialRate)] do',
                f"      match (expression q z).gridEnclosureWith? [sval,dsval] 96 {degree} with",
                '      | none => IO.println s!"{tag} {z} unresolved"',
                '      | some J => IO.println s!"{tag} {z} {J.lo} {J.hi}"',
            ]
        )
        + "\n"
    )
    out = check(output.with_name(output.stem + "-measure.lean"), header + runtime, modules, True)
    model = {}
    rows = {z: {"z": str(z), "resolved": True} for z in epochs}
    theorems = list(base)
    axioms = []
    for line in out.splitlines():
        parts = line.split()
        if parts[0] in ("q", "child", "radius"):
            tag = parts[0]
            model[tag] = bounds = box(*map(Q, parts[1:]))
            rhs = {"q": "qExpr", "child": "childExpr", "radius": "childRadiusExpr"}[tag]
            name = f"host_fit_{tag}"
            theorems.append(
                f"theorem {name} : {I(bounds)}.Contains ({rhs}.eval (fun _ => 0)) :=\n"
                f"  FloatLib.Numerics.Interval.Expr.contains_gridEnclosure? "
                f"(show {rhs}.gridEnclosure? 128 40 = some {I(bounds)} from by decide +kernel)"
            )
            axioms.append(f"#print axioms {name}")
            continue
        key, z = parts[0], Q(parts[1])
        if parts[2] == "unresolved":
            rows[z].update(resolved=False)
            rows[z][key] = None
            continue
        rows[z][key] = val = box(*map(Q, parts[2:]))
        case = ".m200" if key == "M200" else ".virialRate"
        name = f"host_{key}_{epochs.index(z)}"
        arg = "" if paired else " (hDS : dsval.Contains deltaS)"
        theorems.append(
            f"theorem {name}{arg} : {I(val)}.Contains (target {case} {literal(z)}) :=\n"
            f"  eval_sound (bits:=96) (degree:={degree}) hSlope hDS (by decide +kernel)"
        )
        axioms.append(f"#print axioms {name}")
    if any("M200" not in r or "dM_dz" not in r for r in rows.values()):
        raise ValueError("Incomplete host history response")
    proof = output.with_suffix(".lean")
    log = check(proof, header + "\n".join(theorems + axioms) + "\n", modules)
    if log.count("depends on axioms: [propext, Classical.choice, Quot.sound]") != len(axioms):
        raise RuntimeError("Unexpected host proof axiom report")
    return {
        "model": "SASHIMI cold host M200 and virial-weighted signed dM/dz, zi=0, M0=1e12 solar masses",
        "units": {"M200": "Msun", "dM_dz": "Msun per unit redshift", "radius": "Mpc/h"},
        "q_Msun_convention": "q=4.137/formationRedshift(M0)^0.9476; child mass M0/q; q is dimensionless",
        "fit": model,
        "delta_S": delta,
        "mesh": sweep["mesh"],
        "kernel_replayed": bool(paired),
        "status": "Closed real host theorems"
        if paired
        else "Per-epoch output, exact fitted radius and Heath derivative kernel-checked; ΔS based on runtime continuum integration",
        "witness": proof.name,
        "rows": list(rows.values()),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sweep", "growth", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--paired", type=Path)
    parser.add_argument("--modules", type=Path, nargs="+", default=[])
    args = parser.parse_args()
    result = evaluate(args.sweep, args.growth, args.output, args.paired, args.modules)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
