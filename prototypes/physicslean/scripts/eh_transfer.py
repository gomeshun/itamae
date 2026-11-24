"""Evaluate Colossus EH98 k samples with native Lean and certify selected points.

python scripts/eh_transfer.py colossus.json transfer.json
lake env lean transfer.lean
"""

import argparse
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import struct
import subprocess
import tempfile

if __package__:
    from .table_replay import check_proof, literal
else:
    from table_replay import check_proof, literal

PROJECT = Path(__file__).resolve().parents[1]


def calculate(reference, output, bits, degree, certified_k_max):
    grid = [Fraction(str(item["k_h_mpc"])) for item in reference["transfer"]]
    if len(grid) != len(set(grid)) or not grid or min(grid) <= 0 or min(bits, degree) <= 0:
        raise ValueError("positive, unique k samples and positive precision required")
    imports = "import SashimiEisensteinHu\nopen SashimiEisensteinHu FloatLib.Numerics\n"
    params = "0.674 0.315 0.049 2.7255"
    source = (
        imports
        + f"""
def main : IO Unit := do
  let a := prepareFloat {params}
  let some bs := coefficientBoxes? {params} {bits} {degree} | throw (IO.userError "coefficient evaluation failed")
  for b in bs do IO.println s!"BOX {{b.lo}} {{b.hi}}"
  for k in ([{", ".join(map(literal, grid))}] : List ℚ) do
    let t := transferFloat a (Float.ofInt k.num / Float.ofNat k.den)
    IO.println s!"NATIVE {{k}} {{t.toBits.toNat}}"
    if k ≤ {literal(certified_k_max)} then
      let some i := rangeFast? bs (.point k) {bits} {degree} | throw (IO.userError s!"unresolved T({{k}})")
      IO.println s!"INTERVAL {{k}} {{i.lo}} {{i.hi}}"
"""
    )
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "Eval.lean"
        path.write_text(source)
        run = subprocess.run(
            ["lake", "env", "lean", "--run", str(path)], cwd=PROJECT, capture_output=True, text=True
        )
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    boxes, rows = [], {}
    for line in run.stdout.splitlines():
        fields = line.split()
        if fields[:1] == ["BOX"]:
            boxes.append(tuple(map(Fraction, fields[1:])))
        elif fields[:1] == ["NATIVE"]:
            k = Fraction(fields[1])
            (value,) = struct.unpack(">d", int(fields[2]).to_bytes(8, "big"))
            rows[k] = dict(
                k_h_mpc=str(k), T=value, native_encoding="binary64", kernel_replayed=False
            )
        elif fields[:1] == ["INTERVAL"]:
            rows[Fraction(fields[1])].update(lo=fields[2], hi=fields[3])
    if len(boxes) != 11 or list(rows) != grid:
        raise ValueError("missing or misordered EH samples")

    def box(lo, hi):
        return f"⟨{literal(lo)}, {literal(hi)}⟩"

    proof = [
        imports,
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 20000000",
        "set_option Elab.async false",
        "private def bs : List RationalInterval := [" + ", ".join(box(*b) for b in boxes) + "]",
        f"private theorem certifiedBoxes : coefficientBoxes? {params} {bits} {degree} = some bs := by decide +kernel",
        f"private noncomputable def model := prepareReal {params}",
    ]
    claims, terms = [], []
    for k, row in rows.items():
        if k > certified_k_max:
            continue
        lo, hi = (Fraction(row[name]) for name in ("lo", "hi"))
        if lo > hi:
            raise ValueError("reversed EH enclosure")
        name = str(len(terms))
        proof += [
            f"private def box{name} : RationalInterval := {box(lo, hi)}",
            f"private theorem point{name} : box{name}.Contains (transferReal model {literal(k)}) :=",
            f"  rangeFast_sound (RationalInterval.contains_point {literal(k)})",
            f"    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := {bits}) (degree := {degree})",
            "    (by decide +kernel)",
        ]
        claims.append(f"box{name}.Contains (transferReal model {literal(k)})")
        terms.append(f"point{name}")
        row["kernel_replayed"] = True
    if not claims:
        raise ValueError("the certification cutoff must cover at least one requested k")
    term = "⟨" + ", ".join(terms) + "⟩" if len(terms) > 1 else terms[0]
    proof.append(
        "theorem concrete : " + " ∧\n  ".join(claims) + " := " + term + "\n#print axioms concrete\n"
    )
    proof = "\n".join(proof)
    check_proof(proof, PROJECT, max_memory_mb=4096, timeout_seconds=600)
    certificate = output.with_suffix(".lean")
    certificate.write_text(proof)
    output.write_text(
        json.dumps(
            dict(
                model="Colossus 1.4.0 modelEisenstein98 with baryons (Eisenstein & Hu 1998)",
                k_unit="h/Mpc",
                parameters=dict(h=".674", omega_m=".315", omega_b=".049", Tcmb_K="2.7255"),
                native="Lean Float model, returned as IEEE binary64 bits",
                bits=bits,
                degree=degree,
                certificate=certificate.name,
                certified_k_max=str(certified_k_max),
                rows=list(rows.values()),
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--bits", type=int, default=96)
    p.add_argument("--degree", type=int, default=40)
    p.add_argument("--certified-k-max", type=Fraction, default=Fraction(10**8))
    a = p.parse_args()
    calculate(
        json.loads(a.reference.read_text(), parse_float=Decimal),
        a.output,
        a.bits,
        a.degree,
        a.certified_k_max,
    )


if __name__ == "__main__":
    main()
