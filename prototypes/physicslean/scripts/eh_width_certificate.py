"""Generate and kernel-check rational width comparisons from the matched-model publication CSV."""

import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
import io
from pathlib import Path
import re
import subprocess
import tempfile
import time

PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results" / "eps-eh"
COUNTS = {"sigma": 73, "S": 73, "dS_dM": 73, "F": 72, "M200": 29, "dM_dz": 29}
FIELDS = (
    "matched_fine_lo_exact",
    "matched_fine_hi_exact",
    "matched_tight_lo_exact",
    "matched_tight_hi_exact",
)
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}

PREAMBLE = """import Mathlib.Data.Rat.Defs

/-! Exact comparisons of the rational endpoints archived in matched-samples.csv.
The source SHA-256 is SOURCE_DIGEST. The fine and tight methods both use 7,680 cells.
This proof concerns the published endpoints. The model and integration theorems
supply the separate connection to the physical quantities. -/
namespace FixedWidthCertificate

inductive Quantity where
  | sigma | S | dS_dM | F | M200 | dM_dz
deriving DecidableEq

/-- The coordinate is physical solar mass for sigma, S, dS_dM and F.
For M200 and dM_dz it is redshift. Each method gives lower and upper bounds. -/
structure Row where
  quantity : Quantity
  coordinate : ℚ
  fineLo : ℚ
  fineHi : ℚ
  tightLo : ℚ
  tightHi : ℚ

def narrower (r : Row) := decide (r.tightHi - r.tightLo < r.fineHi - r.fineLo)
def hostAtZero (r : Row) := decide (r.quantity = .M200 ∧ r.coordinate = 0)

def Valid (r : Row) : Prop :=
  r.fineLo ≤ r.tightLo ∧ r.tightLo ≤ r.tightHi ∧ r.tightHi ≤ r.fineHi ∧
  (if hostAtZero r then
    r.fineLo = 10^12 ∧ r.fineHi = 10^12 ∧ r.tightLo = 10^12 ∧ r.tightHi = 10^12
   else r.tightHi - r.tightLo < r.fineHi - r.fineLo)

instance (r : Row) : Decidable (Valid r) := by
  unfold Valid
  split <;> infer_instance

private def Checked (rows : List Row) (size strict anchors : Nat) : Prop :=
  rows.all (fun r => decide (Valid r)) = true ∧ rows.length = size ∧
    rows.countP narrower = strict ∧ rows.countP hostAtZero = anchors

private instance (rows : List Row) (n s a : Nat) : Decidable (Checked rows n s a) := by
  unfold Checked
  infer_instance

private theorem join {x y : List Row} {nx sx ax ny sy ay : Nat}
    (hx : Checked x nx sx ax) (hy : Checked y ny sy ay) :
    Checked (x ++ y) (nx + ny) (sx + sy) (ax + ay) := by
  obtain ⟨xv,xn,xs,xa⟩ := hx
  obtain ⟨yv,yn,ys,ya⟩ := hy
  exact ⟨by simp [List.all_append,xv,yv], by simp [xn,yn],
    by simp [List.countP_append,xs,ys], by simp [List.countP_append,xa,ya]⟩

set_option maxRecDepth 100000
set_option maxHeartbeats 4000000
"""

CONCLUSION = """
/-- All 349 distinct quantity-coordinate pairs are nested. Exactly 348 have
strictly smaller width. Both intervals for M200 at z=0 are the singleton 10^12. -/
theorem publication :
    (entries.map (fun r => (r.quantity, r.coordinate))).Nodup ∧
    entries.length = 349 ∧ entries.countP narrower = 348 ∧ entries.countP hostAtZero = 1 ∧
    ∀ r ∈ entries, Valid r := by
  refine ⟨?_, allChecked.2.1, allChecked.2.2.1, allChecked.2.2.2, ?_⟩
  · decide +kernel
  · intro r hr
    exact of_decide_eq_true ((List.all_eq_true.mp allChecked.1) r hr)

#print axioms publication
end FixedWidthCertificate
"""


def read_rows(content):
    rows = []
    seen = set()
    for line, record in enumerate(csv.DictReader(io.StringIO(content.decode("utf-8"))), start=2):
        tag = record["quantity"]
        if tag not in COUNTS:
            if any(record[field] for field in FIELDS):
                raise ValueError(f"Unrecognized width comparison on CSV line {line}")
            continue
        unit = "z" if tag in {"M200", "dM_dz"} else "Msun"
        if record["unit_coordinate"] != unit or record["matched_tight_cells"] != "7680":
            raise ValueError(f"Unexpected coordinate unit or cell count on CSV line {line}")
        coordinate = Fraction(record["coordinate"])
        key = (tag, coordinate)
        if key in seen:
            raise ValueError(f"Duplicate CSV key {key}")
        seen.add(key)
        rows.append((line, tag, coordinate, *(Fraction(record[field]) for field in FIELDS)))
    if Counter(row[1] for row in rows) != COUNTS:
        raise ValueError("The six quantity counts differ from the published 349-row comparison")
    return rows


def literal(value):
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def generate(content):
    rows = read_rows(content)
    digest = hashlib.sha256(content).hexdigest()
    lines = [PREAMBLE.replace("SOURCE_DIGEST", digest)]
    nodes = []
    for start in range(0, len(rows), 16):
        block = rows[start : start + 16]
        name = f"block{len(nodes)}"
        anchors = sum(tag == "M200" and coord == 0 for _, tag, coord, *_ in block)
        lines.extend(
            [
                f"-- Matched CSV lines {block[0][0]} through {block[-1][0]}.",
                f"private def {name} : List Row := [",
            ]
        )
        for i, (_, tag, *values) in enumerate(block):
            suffix = "," if i + 1 < len(block) else "]"
            lines.append(f"  ⟨.{tag}, " + ", ".join(map(literal, values)) + "⟩" + suffix)
        lines.extend(
            [
                f"private theorem h{name} : Checked {name} {len(block)} {len(block) - anchors} {anchors} := by decide +kernel",
                "",
            ]
        )
        nodes.append((name, "h" + name))

    # Join checked chunks in exactly the same tree as the archived list.
    while len(nodes) > 1:
        merged = []
        for i in range(0, len(nodes), 2):
            x, hx = nodes[i]
            if i + 1 == len(nodes):
                merged.append((x, hx))
            else:
                y, hy = nodes[i + 1]
                merged.append((f"({x} ++ {y})", f"(join {hx} {hy})"))
        nodes = merged
    expression, proof = nodes[0]
    lines += [
        f"def entries : List Row := {expression}",
        f"private theorem allChecked : Checked entries 349 348 1 := {proof}",
        CONCLUSION,
    ]
    return "\n".join(lines)


def replay(path):
    started = time.monotonic()
    run = subprocess.run(
        ["lake", "env", "lean", "-j1", str(path)],
        cwd=PROJECT,
        capture_output=True,
        text=True,
    )
    if run.returncode:
        raise RuntimeError((run.stdout + run.stderr)[-4000:])
    audit = re.search(
        r"^'FixedWidthCertificate\.publication' depends on axioms: \[([^\]]*)\]$",
        run.stdout,
        re.MULTILINE,
    )
    if audit is None or not set(audit[1].split(", ")) <= ALLOWED_AXIOMS:
        raise RuntimeError("Missing or unexpected publication axiom audit\n" + run.stdout[-1200:])
    return time.monotonic() - started, audit[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=RESULTS / "matched-samples.csv")
    parser.add_argument("--output", type=Path, default=RESULTS / "fixed-widths.lean")
    parser.add_argument(
        "--check", action="store_true", help="compare to and replay the existing witness"
    )
    args = parser.parse_args()
    output = args.output.resolve()
    source = generate(args.csv.read_bytes()).encode("utf-8")
    if args.check:
        if output.read_bytes() != source:
            raise ValueError("The archived witness differs from its regenerated CSV source")
        seconds, audit = replay(output)
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
            candidate = Path(tmp) / output.name
            candidate.write_bytes(source)
            seconds, audit = replay(candidate)
            candidate.replace(output)
    print(audit)
    print(
        f"349 nested pairs, 348 strictly narrower, one singleton. Kernel check took {seconds:.3f} s."
    )
    print(f"Lean witness SHA-256 {hashlib.sha256(source).hexdigest()}")


if __name__ == "__main__":
    main()
