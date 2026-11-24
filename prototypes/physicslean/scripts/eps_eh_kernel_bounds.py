"""Kernel-check Yang conditional densities for the full matched EH98 spectrum.

Use eh_full_grid.py's output and eps_eh_growth.py's Heath intervals. Variance
differences retain the common sigma8 normalization. With --modules, the named
EhMass{index} modules and WindowCertificate close the finite-band assumptions
for those children; their module search directories are supplied in order.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

if __package__:
    from .eps_plot_kernel import box
    from .spectrum_table import PROJECT
    from .table_replay import check_proof, literal
else:
    from eps_plot_kernel import box
    from spectrum_table import PROJECT
    from table_replay import check_proof, literal

HOST = 10**12


def q(x):
    return literal(Fraction(x))


def difference(left, right):
    return {
        "lo": str(Fraction(left["lo"]) - Fraction(right["hi"])),
        "hi": str(Fraction(left["hi"]) - Fraction(right["lo"])),
    }


def inputs(full, growth):
    if (
        list(map(Fraction, full["finite_band"])) != [Fraction("1e-8"), Fraction("1e10")]
        or tuple(Fraction(full["anchor"]["radius"][k]) for k in ("lo", "hi")) != (8, 8)
        or Fraction(full["anchor"]["raw"]["lo"]) <= 0
    ):
        raise ValueError(
            "expected matched analytic EH full variance and a positive radius-8 anchor"
        )
    if (
        Fraction(growth["parameters"]["omega_m"]) != Fraction(".315")
        or Fraction(growth["parameters"]["omega_lambda"]) != Fraction(".685")
        or Fraction(growth["parameters"]["radiation"]) != 0
    ):
        raise ValueError("expected flat matter/Lambda Heath growth at Omega_m=0.315")
    epochs = {Fraction(r["z"]): r for r in growth["rows"]}
    if len(epochs) != len(growth["rows"]) or not {0, 1} <= epochs.keys():
        raise ValueError("missing or duplicate growth redshifts")
    rows = full["rows"]
    masses = [Fraction(r["mass_msun"]) for r in rows]
    if sorted(set(masses)) != masses or len(rows) < 2 or masses[0] <= 0 or masses[-1] != HOST:
        raise ValueError("need ascending unique child masses and a 1e12 Msun host")
    return rows[:-1], rows[-1], epochs


def declarations(host, anchor, epochs, degree):
    header = [
        "open FloatLib.Numerics SashimiEhYang SashimiEhFullGrid SashimiEisensteinIntegral",
        "set_option maxHeartbeats 8000000",
        "set_option maxRecDepth 100000",
        "set_option Elab.async false",
        "private noncomputable def barrier (z : ℚ) := SashimiHeath.barrier ((63/200 : ℚ) : ℝ) z",
        f"private def rp : ℚ := {q(host['radius']['lo'])}",
        f"private def p : RationalInterval := {box(host['raw'])}",
        f"private def a : RationalInterval := {box(anchor)}",
    ]
    checks = []
    for z in (0, 1):
        parts = ", ".join(box(epochs[z][key]) for key in ("D", "dD_dz", "delta_c"))
        header.append(f"private def b{z} : SashimiHeath.Row := ⟨{parts}⟩")
        checks.append(
            f"private theorem bh{z} : b{z}.delta_c.Contains (barrier {z}) := "
            f"(SashimiHeath.check_sound (m := 63/200) (z := {z}) (n := {degree}) (out := b{z}) (by decide +kernel)).2.2"
        )
    return header, checks


IMPORTS = ["import SashimiEhYang", "import SashimiHeath"]


def cell(i, row):
    if "paired" in row:
        return [f"private def ds{i} : RationalInterval := {box(row['paired']['delta_S'])}"]
    return [
        f"private def rc{i} : ℚ := {q(row['radius']['lo'])}",
        f"private def c{i} : RationalInterval := {box(row['raw'])}",
        f"private def ds{i} : RationalInterval := varianceDifference rc{i} rp c{i} p a",
    ]


def evaluate(header, children, bits, degree, split_depth):
    lines = [*IMPORTS, *header]
    for i, row in enumerate(children):
        lines.extend(cell(i, row))
    cells = ", ".join(f"ds{i}" for i in range(len(children)))
    lines.extend(
        [
            "def main : IO Unit := do",
            f"  for d in ([{cells}] : List RationalInterval) do",
            # Near the maximum in ΔS, exp and (ΔS)^(-3/2) vary in opposite directions.
            f"    let n : Nat := if 0 < d.lo ∧ d.lo ≤ 4*(b1.delta_c.hi-b0.delta_c.lo)^2 ∧ 3*d.lo ≤ 2*d.hi then {split_depth} else 0",
            f"    match densitySplit? b1.delta_c b0.delta_c {bits} {degree} n d with",
            '    | none => IO.println s!"{n} {d.lo} {d.hi} none"',
            '    | some v => IO.println s!"{n} {d.lo} {d.hi} {v.lo} {v.hi}"',
        ]
    )
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "Evaluate.lean"
        path.write_text("\n".join(lines) + "\n")
        run = subprocess.run(
            ["lake", "env", "lean", "--run", str(path)], cwd=PROJECT, text=True, capture_output=True
        )
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    if len(run.stdout.splitlines()) != len(children):
        raise RuntimeError("unexpected Lean result count")
    return [line.split() for line in run.stdout.splitlines()]


def result_claim(i, mass):
    return f"out{i}.Contains (EpsKernel.density (barrier 1 - barrier 0) (S {mass} - S {q(HOST)}))"


def case(i, row, output, bits, degree):
    mass = q(row["mass_msun"])
    common = cell(i, row)
    common.append(
        f"private theorem endpoints{i} : ds{i} = ({box(output['delta_S'])} : RationalInterval) := by decide +kernel"
    )
    split = output["variance_split_depth"]
    computation = f"densitySplit? b1.delta_c b0.delta_c {bits} {degree} {split} ds{i}"
    if "F" not in output:
        claim = f"{computation} = none"
        return [*common, f"private theorem case{i} : {claim} := by decide +kernel"], claim
    common.append(f"private def out{i} : RationalInterval := {box(output['F'])}")
    if "paired" in row:
        claim = f"ds{i}.Contains (S {mass} - S {q(HOST)}) → {result_claim(i, mass)}"
        common += [
            f"private theorem case{i} : {claim} := fun h =>",
            f"  densitySplit_sound (bits := {bits}) (degree := {degree}) (depth := {split}) bh1 bh0 h (by decide +kernel)",
        ]
        return common, claim
    claim = (
        f"(rc{i}:ℝ) ≤ radius matched {mass} → (rp:ℝ) ≤ radius matched {q(HOST)} → "
        f"c{i}.Contains (raw matched lowEdge highEdge (radius matched {mass})) → "
        f"p.Contains (raw matched lowEdge highEdge (radius matched {q(HOST)})) → "
        f"a.Contains (raw matched lowEdge highEdge 8) → {result_claim(i, mass)}"
    )
    common.extend(
        [
            f"private theorem case{i} : {claim} := by",
            "  intro hc hp cn pn an",
            "  have variance := varianceDifference_sound hc hp cn pn an (by decide +kernel)",
            f"  exact densitySplit_sound (bits := {bits}) (degree := {degree}) (depth := {split}) bh1 bh0 variance (by decide +kernel)",
        ]
    )
    return common, claim


def conjunction(claims, names, name):
    proposition = " ∧\n".join(f"({p})" for p in claims)
    term = names[0] if len(names) == 1 else "⟨" + ", ".join(names) + "⟩"
    return [f"theorem {name} : {proposition} := {term}", f"#print axioms {name}"]


def close_cases(modules, indices, children, rows, header, bits, degree, cert):
    if any(i < 0 or i >= len(rows) or rows[i]["status"] != "enclosed" for i in indices):
        raise ValueError("closed children must be resolved indices in the input grid")
    extra = [
        f"import EhPaired{i}" if "paired" in children[i] else f"import EhMass{i}" for i in indices
    ]
    lines = [*IMPORTS, *extra, *header]
    lines.extend(
        [
            "private theorem lowEdge_eq : lowEdge = (1/100000000 : ℚ) := by decide +kernel",
            "private theorem highEdge_eq : highEdge = (10000000000/1 : ℚ) := by decide +kernel",
        ]
    )
    claims = []
    for i in indices:
        local, _ = case(i, children[i], rows[i], bits, degree)
        lines.extend(local)
        claim = result_claim(i, q(children[i]["mass_msun"]))
        if "paired" in children[i]:
            row = children[i]
            lines.extend(
                [
                    f"private theorem radius{i} : ({box(row['radius'])} : RationalInterval).Contains (radius matched {q(row['mass_msun'])}) := by",
                    f"  apply radiusBox_sound (bits := {row['radius_bits']}) (degree := {row['radius_degree']}); decide +kernel",
                    f"private theorem closed{i} : {claim} := case{i} (ehPairedDelta{i} radius{i} (radiusBox_sound hEhRadius))",
                ]
            )
        else:
            lines.append(
                f"private theorem closed{i} : {claim} := case{i} ehRadius{i}.1 "
                f"(radiusBox_sound hEhRadius).1 (by rw [lowEdge_eq, highEdge_eq]; exact ehRaw{i}) "
                "(by rw [lowEdge_eq, highEdge_eq]; exact concrete.1) "
                "(by rw [lowEdge_eq, highEdge_eq]; exact concrete.2.1)"
            )
        claims.append(claim)
    lines.extend(conjunction(claims, [f"closed{i}" for i in indices], "yangComplete"))
    cert.write_text("\n".join(lines) + "\n")
    lean_path = subprocess.check_output(
        ["lake", "env", "printenv", "LEAN_PATH"], cwd=PROJECT, text=True
    ).strip()
    lean = subprocess.check_output(["lake", "env", "which", "lean"], cwd=PROJECT, text=True).strip()
    env = dict(
        os.environ, LEAN_PATH=os.pathsep.join([*(str(d.resolve()) for d in modules), lean_path])
    )
    run = subprocess.run(
        [lean, "-j1", "-M5550", str(cert.resolve())],
        cwd=PROJECT,
        env=env,
        capture_output=True,
        text=True,
    )
    log = cert.with_suffix(".log")
    log.write_text(run.stdout + run.stderr)
    if (
        run.returncode
        or "'yangComplete' depends on axioms: [propext, Classical.choice, Quot.sound]"
        not in run.stdout.splitlines()
    ):
        raise RuntimeError(f"closed EH variance-to-EPS proof or axiom audit failed; see {log}")


def generate(full, growth, output, bits, degree, modules, indices, paired=None, split_depth=5):
    if min(bits, degree, growth["degree"]) < 1 or split_depth < 0:
        raise ValueError("positive precision and degree required")
    children, host, epochs = inputs(full, growth)
    if paired is not None:
        if paired["shape_cells"] != full["mesh"]["cells"]:
            raise ValueError("paired window proof and variance mesh differ")
        seen = set()
        for ref in paired["rows"]:
            i = ref["index"]
            if i in seen or not 0 <= i < len(children):
                raise ValueError("invalid paired child index")
            seen.add(i)
            original = children[i]
            if (
                Fraction(ref["mass_msun"]) != Fraction(original["mass_msun"])
                or ref["radius"] != original["radius"]
                or ref["parent_radius"] != host["radius"]
            ):
                raise ValueError("paired child or parent radius differs from input")
            children[i] = dict(
                original,
                paired=ref,
                radius_bits=full["mesh"]["bits"],
                radius_degree=full["mesh"]["degree"],
            )
    declarations_, barriers = declarations(host, full["anchor"]["raw"], epochs, growth["degree"])
    start = time.monotonic()
    values = evaluate(declarations_, children, bits, degree, split_depth)
    elapsed = time.monotonic() - start
    proof = [*IMPORTS, *declarations_, *barriers]
    claims, rows = [], []
    for i, (child, values_) in enumerate(zip(children, values, strict=True)):
        n, *endpoints = values_
        lo, hi = map(Fraction, endpoints[:2])
        row = {
            "mass_msun": child["mass_msun"],
            "variance_split_depth": int(n),
            "delta_S": {"lo": str(lo), "hi": str(hi)},
            "end_to_end_kernel_replayed": False,
            "difference_method": (
                "paired_window_subpanels" if "paired" in child else "separate_raw_integrals"
            ),
        }
        if endpoints[2:] == ["none"]:
            row.update(
                status="unresolved",
                reason="delta_S_not_strictly_positive" if lo <= 0 else "precision_or_domain",
            )
        else:
            fl, fu = map(Fraction, endpoints[2:])
            if lo <= 0 or fl > fu:
                raise ValueError("invalid EPS interval")
            row.update(status="enclosed", F={"lo": str(fl), "hi": str(fu)})
        rows.append(row)
        lines, claim = case(i, child, row, bits, degree)
        proof.extend(lines)
        claims.append(claim)
    proof.extend(conjunction(claims, [f"case{i}" for i in range(len(rows))], "concrete"))
    cert = output.with_suffix(".lean")
    start = time.monotonic()
    statement = "\n".join(proof) + "\n"
    check_proof(statement, PROJECT, max_memory_mb=5550, timeout_seconds=900)
    cert.write_text(statement)
    replay = time.monotonic() - start
    result = {
        "definition": "Yang fixed two-epoch conditional mass fraction per dS, parent z=0, child z=1, both variances evaluated at z=0.",
        "scope": "Every Heath barrier and density propagation is kernel-replayed. A closed density additionally requires both analytic finite raw integrals, their physical radii, and the common radius-8 anchor; per-row flags record the corresponding replay.",
        "variance_target": full["model"],
        "integration_domain_k_h_per_Mpc": full["integration_domain_k_h_per_Mpc"],
        "finite_band": full["finite_band"],
        "mesh": full["mesh"],
        "difference_method": "Paired-window rows integrate the signed difference using a common spectral enclosure within each bin. Other rows subtract individually tail-padded raw integrals, multiply by sigma8 squared, then divide once by the same tail-padded radius-8 anchor.",
        "host_mass_msun": str(HOST),
        "host_S": host["S"],
        "parent_barrier": epochs[0]["delta_c"],
        "child_barrier": epochs[1]["delta_c"],
        "delta_barrier": difference(epochs[1]["delta_c"], epochs[0]["delta_c"]),
        "bits": bits,
        "degree": degree,
        "evaluation_seconds": elapsed,
        "propagation_replay_seconds": replay,
        "max_variance_split_depth": split_depth,
        "variance_split_rule": "Use 2^depth equal ΔS subintervals when 0 < lower ≤ 4 Δbarrier_upper² and upper ≥ 1.5 lower; use one interval otherwise.",
        "kernel_replayed": True,
        "kernel_replay_scope": "growth, rigorously bounded tails, common variance normalization, and density propagation",
        "certificate": cert.name,
        "rows": rows,
    }
    if modules:
        selected = (
            indices
            if indices is not None
            else [i for i, r in enumerate(rows) if r["status"] == "enclosed"]
        )
        if not selected or len(set(selected)) != len(selected):
            raise ValueError("need one or more distinct closed child indices")
        closed = output.with_name(output.stem + "-closed.lean")
        start = time.monotonic()
        close_cases(
            modules, selected, children, rows, [*declarations_, *barriers], bits, degree, closed
        )
        result["end_to_end_certificate"] = closed.name
        result["end_to_end_replay_seconds"] = time.monotonic() - start
        for i in selected:
            rows[i]["end_to_end_kernel_replayed"] = True
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("full", "growth", "output"):
        parser.add_argument(arg, type=Path)
    parser.add_argument(
        "--split-depth",
        type=int,
        default=5,
        help="Bisection depth for broad near-host ΔS intervals",
    )
    parser.add_argument("--bits", type=int, default=96)
    parser.add_argument("--degree", type=int, default=40)
    parser.add_argument("--modules", type=Path, nargs="+", default=[])
    parser.add_argument("--indices", type=int, nargs="+")
    parser.add_argument("--paired", type=Path, help="Optional child rows from eh_paired_windows.py")
    args = parser.parse_args()
    if args.indices and not args.modules:
        parser.error("--indices requires --modules")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    sources = {k: getattr(args, k).read_bytes() for k in ("full", "growth")}
    if args.paired:
        sources["paired"] = args.paired.read_bytes()
    data = generate(
        json.loads(sources["full"]),
        json.loads(sources["growth"]),
        args.output,
        args.bits,
        args.degree,
        args.modules,
        args.indices,
        json.loads(sources["paired"]) if args.paired else None,
        args.split_depth,
    )
    data["inputs_sha256"] = {k: hashlib.sha256(value).hexdigest() for k, value in sources.items()}
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
