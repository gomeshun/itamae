#!/usr/bin/env python3
"""Benchmark exact certificates, interval reflection and software binary64 against modelH.

Run from the project with a pinned Lean on PATH after lake build:
  python3 scripts/benchmark.py --output results/benchmarks.json
"""

from __future__ import annotations

import argparse
import json
import platform
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time
import timeit
from fractions import Fraction
from pathlib import Path

from hz_cert import decimal_estimate, lean_rat, make_certificate, verify_certificate

PROJECT = Path(__file__).resolve().parents[1]
DIGITS = 8
BITS = 27  # 2^-27 < 10^-8
# Exact inputs, candidate output, posterior epsilon, candidate rounded radicand, forward epsilon.
CASES = [
    (
        "70",
        "3/10",
        "1",
        "8672787140647445/70368744177664",
        "2/100000000000000",
        "1745144855606067/562949953421312",
        "1/1000000000000",
    ),
    (
        "337/5",
        "63/200",
        "2",
        "1797245752001877/8796093022208",
        "5/100000000000000",
        "2586755035970929/281474976710656",
        "1/1000000000000",
    ),
    ("70", "0", "11", "70", "0", "1", "1/1000000000000"),
    ("70", "1", "3", "560", "0", "64", "2/1000000000000"),
    (
        "337/5",
        "63/200",
        "0",
        "2371426678787277/35184372088832",
        "6/1000000000000000",
        "1",
        "1/1000000000000",
    ),
]
AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


def cert_samples():
    seed = 20261124
    rng = random.Random(seed)
    for _ in range(250):
        matter_den = rng.randint(1, 1000)
        h0 = Fraction(rng.randint(1, 10000), rng.randint(1, 1000))
        matter = Fraction(rng.randint(0, matter_den), matter_den)
        z = Fraction(rng.randint(0, 10000), rng.randint(1, 1000))
        cert = make_certificate(h0, matter, z, rng.randint(0, 100))
        verify_certificate(cert)
    return {"passed": 250, "seed": seed, "checker": "Python exact rational arithmetic"}


def microseconds(fn):
    samples = timeit.repeat(fn, repeat=7, number=100)
    return round(statistics.median(samples) * 1_000_000 / 100, 3)


def run_lean(source, directory, *, args=(), reject=False):
    path = directory / "Bench.lean"
    path.write_text(source, encoding="utf-8")
    cmd = ["lake", "env", "lean"]
    if args:
        cmd.append("--run")
    cmd.extend([str(path), *args])
    start = time.perf_counter_ns()
    result = subprocess.run(
        cmd, cwd=PROJECT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8"
    )
    seconds = (time.perf_counter_ns() - start) / 1e9
    if reject:
        if result.returncode == 0 or "case hSqUpper" not in result.stdout:
            raise RuntimeError(f"Expected a false rational inequality to fail:\n{result.stdout}")
    elif result.returncode:
        raise RuntimeError(f"{path.name}\n{source}\n{result.stdout}")
    return seconds, result.stdout


def proof_times(source, directory, repeats):
    samples = []
    for _ in range(repeats):
        elapsed, output = run_lean(source + "\n#print axioms measured\n", directory)
        match = re.search(r"'measured' depends on axioms: \[([^\]]*)\]", output)
        if not match:
            raise RuntimeError(f"Missing axiom report:\n{output}")
        axioms = sorted(item.strip() for item in match[1].split(","))
        if not set(axioms) <= AXIOMS:
            raise RuntimeError(f"Unexpected proof dependencies: {axioms}")
        samples.append(round(elapsed, 6))
    return {
        "seconds": samples,
        "median_seconds": statistics.median(samples),
        "axioms": axioms,
        "source_lines_including_import": len(source.splitlines()),
    }


def enclosure_source(inputs, low, high, method):
    args = " ".join(map(lean_rat, inputs))
    target = f"{lean_rat(low)} ≤ HzCertificate.modelH {args} ∧ HzCertificate.modelH {args} ≤ {lean_rat(high)}"
    if method == "A":
        return f"""import HzCertificate

theorem measured : {target} := by
  apply HzCertificate.modelH_bounds <;> norm_num
"""
    return f"import HzInterval\n\ntheorem measured : {target} := by\n  hz_interval\n"


def signed_source(inputs, output, refined):
    args = " ".join(map(lean_rat, inputs))
    low = Fraction(refined["lower"])
    high = Fraction(refined["upper"])
    lower_error = lean_rat(output - high)
    upper_error = lean_rat(output - low)
    # Bind the signed interval to the decoded executable value, so the result
    # on the left and the candidate endpoints are both checked in Lean.
    qargs = " ".join(f"({v} : ℚ)" for v in inputs)
    val = f"(HzRoundoff.evaluate {qargs})"
    model = f"HzCertificate.modelH {args}"
    return f"""import HzRoundoff

theorem measured :
    {lower_error} ≤ {val}.toReal - {model} ∧
    {val}.toReal - {model} ≤ {upper_error} := by
  have hd : {val}.toRat? = some ({output} : ℚ) := by
    float_decide [HzRoundoff.evaluate, HzRoundoff.hubble64, HzRoundoff.radicand64]
  rw [FloatLib.Floats.ExecFloat.Binary.toReal_of_toRat?_eq_some hd]
  push_cast
  have h : {lean_rat(low)} ≤ {model} ∧ {model} ≤ {lean_rat(high)} := by
    apply HzCertificate.modelH_bounds <;> norm_num
  constructor <;> linarith [h.1, h.2]
"""


def accurate_source(inputs, output, epsilon):
    qargs = " ".join(f"({v} : ℚ)" for v in inputs)
    return f"""import HzRoundoff

theorem measured : HzRoundoff.Accurate {qargs} {lean_rat(epsilon)} := by
  certify_hz64 ({output})
"""


def forward_source(inputs, radicand, epsilon):
    qargs = " ".join(f"({v} : ℚ)" for v in inputs)
    concrete = " ".join(f"({v})" for v in inputs)
    return f"""import HzRoundoff

set_option maxHeartbeats 3000000 in
set_option exponentiation.threshold 4096 in
theorem measured : HzRoundoff.Accurate {qargs} {lean_rat(epsilon)} := by
  certify_hz64_forward {concrete} ({radicand})
"""


# Pass inputs to the executable at runtime. Time the arithmetic + finite decoding
# before formatting or printing any output. The proof timings use separate processes.
RUNTIME = """import HzRoundoff
import HzInterval

def readRat (n d : String) : Rat :=
  mkRat (Int.ofNat n.toNat!) d.toNat!

def main (args : List String) : IO Unit := do
  for row in args do
    let [hNum, hDen, mNum, mDen, zNum, zDen, bitCount] := row.splitOn "," | throw (IO.userError "bad row")
    let h := readRat hNum hDen
    let om := readRat mNum mDen
    let z := readRat zNum zDen
    let bits := bitCount.toNat!
    let t0 ← IO.monoNanosNow
    let some enclosure := HzInterval.enclose h om z bits | throw (IO.userError "enclosure failed")
    let t1 ← IO.monoNanosNow
    let some value := FloatLib.Floats.ExecFloat.Binary.toRat? (HzRoundoff.evaluate h om z)
      | throw (IO.userError "nonfinite output")
    let t2 ← IO.monoNanosNow
    IO.println s!"BENCH {t1-t0} {enclosure.lo} {enclosure.hi} {t2-t1} {value}"
"""


def runtime(directory, inputs, output):
    parts = [str(k) for v in inputs for k in (v.numerator, v.denominator)]
    row = ",".join(parts + [str(BITS)])
    # One untimed warm-up, then nine evaluations inside the same Lean interpreter.
    _, stdout = run_lean(RUNTIME, directory, args=[row] * 10)
    rows = [line.split() for line in stdout.splitlines() if line.startswith("BENCH ")]
    if len(rows) != 10:
        raise RuntimeError(f"Expected 10 runtime samples:\n{stdout}")
    _, _, lo, hi, _, _ = rows[0]
    for _, _, lower, upper, _, value in rows:
        if (lower, upper) != (lo, hi) or Fraction(value) != output:
            raise RuntimeError("Inconsistent executable results or wrong software-float candidate")
    b_us = statistics.median(int(r[1]) for r in rows[1:]) / 1000
    c_us = statistics.median(int(r[4]) for r in rows[1:]) / 1000
    return {
        "precision_bits": BITS,
        "runtime_median_us": b_us,
        "lower": lo,
        "upper": hi,
        "width": str(Fraction(hi) - Fraction(lo)),
    }, c_us


def dyadic_source(inputs, interval):
    qargs = " ".join(f"({v} : ℚ)" for v in inputs)
    h, om, z = inputs
    value = f"(⟨{interval['lower']}, {interval['upper']}⟩ : FloatLib.Numerics.Interval ℚ)"
    return f"""import HzInterval

theorem measured :
    (({value}.lo : ℝ) ≤ HzCertificate.modelH {qargs} ∧
      HzCertificate.modelH {qargs} ≤ {value}.hi) ∧
    {value}.hi - {value}.lo ≤ 1 / (2 : ℚ) ^ {BITS} := by
  exact HzInterval.enclose_sound (h₀ := ({h} : ℚ)) (Ω := ({om} : ℚ)) (z := ({z} : ℚ))
    (bits := {BITS}) (I := {value}) (by decide +kernel)
"""


def benchmark_case(case, directory, repeats):
    h0, matter, z, rounded, epsilon, radicand, forward_epsilon = map(Fraction, case)
    inputs = (h0, matter, z)

    def generate():
        return make_certificate(*inputs, DIGITS)

    cert = generate()
    low, high = Fraction(cert["lower"]), Fraction(cert["upper"])
    refined = make_certificate(*inputs, 20)
    row = {
        "inputs": {"h0": str(h0), "omega_m": str(matter), "z": str(z)},
        "model_square_exact": cert["h_squared"],
        "decimal_estimate": decimal_estimate(cert),
        "A": {
            "lower": str(low),
            "upper": str(high),
            "width": str(high - low),
            "generate_including_self_check_median_us": microseconds(generate),
            "verify_only_median_us": microseconds(lambda: verify_certificate(cert)),
            "proof": proof_times(enclosure_source(inputs, low, high, "A"), directory, repeats),
        },
    }
    dyadic, runtime_us = runtime(directory, inputs, rounded)
    row["B"] = {
        "dyadic": dyadic,
        "dyadic_kernel_check": proof_times(dyadic_source(inputs, dyadic), directory, 1),
        "same_decimal_endpoints_as_A": proof_times(
            enclosure_source(inputs, low, high, "B"), directory, repeats
        ),
    }
    row["C"] = {
        "binary64_exact_rational": str(rounded),
        "absolute_error_at_most": str(epsilon),
        "symmetric_enclosure_lower": str(rounded - epsilon),
        "symmetric_enclosure_upper": str(rounded + epsilon),
        "symmetric_enclosure_width": str(2 * epsilon),
        "runtime_median_us": runtime_us,
        "proof": proof_times(accurate_source(inputs, rounded, epsilon), directory, repeats),
        "forward_roundoff": {
            "rounded_radicand": str(radicand),
            "absolute_error_at_most": str(forward_epsilon),
            "symmetric_enclosure_width": str(2 * forward_epsilon),
            "proof": proof_times(
                forward_source(inputs, radicand, forward_epsilon), directory, repeats
            ),
        },
        "signed_error_lower": str(rounded - Fraction(refined["upper"])),
        "signed_error_upper": str(rounded - Fraction(refined["lower"])),
        "signed_error_kernel_check": proof_times(
            signed_source(inputs, rounded, refined), directory, 1
        ),
    }
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PROJECT / "results/benchmarks.json")
    parser.add_argument(
        "--repeats", type=int, default=3, help="fresh Lean processes per primary proof"
    )
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    lake = subprocess.run(
        ["lake", "env", "lean", "--version"],
        cwd=PROJECT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    manifest = json.loads((PROJECT / "lake-manifest.json").read_text())
    result = {
        "units": "km s^-1 Mpc^-1",
        "decimal_digits": DIGITS,
        "lean": lake,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "pinned_revisions": {p["name"]: p["rev"] for p in manifest["packages"] if "rev" in p},
        "methodology": {
            "proof_wall_time": (
                "Fresh lake env lean process on standalone theorem, with compiled imports. "
                "Includes startup, imports, elaboration, kernel checking and #print axioms."
            ),
            "proof_repeats": args.repeats,
            "A_and_B": (
                "Exactly the same decimal enclosure goal. B uses hz_interval (64 fractional bits). "
                "B also evaluates and kernel-checks an independently generated 27-bit enclosure."
            ),
            "C": (
                "Proof of finiteness and absolute error against exact rational inputs. "
                "The tight posterior bound checks an exact decoded value. The forward bound applies "
                "evaluate_roundoff and discharges its range, sign and error-budget conditions. "
                "Signed error (software result minus modelH) uses a 20-place rational enclosure."
            ),
            "timers": (
                "A Python timings are medians of seven 100-call batches. Candidate generation includes "
                "a self-check. B and C are medians of nine evaluations after one warm-up in the Lean "
                "interpreter, with inputs supplied at runtime. C includes input rounding and finite decoding."
            ),
        },
        "extra_certificates": cert_samples(),
        "cases": [],
    }
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        model = re.search(r"^model name\s*:\s*(.+)$", cpuinfo.read_text(), re.MULTILINE)
        if model:
            result["cpu_model"] = model[1]
    quota = Path("/sys/fs/cgroup/cpu.max")
    if quota.is_file():
        result["linux_cpu_quota_and_period_us"] = quota.read_text().strip()

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        for case in CASES:
            print(
                f"Checking h0={case[0]}, matter={case[1]}, z={case[2]}", file=sys.stderr, flush=True
            )
            result["cases"].append(benchmark_case(case, directory, args.repeats))
        # Inject an impossible upper bound after rendering. This bypasses the
        # generator and exercises Lean on exactly the template of a passing case.
        h, om, z = map(Fraction, CASES[0][:3])
        low = Fraction(result["cases"][0]["A"]["lower"])
        negative_source = enclosure_source((h, om, z), low, low, "A")
        elapsed, _ = run_lean(negative_source, directory, reject=True)
        result["corrupted_witness"] = {
            "case_index": 0,
            "replacement_upper": str(low),
            "lean_rejected": True,
            "wall_seconds": round(elapsed, 6),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
