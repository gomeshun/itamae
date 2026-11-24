"""Apply the same Yang and host-assembly formulas to EH quadrature and reference inputs."""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess


PROJECT = Path(__file__).resolve().parents[1]


def endpoints(row):
    if isinstance(row, dict):
        return Q(row["lo"]), Q(row["hi"])
    value = Q(row)
    return value, value


def evaluate(variance, growth, colossus):
    if not growth["kernel_replayed"]:
        raise ValueError("Heath growth intervals must be kernel-checked")
    reference = {"matched": variance["mass_rows"]}
    for backend in ("spline", "direct"):
        reference[backend] = [
            {**r[backend], "mass_msun": r["mass_msun"]} for r in colossus["variance"]
        ]
    results = {}
    lines = []
    records = []
    for backend, masses in reference.items():
        if backend == "matched":
            bg = {Q(r["z"]): r for r in growth["rows"]}
            at_q = variance["host_progenitor"]["S"]
        else:
            bg = {Q(r["z"]): r[backend] for r in colossus["background"]}
            at_q = colossus["host_progenitor"][backend]["S"]
        at_parent = next(r["S"] for r in masses if Q(r["mass_msun"]) == 10**12)
        series = {"background": [], "conditional": []}
        results[backend] = series
        for z in sorted(bg):
            row = {"z": str(z)}
            series["background"].append(row)
            args = [z, *endpoints(bg[0]["dD_dz"]), Q(at_q), Q(at_parent)]
            lines.append("host " + " ".join(map(str, args)))
            records.append((row, ("M200", "dM_dz")))
        for mass in masses:
            m = Q(mass["mass_msun"])
            if m < 10**12:
                row = {"mass_msun": str(m)}
                series["conditional"].append(row)
                args = [
                    *endpoints(bg[1]["delta_c"]),
                    *endpoints(bg[0]["delta_c"]),
                    Q(mass["S"]),
                    Q(at_parent),
                ]
                lines.append("F " + " ".join(map(str, args)))
                records.append((row, ("F",)))
    run = subprocess.run(
        ["lake", "env", "lean", "--run", "scripts/EhDownstream.lean"],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        cwd=PROJECT,
    )
    if run.returncode:
        raise RuntimeError(run.stderr + run.stdout)
    for (row, names), line in zip(records, run.stdout.splitlines(), strict=True):
        values = list(map(Q, line.split()))
        if len(values) != 2 * len(names):
            raise ValueError("Incomplete downstream response")
        for name, lo, hi in zip(names, values[::2], values[1::2], strict=True):
            if lo > hi:
                raise ValueError("Reversed downstream interval")
            row[name] = float((lo + hi) / 2)
    return {
        "description": "Midpoints for downstream formulas conditional on supplied numerical variance values. Model quadrature error is additional.",
        **results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("variance", "growth", "colossus", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    result = evaluate(
        *(json.loads(p.read_text()) for p in (args.variance, args.growth, args.colossus))
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
