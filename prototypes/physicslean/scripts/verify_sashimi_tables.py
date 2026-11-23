"""Run W and FDM table regressions, recording exact enclosures and elapsed times.

Choose the family or the independently generated z=0 CDM input explicitly. Their content hashes also pin the h
convention, tabulated ordinates and extrapolation already present in the data.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

from spectrum_table import evaluate

WMAP_HASH = "417a0e8e9feafae8c6afb4b89e63aef80718c9641ae95994d25a6cd97438f8c9"
PLANCK_HASH = "91269da31616bf5f8d506ad138f157c8eceb284d4cd328cd99e730b882a8dd61"


def measure(path, **params):
    start = time.monotonic()
    row = evaluate(path, **params)
    if row["status"] != "enclosed":
        raise RuntimeError(f"Unresolved regression for {params}: {row!r}")
    return {"seconds": round(time.monotonic() - start, 3), **row}


def verify(w, f, reference="family"):
    f_hash = {
        "family": PLANCK_HASH,
        "independent-camb": "457dfc7caa0a4310f90aa0036a2a02d6ebd50ec62fbd15e94504cb0d06bcf84f",
    }[reference]
    for path, expected in ((w, WMAP_HASH), (f, f_hash)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"unexpected source table {path}")
    wparams = dict(
        layout="wmap7",
        interpolation="linear",
        placement="continuous",
        parameter=".07",
        bits=100,
        degree=16,
        depth=1,
    )
    # Common unnormalized spectrum, alpha=0.07 Mpc; q is the only changed input.
    rows = [
        measure(w, **wparams, transfer=model, k="7e-13", k_max=7, tolerance="0.9")
        for model in ("wdm-q5", "wdm-q10")
    ]
    if Fraction(rows[1]["upper"]) >= Fraction(rows[0]["lower"]):
        raise RuntimeError("expected to resolve convention dependence across this band")
    fparams = dict(
        layout="camb",
        h=".674",
        interpolation="log-log",
        placement="at-knots",
        k=".674",
        bits=80,
        degree=20,
        depth=1,
        tolerance="1e-5",
    )
    rows += [
        measure(f, **fparams, transfer="cdm"),
        measure(f, **fparams, transfer="fdm", parameter=1),
    ]
    if Fraction(rows[3]["upper"]) >= Fraction(rows[2]["lower"]):
        raise RuntimeError("expected to resolve FDM power suppression at the tabulated point")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wmap7", type=Path)
    parser.add_argument("planck", type=Path)
    parser.add_argument("--reference", choices=["family", "independent-camb"], default="family")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(verify(args.wmap7, args.planck, args.reference), indent=2) + "\n"
    if args.output:
        args.output.write_text(report)
    else:
        print(report, end="")
