"""Bound sigma8-normalized variance using one shared power interval per cell.

Use a sweep JSON from eh_sweep.py. Exact rational endpoint tests follow a
floating-point search. Kernel witness generation is a separate replay step.
"""

import argparse
from fractions import Fraction as Q
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from eh_full_grid import box, budget, pair, promote
from eh_plot_variance import PROJECT

SIGMA8_SQUARED = Q(811, 1000) ** 2


def write_request(sweep, ids, path):
    rows = sweep["rows"] + [sweep["host_progenitor"]]
    left = Q(sweep["k_band"][0])
    fixed = sweep["mesh"].get("enclosure_method") == "fixed_cell_enclosure"
    with path.open("w") as f:
        for item in sweep["shared_bins"]:
            lo, hi = Q(item["left"]), Q(item["right"])
            p, q = pair(item["weight"])
            if left != lo or hi <= lo or p < 0:
                raise ValueError("ordered positive-width bins with nonnegative weights required")
            if fixed:
                al, au = pair(item["integral"])
                if al < 0:
                    raise ValueError("expected nonnegative integrated weight")
                print("f", lo, hi, al, au, file=f)
            else:
                if "integral" in item:
                    raise ValueError(
                        "integrated weights require the matching fixed-cell window method"
                    )
                print("p", lo, hi, p, q, file=f)
            left = hi
        if left != Q(sweep["k_band"][1]) or list(map(Q, sweep["k_band"])) != [Q("1e-8"), Q("1e10")]:
            raise ValueError("unexpected matched spectrum band")
        for i in ids:
            print("r", i, *pair(rows[i]["radius"]), file=f)


def read_cells(path, ids, size):
    energy, anchors = [], []
    with path.open() as f:
        for _ in range(size):
            tag, el, eu, al, au = f.readline().split()
            e, a = (int(el), int(eu)), (int(al), int(au))
            if tag != "b" or not (0 <= e[0] <= e[1] and 0 <= a[0] <= a[1]):
                raise ValueError("invalid common-weight or anchor enclosure")
            energy.append(e)
            anchors.append(a)
        for i in ids:
            windows = []
            for _ in range(size):
                tag, j, sl, su = f.readline().split()
                s = int(sl), int(su)
                if tag != "w" or int(j) != i or not 0 <= s[0] <= s[1]:
                    raise ValueError("invalid or incomplete mass window enclosure")
                windows.append(s)
            yield i, energy, windows, anchors
        if f.read().strip():
            raise ValueError("unexpected trailing cells")


def residual(cells, qn, bits, qbits):
    """Return both endpoints of the exact four-corner rational rectangle sum."""
    scale = 2**qbits
    lo = hi = 0
    for (el, eu), (nl, nu), (dl, du) in cells:
        a, b = nl * scale - qn * du, nu * scale - qn * dl
        lo += min(el * a, eu * a)
        hi += max(el * b, eu * b)
    den = 2 ** (2 * bits + qbits)
    return Q(lo, den), Q(hi, den)


def propose(els, wins, dns, bits, tol_n, tol_d, limit, lower):
    # Binary64 only proposes an endpoint. The exact sign is checked by residual().
    e, w, d = (np.asarray(v, dtype=np.float64).T / float(2**bits) for v in (els, wins, dns))
    a, b = 0.0, float(limit) * 1.01
    for _ in range(48):
        q = (a + b) / 2
        if lower:
            diff = w[0] - q * d[1]
            positive = np.minimum(e[0] * diff, e[1] * diff).sum() >= q * float(tol_d)
        else:
            diff = w[1] - q * d[0]
            positive = np.maximum(e[0] * diff, e[1] * diff).sum() + float(tol_n) > 0
        if positive:
            a = q
        else:
            b = q
    return a if lower else b


def evaluate(sweep, ids, cells, qbits):
    bits = sweep["mesh"]["bits"]
    rows = sweep["rows"] + [sweep["host_progenitor"]]
    scale = 2**qbits
    results = []
    for i, e, w, d in read_cells(cells, ids, len(sweep["shared_bins"])):
        row = rows[i]
        rad = pair(row["radius"])[0]
        tn, td = budget(rad), budget(Q(8))
        baseline = promote(row["raw"], sweep["anchor"]["raw"], rad, bits)["S"]
        bl, bu = pair(baseline)
        correlated = list(zip(e, w, d))
        qs, signs = [], []
        for lower in (True, False):
            cand = propose(e, w, d, bits, tn, td, bu / SIGMA8_SQUARED, lower)
            qi = max(0, math.floor(cand * scale) - 1) if lower else math.ceil(cand * scale) + 1
            for _ in range(32):
                endpoints = residual(correlated, qi, bits, qbits)
                q = Q(qi, scale)
                admissible = endpoints[0] >= q * td if lower else endpoints[1] + tn <= 0
                if admissible:
                    break
                qi += -1 if lower else 1
                if qi < 0:
                    raise ValueError("positive lower bound unresolved")
            else:
                raise ValueError("proposed endpoint failed the exact rational sign test")
            qs.append(q)
            signs.append(box(*endpoints))
        lo, hi = (q * SIGMA8_SQUARED for q in qs)
        results.append(
            dict(
                index=i,
                mass_msun=row["mass_msun"],
                radius=row["radius"],
                ratio=box(*qs),
                residual_at_lower=signs[0],
                residual_at_upper=signs[1],
                S=box(lo, hi),
                independent_S=baseline,
                width_relative_to_independent=float((hi - lo) / (bu - bl)),
            )
        )
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sweep", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--indices", type=int, nargs="+")
    p.add_argument("--qbits", type=int, default=36)
    a = p.parse_args()
    sweep = json.loads(a.sweep.read_text())
    size = len(sweep["rows"]) + 1
    ids = list(range(size)) if a.indices is None else a.indices
    if len(set(ids)) != len(ids) or not all(0 <= i < size for i in ids) or a.qbits < 1:
        p.error("invalid mass index or precision")
    with tempfile.TemporaryDirectory() as temp:
        request = Path(temp) / "input"
        cells = Path(temp) / "cells"
        write_request(sweep, ids, request)
        with request.open() as f, cells.open("w") as out:
            subprocess.run(
                [
                    "lake",
                    "env",
                    "lean",
                    "--run",
                    "scripts/EhRatioCells.lean",
                    str(sweep["mesh"]["bits"]),
                    str(sweep["mesh"]["degree"]),
                ],
                stdin=f,
                stdout=out,
                cwd=PROJECT,
                check=True,
            )
        results = evaluate(sweep, ids, cells, a.qbits)
    report = {
        "mesh": sweep["mesh"],
        "sigma8": "811/1000",
        "qbits": a.qbits,
        "scope": "Full positive-k variance from shared power and window ranges with separate numerator/anchor tail allowances. A runtime result requires upstream spectrum and concrete residual replay for a closed kernel witness.",
        "rows": results,
    }
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
