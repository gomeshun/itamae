"""Certify signed dS/dM for the improper normalized EH98 variance.

Consumes the shared spectrum certificate and per-mass signed numerator modules
from eh_derivative_certificate.py. The optional finite-evaluation mode checks
the specialized implication with explicit finite-integral premises.
"""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from time import monotonic

if __package__:
    from .eh_full_grid import (
        C1,
        HIGH,
        LOW,
        budget,
        check_source,
        interval_expr,
        literal,
        pair,
        rounded,
    )
else:
    from eh_full_grid import (
        C1,
        HIGH,
        LOW,
        budget,
        check_source,
        interval_expr,
        literal,
        pair,
        rounded,
    )


def propose(mass, radius, finite, anchor_lo, bits):
    rmin = pair(radius)[0] / 2
    if mass <= 0 or rmin <= 0 or HIGH * rmin < 1 or anchor_lo <= 0 or bits <= 0:
        raise ValueError("Mass, radius, precision or positive denominator check failed")
    fl, fu = pair(finite)
    tail = 36 * LOW**3 + 84 * C1**2 / (rmin**3 * HIGH)
    delta = (Q(811, 1000) ** 2 * tail / mass + max(abs(fl), abs(fu)) * budget(Q(8))) / anchor_lo
    return rounded(fl - delta, fu + delta, 2**bits), rmin, tail


def certify(sweep_path, shared, slopes, destination, requested_masses=(), bits=96):
    sweep = json.loads(sweep_path.read_text())
    manifest = json.loads((slopes / "derivative-manifest.json").read_text())
    if list(map(Q, manifest["k_band"])) != [LOW, HIGH] or list(map(Q, sweep["k_band"])) != [
        LOW,
        HIGH,
    ]:
        raise ValueError("Full derivative proof requires the matched 1e-8 .. 1e10 h/Mpc band")
    if bits <= 0:
        raise ValueError("Endpoint precision must be positive")
    rows = {Q(row["mass_msun"]): row for row in sweep["rows"]}
    available = {Q(row["mass_msun"]): row for row in manifest["masses"]}
    masses = list(dict.fromkeys(map(Q, requested_masses))) or list(available)
    if not masses or any(mass not in available or mass not in rows for mass in masses):
        raise ValueError("No matching finite derivative certificate for requested mass")
    anchor_lo, _ = pair(sweep["anchor"]["raw"])
    if anchor_lo <= 0:
        raise ValueError("Expected a positive lower bound for the σ8 anchor")
    n_chunks = (manifest["cells"] - 1) // manifest["proof_chunk_size"] + 1
    root = 2 * n_chunks - 2
    destination.mkdir(parents=True, exist_ok=True)
    start = monotonic()
    results = []
    proofs = []
    imports = []

    def check_batch():
        path = destination / f"FullDerivative{imports[0]}_{imports[-1]}.lean"
        source = "".join(f"import SlopeFinal{i}\n" for i in imports)
        source += "import SashimiEhFullDerivative\nopen FloatLib.Numerics SashimiEhFullDerivative\n"
        source += "set_option maxRecDepth 100000\n"
        source += f"private theorem leftEdge : SashimiEhFullGrid.lowEdge = {literal(LOW)} := by decide +kernel\n"
        source += f"private theorem rightEdge : SashimiEhFullGrid.highEdge = {literal(HIGH)} := by decide +kernel\n"
        source += "".join(proofs)
        check_source(path, source, [slopes, shared])
        imports.clear()
        proofs.clear()

    for mass in masses:
        finite = available[mass]
        idx = finite["index"]
        if not (slopes / f"SlopeFinal{idx}.olean").is_file():
            raise FileNotFoundError(f"Missing finite slope module SlopeFinal{idx}.olean")
        box, rmin, tail = propose(mass, rows[mass]["radius"], finite["dS_dM"], anchor_lo, bits)
        theorem = f"fullCompleteDerivative{idx}"
        source = (
            f"theorem {theorem} : {interval_expr(box)}.Contains "
            f"(deriv fullS ({literal(mass)} : ℝ)) := by\n"
            f"  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) "
            f"(by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope{idx}_radiusCorrect "
            f"(by decide +kernel) (by simpa only [leftEdge] using slope{idx}_sum{root}) "
            f"(by simpa only [leftEdge,rightEdge] using concrete.2.1)\n"
            "  decide +kernel\n"
            f"#print axioms {theorem}\n"
        )
        imports.append(idx)
        proofs.append(source)
        if len(proofs) == 6:
            check_batch()
        results.append(
            {
                "mass_msun": str(mass),
                "dS_dM": box,
                "signed_raw_tail_budget": str(tail),
                "radius_lower_bound_mpc_over_h": str(rmin),
                "kernel_replayed": True,
            }
        )
    if proofs:
        check_batch()
    return {
        "target": "deriv of σ8-normalized improper EH98 top-hat variance by physical mass",
        "unit": "1/Msun",
        "integration_domain": "0 < k < infinity (h/Mpc)",
        "kernel_replayed": True,
        "cells": manifest["cells"],
        "replayed_count": len(results),
        "mesh_sha256": manifest["mesh_sha256"],
        "seconds": monotonic() - start,
        "rows": results,
    }


def promote_evaluation(sweep_path, finite_path, destination, requested_masses=(), bits=96):
    """Check each tail implication while leaving its finite-band inputs explicit."""
    sweep = json.loads(sweep_path.read_text())
    finite = json.loads(finite_path.read_text())
    rows = {Q(row["mass_msun"]): row for row in sweep["rows"]}
    slopes = {Q(row["mass_msun"]): row for row in finite["rows"] if row["status"] == "enclosed"}
    anchor = sweep["anchor"]["raw"]
    if (
        pair(anchor) != pair(finite["anchor"])
        or any(list(map(Q, doc["k_band"])) != [LOW, HIGH] for doc in (finite, sweep))
        or len(sweep["shared_bins"]) != finite["cells"]
    ):
        raise ValueError("Finite numerator and anchor metadata disagree")
    masses = list(dict.fromkeys(map(Q, requested_masses))) or list(slopes)
    if not masses or any(m not in rows or m not in slopes for m in masses):
        raise ValueError("Requested a missing or unresolved finite slope")
    start = monotonic()
    al, au = pair(anchor)
    source = [
        "import SashimiEhFullDerivative",
        "open FloatLib.Numerics SashimiEhFullDerivative SashimiEhFullGrid "
        "SashimiEisensteinIntegral SashimiEhDerivativeBounds",
    ]
    result = []
    for idx, m in enumerate(masses):
        radius = rows[m]["radius"]
        d = slopes[m]["dS_dM"]
        output, rmin, tail = propose(m, radius, d, al, bits)
        lo, hi = pair(d)
        # The runtime uses exact division by a positive anchor. The extremizing
        # denominator endpoint depends on the sign of the signed numerator.
        fl = lo * m * (al if lo < 0 else au) / Q(811, 1000) ** 2
        fu = hi * m * (au if hi < 0 else al) / Q(811, 1000) ** 2
        numerator = {"lo": str(fl), "hi": str(fu)}
        theorem = f"fullCompleteConditionalDerivative{idx}"
        q = literal(m)
        source.extend(
            [
                f"theorem {theorem}",
                f"    (hr : {interval_expr(radius)}.Contains (radius matched {q}))",
                f"    (hn : {interval_expr(numerator)}.Contains (rawSlope matched (radius matched {q}) lowEdge highEdge))",
                f"    (ha : {interval_expr(anchor)}.Contains (raw matched lowEdge highEdge 8)) :",
                f"    {interval_expr(output)}.Contains (deriv fullS ({q} : ℝ)) :=",
                "  sound hr hn ha (by decide +kernel)",
                f"#print axioms {theorem}",
            ]
        )
        result.append(
            {
                "mass_msun": str(m),
                "dS_dM": output,
                "raw_signed_finite_input": numerator,
                "signed_raw_tail_budget": str(tail),
                "radius_lower_bound_mpc_over_h": str(rmin),
                "kernel_replayed": False,
                "conditional_implication_kernel_checked": True,
            }
        )
    destination.mkdir(parents=True, exist_ok=True)
    check_source(destination / "ConditionalFullDerivative.lean", "\n".join(source) + "\n", [])
    return {
        "target": "deriv of σ8-normalized improper EH98 top-hat variance by physical mass",
        "unit": "1/Msun",
        "integration_domain": "0 < k < infinity (h/Mpc)",
        "cells": finite["cells"],
        "kernel_replayed": False,
        "conditional_implications_checked": len(result),
        "mesh_sha256": finite["mesh_sha256"],
        "seconds": monotonic() - start,
        "verification": "The Lean kernel checks each implication from its finite radius, signed numerator "
        "and σ8 anchor enclosures. The finite-cell evaluator ran without whole-grid kernel replay.",
        "rows": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sweep", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--shared", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--slopes", type=Path)
    mode.add_argument("--finite-evaluation", type=Path)
    parser.add_argument("--certificates", type=Path, required=True)
    parser.add_argument("--mass", action="append", type=Q, default=[])
    parser.add_argument("--bits", type=int, default=96)
    args = parser.parse_args()
    if args.slopes:
        if not args.shared:
            parser.error("Full replay requires --shared")
        result = certify(
            args.sweep, args.shared, args.slopes, args.certificates, args.mass, args.bits
        )
    else:
        if args.shared:
            parser.error("--shared is for --slopes replay mode")
        result = promote_evaluation(
            args.sweep, args.finite_evaluation, args.certificates, args.mass, args.bits
        )
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
