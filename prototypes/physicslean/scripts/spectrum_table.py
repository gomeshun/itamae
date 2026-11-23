"""Evaluate a supplied SASHIMI spectrum with the kernel-proved Lean enclosure routines.

All public wavenumbers are in 1/Mpc, radii in Mpc and powers in Mpc³.
No table is extrapolated. The caller chooses the power-transfer convention.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import subprocess
import sys


if __package__:
    from .table_replay import check_table
else:
    from table_replay import check_table

PROJECT = Path(__file__).resolve().parents[1]


def number(value):
    if type(value) not in (int, str, Decimal, Fraction):
        raise ValueError("use a decimal string or an exact rational")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(f"invalid finite rational {value!r}") from exc


def _read_table(path, layout, h):
    if layout not in ("wmap7", "camb", "canonical"):
        raise ValueError(f"unknown table format {layout!r}")
    if h is not None:
        h = number(h)
    data = Path(path).read_bytes()
    lines = [line.split("#", 1)[0].strip() for line in data.decode().splitlines()]
    lines = [line for line in lines if line]
    metadata = {"sha256": hashlib.sha256(data).hexdigest(), "format": layout}
    if layout == "wmap7":
        header = [number(s) for s in lines[:5]]
        if len(header) != 5:
            raise ValueError("WMAP7 table needs five header scalars")
        if h is not None and h != header[2]:
            raise ValueError("supplied h disagrees with the WMAP7 header")
        h = header[2]
        metadata["header"] = dict(
            zip(("omega_b", "omega_m", "h", "ns", "sigma8"), map(str, header))
        )
        lines = lines[5:]
    if layout == "canonical":
        if h is not None:
            raise ValueError("canonical input is already in Mpc units")
        h = Fraction(1)
    if h is None or h <= 0:
        raise ValueError("positive --h required for CAMB k in h/Mpc")
    raw_rows = [tuple(map(number, line.split())) for line in lines]
    rows = []
    for k, power in raw_rows:
        if layout == "wmap7":
            k, power = power, k
        rows.append((k * h, power / h**3))
    if len(rows) < 2 or any(k <= 0 or power < 0 for k, power in rows):
        raise ValueError("need at least two positive wavenumbers and nonnegative powers")
    if any(left[0] >= right[0] for left, right in zip(rows, rows[1:])):
        raise ValueError("wavenumbers must increase strictly")
    metadata.update(
        h=str(h), rows=len(rows), domain_Mpc_inverse=[str(rows[0][0]), str(rows[-1][0])]
    )
    return rows, metadata, raw_rows


def load_table(path, layout, h):
    rows, metadata, _ = _read_table(path, layout, h)
    return rows, metadata


def evaluate(
    table,
    *,
    layout,
    h=None,
    k,
    k_max=None,
    tolerance,
    interpolation,
    transfer,
    placement,
    parameter=None,
    radius=Fraction(0),
    order=0,
    depth=2,
    bits=128,
    degree=24,
    sigma8=None,
    anchor_radius=None,
    anchor_k_min=None,
    anchor_k_max=None,
    anchor_depth=None,
    wdm_mass_kev=None,
    omega_wdm=None,
    thermal_h=None,
    replay=False,
):
    """Return rational endpoints or an explicit unresolved status for the entire requested band."""
    k, tolerance, radius = map(number, (k, tolerance, radius))
    k_max = None if k_max is None else number(k_max)
    parameter = None if parameter is None else number(parameter)
    if any(type(value) is not int for value in (order, depth, bits, degree)):
        raise ValueError("order, depth, bits and degree must be integers")
    if type(replay) is not bool:
        raise ValueError("replay must be a boolean")
    rows, metadata, raw_rows = _read_table(Path(table), layout, h)
    if interpolation not in ("linear", "log-log") or placement not in ("continuous", "at-knots"):
        raise ValueError("unknown interpolation or transfer placement")
    if transfer not in ("cdm", "wdm-q5", "wdm-q10", "fdm", "gaussian"):
        raise ValueError("unknown transfer convention")
    if transfer == "cdm" and parameter is not None:
        raise ValueError("CDM has no transfer parameter")
    particle_inputs = None
    wire_parameter = parameter or 0
    if any(v is not None for v in (wdm_mass_kev, omega_wdm, thermal_h)):
        if transfer not in ("wdm-q5", "wdm-q10") or parameter is not None:
            raise ValueError(
                "particle inputs require q5 or q10 and replace the supplied transfer scale"
            )
        if wdm_mass_kev is None or omega_wdm is None:
            raise ValueError("supply both --wdm-mass-kev and --omega-wdm")
        if layout == "canonical" and thermal_h is None:
            raise ValueError("canonical particle inputs require --thermal-h")
        table_h = number(metadata["h"])
        thermal_h = table_h if thermal_h is None else number(thermal_h)
        if layout != "canonical" and thermal_h != table_h:
            raise ValueError("thermal h disagrees with the table unit conversion")
        mass, abundance = map(number, (wdm_mass_kev, omega_wdm))
        if min(mass, abundance, thermal_h) <= 0:
            raise ValueError("particle mass, density and thermal h must be positive")
        wire_parameter = ",".join(map(str, (mass, abundance, thermal_h)))
        particle_inputs = dict(
            mass_keV=str(mass),
            omega_wdm=str(abundance),
            h=str(thermal_h),
            alpha_unit="Mpc",
            formula="wdmAlpha(mass_keV, omega_wdm, h) / h",
        )
    elif transfer != "cdm" and (parameter is None or parameter <= 0):
        raise ValueError("the selected transfer requires a positive parameter")
    if min(order, depth, bits, degree, radius, tolerance) < 0:
        raise ValueError("negative radius, precision, order or tolerance")
    if interpolation == "log-log" and any(p <= 0 for _, p in rows):
        raise ValueError("log-log interpolation requires positive table powers")
    upper = k if k_max is None else k_max
    if not (rows[0][0] <= k <= upper <= rows[-1][0]):
        raise ValueError("query outside the tabulated domain")
    action = "power" if k_max is None else "moment"
    if action == "power" and (radius or order):
        raise ValueError("radius and order apply only to moments")
    args = [
        action,
        interpolation,
        placement,
        transfer,
        wire_parameter,
        radius,
        k,
        upper,
        order,
        depth,
        bits,
        degree,
        tolerance,
    ]
    normalization = None
    if sigma8 is None:
        if any(v is not None for v in (anchor_radius, anchor_k_min, anchor_k_max, anchor_depth)):
            raise ValueError("anchor settings require --sigma8")
    else:
        sigma8 = number(sigma8)
        if anchor_radius is None and layout == "canonical":
            raise ValueError("canonical data needs --anchor-radius in Mpc (8/h for sigma8)")
        anchor_radius = (
            8 / number(metadata["h"]) if anchor_radius is None else number(anchor_radius)
        )
        anchor_k_min = rows[0][0] if anchor_k_min is None else number(anchor_k_min)
        anchor_k_max = rows[-1][0] if anchor_k_max is None else number(anchor_k_max)
        anchor_depth = depth if anchor_depth is None else anchor_depth
        if type(anchor_depth) is not int or sigma8 <= 0 or anchor_radius < 0 or anchor_depth < 0:
            raise ValueError("sigma8 must be positive and anchor radius and depth nonnegative")
        if not (rows[0][0] <= anchor_k_min < anchor_k_max <= rows[-1][0]):
            raise ValueError("anchor band must be nonempty and inside the table")
        args.extend([sigma8, anchor_radius, anchor_k_min, anchor_k_max, anchor_depth])
        normalization = dict(
            target_sigma=str(sigma8),
            radius_Mpc=str(anchor_radius),
            window="top-hat" if anchor_radius else "sharp-k",
            band_Mpc_inverse=[str(anchor_k_min), str(anchor_k_max)],
            full_table=anchor_k_min == rows[0][0] and anchor_k_max == rows[-1][0],
            depth=anchor_depth,
        )
    run = subprocess.run(
        ["lake", "env", "lean", "--run", "scripts/SpectrumTable.lean", *map(str, args)],
        cwd=PROJECT,
        input="".join(f"{x} {y}\n" for x, y in rows),
        text=True,
        capture_output=True,
        check=True,
    )
    result = dict(
        table=metadata,
        quantity=action,
        interpolation=interpolation,
        placement=placement,
        power_transfer=transfer,
        parameter=None if parameter is None else str(parameter),
        k_band_Mpc_inverse=[str(k), str(upper)],
        radius_Mpc=str(radius),
        moment_order=order,
        unit="Mpc^3" if action == "power" else ("1" if order == 0 else f"Mpc^-{2 * order}"),
        absolute_width_budget=str(tolerance),
        bits=bits,
        degree=degree,
        depth=depth,
        kernel_replayed=False,
    )
    if particle_inputs is not None:
        result["wdm_particle"] = particle_inputs
    if normalization is not None:
        result["normalization"] = normalization
    last = run.stdout.strip().splitlines()[-1].split()
    if last == ["unresolved"]:
        return dict(result, status="unresolved")
    if len(last) != 3 or last[0] != "enclosed":
        raise RuntimeError(f"unexpected output: {run.stdout}")
    lo, hi = map(Fraction, last[1:])
    if not lo <= hi or hi - lo > tolerance:
        raise RuntimeError("inconsistent output width")
    if replay:
        check_table(args, raw_rows, layout, number(metadata["h"]), lo, hi, PROJECT)
        result["kernel_replayed"] = True
    return dict(result, status="enclosed", lower=str(lo), upper=str(hi), width=str(hi - lo))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument(
        "--format", dest="layout", choices=["wmap7", "camb", "canonical"], required=True
    )
    parser.add_argument("--h", type=number)
    parser.add_argument(
        "--k", type=number, required=True, help="query or lower integration limit in 1/Mpc"
    )
    parser.add_argument("--k-max", type=number, help="upper finite integration limit in 1/Mpc")
    parser.add_argument("--tolerance", type=number, required=True, help="absolute enclosure width")
    parser.add_argument("--interpolation", choices=["linear", "log-log"], required=True)
    parser.add_argument(
        "--transfer", choices=["cdm", "wdm-q5", "wdm-q10", "fdm", "gaussian"], required=True
    )
    parser.add_argument("--placement", choices=["continuous", "at-knots"], required=True)
    parser.add_argument(
        "--replay",
        action="store_true",
        help="reproduce the endpoints from source rows in the Lean kernel",
    )
    parser.add_argument(
        "--wdm-mass-kev", type=number, help="thermal particle mass in keV for WDM q5 or q10"
    )
    parser.add_argument(
        "--omega-wdm", type=number, help="WDM density fraction for the particle fit"
    )
    parser.add_argument(
        "--thermal-h", type=number, help="h for the particle fit (required for canonical tables)"
    )
    parser.add_argument(
        "--parameter", type=number, help="WDM alpha [Mpc], FDM m22, or Gaussian k_fs [1/Mpc]"
    )
    parser.add_argument(
        "--radius",
        type=number,
        default=Fraction(0),
        help="top-hat radius [Mpc], zero for unwindowed band",
    )
    parser.add_argument(
        "--order", type=int, default=0, help="spectral moment j; j=0 gives variance"
    )
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--bits", type=int, default=128)
    parser.add_argument("--sigma8", type=number, help="rescale by sigma8²/raw anchor variance")
    parser.add_argument(
        "--anchor-radius",
        type=number,
        help="top-hat radius in physical Mpc (default 8/h); zero uses a sharp-k cutoff",
    )
    parser.add_argument(
        "--anchor-k-min",
        type=number,
        help="lower normalization limit in 1/Mpc (default first knot)",
    )
    parser.add_argument(
        "--anchor-k-max", type=number, help="upper normalization limit in 1/Mpc (default last knot)"
    )
    parser.add_argument(
        "--anchor-depth", type=int, help="independent subdivision depth for the normalization"
    )
    parser.add_argument("--degree", type=int, default=24)
    try:
        result = evaluate(**vars(parser.parse_args()))
    except (ValueError, ZeroDivisionError) as exc:
        parser.error(str(exc))
    except RuntimeError as exc:
        print(exc, file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as exc:
        print(exc.stderr or exc.stdout, file=sys.stderr)
        return 1
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "enclosed" else 2


if __name__ == "__main__":
    sys.exit(main())
