"""Mass-dependent finite-band S and dS/dM for a supplied linear spectrum.

Canonical input units are M_sun, M_sun/Mpc³ and 1/Mpc. The supplied growth
factor is relative to the spectrum epoch. Alternatively --redshift and
--omega-m evaluate the real CPT fit for an unscaled z=0 input spectrum. Both
paths preserve the supplied table normalization. Historical F single-D mass
derivative scaling is available via --derivative-growth legacy-f-linear.
"""

import argparse
import json
from fractions import Fraction
from pathlib import Path
import subprocess
from spectrum_table import PROJECT, load_table, number


def evaluate_mass(
    table,
    *,
    layout,
    h=None,
    quantity,
    interpolation,
    placement,
    transfer,
    density,
    mass,
    k_min,
    k_max,
    tolerance,
    growth=None,
    parameter=None,
    wdm_mass_kev=None,
    omega_wdm=None,
    thermal_h=None,
    calibration=1,
    depth=2,
    bits=80,
    degree=16,
    redshift=None,
    omega_m=None,
    omega_lambda=None,
    derivative_growth="squared",
    raw_tolerance=None,
    growth_tolerance=None,
):
    rows, source = load_table(Path(table), layout, h)
    a, b, rho, m, tol = map(number, (k_min, k_max, density, mass, tolerance))
    if (growth is None) == (redshift is None):
        raise ValueError("choose exactly one of growth and redshift")
    if quantity not in (
        "sharp",
        "sharp-derivative",
        "sharp-sigma",
        "top-hat",
        "top-hat-derivative",
        "top-hat-sigma",
    ):
        raise ValueError("unknown quantity")
    derivative = quantity.endswith("-derivative")
    if derivative_growth not in ("squared", "legacy-f-linear") or (
        not derivative and derivative_growth != "squared"
    ):
        raise ValueError("legacy-f-linear is available for a signed mass derivative only")
    if derivative:
        output = "legacy-f-linear" if derivative_growth == "legacy-f-linear" else "derivative"
    else:
        output = "sigma" if quantity.endswith("-sigma") else "variance"
    d = number(growth) if growth is not None else 1
    raw_tol = number(raw_tolerance) if raw_tolerance is not None else tol
    post = (
        redshift is not None or output in ("legacy-f-linear", "sigma") or raw_tolerance is not None
    )
    z = om = ol = gt = None
    if redshift is not None:
        if omega_m is None:
            raise ValueError("CPT growth requires omega_m")
        z, om = map(number, (redshift, omega_m))
        ol = 1 - om if omega_lambda is None else number(omega_lambda)
        gt = Fraction(1, 10000) if growth_tolerance is None else number(growth_tolerance)
        if z < 0 or not 0 < om <= 1 or ol < 0 or gt < 0:
            raise ValueError("invalid CPT redshift, density parameters or growth tolerance")
    elif omega_m is not None or omega_lambda is not None or growth_tolerance is not None:
        raise ValueError("CPT parameters require redshift")
    if raw_tol < 0:
        raise ValueError("raw tolerance must be nonnegative")
    if interpolation not in ("linear", "log-log") or placement not in ("continuous", "at-knots"):
        raise ValueError("unknown interpolation or transfer placement")
    if interpolation == "log-log" and any(p <= 0 for _, p in rows):
        raise ValueError("log-log interpolation requires positive input powers")
    if transfer not in ("cdm", "fdm", "wdm-q5", "wdm-q10"):
        raise ValueError("unknown transfer")
    if not (rows[0][0] <= a < b <= rows[-1][0]) or min(rho, m) <= 0 or min(d, tol) < 0:
        raise ValueError("invalid band, density, mass, growth, or width")
    if any(type(x) is not int or x < 0 for x in (depth, bits, degree)):
        raise ValueError("depth, bits, degree must be nonnegative integers")
    if transfer == "cdm" and parameter is not None:
        raise ValueError("CDM takes no parameter")
    if any(v is not None for v in (wdm_mass_kev, omega_wdm, thermal_h)):
        if transfer not in ("wdm-q5", "wdm-q10") or parameter is not None:
            raise ValueError("WDM particle inputs replace --parameter, the transfer scale in Mpc")
        if wdm_mass_kev is None or omega_wdm is None:
            raise ValueError("thermal WDM requires --wdm-mass-kev and --omega-wdm")
        if layout == "canonical" and thermal_h is None:
            raise ValueError("canonical particle inputs require --thermal-h")
        table_h = number(source["h"])
        hh = table_h if thermal_h is None else number(thermal_h)
        if layout != "canonical" and hh != table_h:
            raise ValueError("thermal h disagrees with the table unit conversion")
        p, omega = map(number, (wdm_mass_kev, omega_wdm))
        if min(p, omega, hh) <= 0:
            raise ValueError("thermal WDM mass, Omega_WDM and h must be positive")
    else:
        if transfer != "cdm" and (parameter is None or number(parameter) <= 0):
            raise ValueError("positive transfer parameter required")
        p = number(parameter) if parameter is not None else 0
        omega, hh = 0, 0
    supplied_calibration = str(calibration)
    if calibration == "sashimi-w":
        calibration = "w-family-legacy"
    elif calibration not in ("w-family-legacy", "w-main"):
        calibration = number(calibration)
        if calibration <= 0:
            raise ValueError("positive sharp-k calibration required")
    if quantity.startswith("top-hat") and calibration != 1:
        raise ValueError("top-hat uses the physical mass radius")
    op = quantity.removesuffix("-sigma")
    inputs = [
        op,
        interpolation,
        placement,
        transfer,
        p,
        omega,
        hh,
        calibration,
        rho,
        m,
        d,
        a,
        b,
        tol,
        depth,
        bits,
        degree,
    ]
    if post:
        inputs.extend(
            [
                "cpt" if redshift is not None else "fixed",
                output,
                om or 0,
                ol or 0,
                z or 0,
                gt or 0,
                raw_tol,
            ]
        )
    run = subprocess.run(
        ["lake", "env", "lean", "--run", "scripts/MassVariance.lean", *map(str, inputs)],
        cwd=PROJECT,
        text=True,
        input="".join(f"{k} {v}\n" for k, v in rows),
        capture_output=True,
    )
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    parts = run.stdout.strip().splitlines()[-1].split()
    result = dict(
        table=source,
        quantity=quantity,
        unit="1/Msun" if "derivative" in quantity else "1",
        mass_Msun=str(m),
        density_Msun_per_Mpc3=str(rho),
        growth_factor=None if z is not None else str(d),
        growth_model="cpt" if z is not None else "supplied",
        derivative_growth=derivative_growth if derivative else None,
        redshift=None if z is None else str(z),
        growth_omega_m=None if om is None else str(om),
        growth_omega_lambda=None if ol is None else str(ol),
        reference_spectrum_epoch="z=0, unit growth" if z is not None else "caller-defined",
        growth_tolerance=None if gt is None else str(gt),
        raw_tolerance=str(raw_tol) if post else None,
        growth_power=1 if output in ("legacy-f-linear", "sigma") else 2,
        k_band_Mpc_inverse=list(map(str, (a, b))),
        tail_model="none; finite band",
        normalization="as supplied",
        calibration=str(calibration),
        requested_calibration=supplied_calibration,
        interpolation=interpolation,
        placement=placement,
        transfer=transfer,
        parameter=None if parameter is None else str(p),
        wdm_mass_kev=str(p) if omega else None,
        omega_wdm=str(omega) if omega else None,
        thermal_h=str(hh) if hh else None,
        bits=bits,
        depth=depth,
        degree=degree,
        tolerance=str(tol),
        verification="proved evaluator; Lean runtime execution",
    )
    if parts == ["unresolved"]:
        return dict(result, status="unresolved")
    if len(parts) != 3 or parts[0] != "enclosed":
        raise RuntimeError(f"invalid Lean response {run.stdout!r}")
    lo, hi = map(Fraction, parts[1:])
    if lo > hi or hi - lo > tol:
        raise RuntimeError("invalid returned width")
    return dict(result, status="enclosed", lower=str(lo), upper=str(hi), width=str(hi - lo))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument(
        "--format", dest="layout", choices=["wmap7", "camb", "canonical"], required=True
    )
    parser.add_argument("--h", type=number)
    parser.add_argument(
        "--quantity",
        choices=[
            "sharp",
            "top-hat",
            "sharp-sigma",
            "top-hat-sigma",
            "sharp-derivative",
            "top-hat-derivative",
        ],
        required=True,
    )
    parser.add_argument("--interpolation", choices=["linear", "log-log"], required=True)
    parser.add_argument("--placement", choices=["continuous", "at-knots"], required=True)
    parser.add_argument("--transfer", choices=["cdm", "wdm-q5", "wdm-q10", "fdm"], required=True)
    parser.add_argument(
        "--parameter",
        type=number,
        help="WDM alpha [Mpc] or FDM m22",
    )
    parser.add_argument("--wdm-mass-kev", type=number, help="thermal WDM particle mass in keV")
    parser.add_argument("--omega-wdm", type=number)
    parser.add_argument("--thermal-h", type=number)
    parser.add_argument(
        "--calibration",
        default="1",
        help="positive sharp-k c, w-main (5/2), or w-family-legacy ((5/2)(9π/2)^(1/3)); sashimi-w aliases w-family-legacy",
    )
    epoch = parser.add_mutually_exclusive_group(required=True)
    epoch.add_argument(
        "--growth", type=number, help="exact supplied growth relative to the input power"
    )
    epoch.add_argument("--redshift", type=number, help="apply CPT to a reference z=0 spectrum")
    parser.add_argument("--omega-m", type=number, help="CPT present-day matter fraction")
    parser.add_argument("--omega-lambda", type=number, help="default 1 minus omega-m")
    parser.add_argument(
        "--raw-tolerance",
        type=number,
        help="budget for unscaled integral; defaults to final tolerance",
    )
    parser.add_argument(
        "--growth-tolerance", type=number, help="CPT factor width budget; default 1e-4"
    )
    parser.add_argument(
        "--derivative-growth", choices=["squared", "legacy-f-linear"], default="squared"
    )
    for flag in ("density", "mass", "k-min", "k-max", "tolerance"):
        parser.add_argument("--" + flag, type=number, required=True)
    for flag, default in (("depth", 2), ("bits", 80), ("degree", 16)):
        parser.add_argument("--" + flag, type=int, default=default)
    try:
        result = evaluate_mass(**vars(parser.parse_args()))
    except (ValueError, ZeroDivisionError) as err:
        parser.error(str(err))
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "enclosed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
