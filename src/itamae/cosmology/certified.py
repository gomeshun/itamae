"""Exact interval-valued cosmology using the companion PhysicsLean Lake project."""

from __future__ import annotations

import json
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Literal, TypeAlias


Exact: TypeAlias = str | int | Fraction | Decimal


@dataclass(frozen=True, slots=True)
class EnclosureResult:
    """Closed rational interval or an explicit unresolved result, in the stated unit."""

    status: Literal["enclosed", "unresolved"]
    unit: str
    tolerance: Fraction
    lower: Fraction | None
    upper: Fraction | None
    kernel_replayed: bool
    metadata: Mapping[str, Any]

    @property
    def width(self) -> Fraction | None:
        if self.upper is None or self.lower is None:
            return None
        return self.upper - self.lower


def _number(x: Exact) -> Fraction:
    if type(x) not in (str, int, Fraction, Decimal):
        raise ValueError("Use an integer, Decimal, Fraction or decimal/fraction string")
    try:
        return Fraction(x)
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(f"Invalid finite rational {x!r}") from exc


def _natural(name: str, n: int) -> str:
    if type(n) is not int or n < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return str(n)


def _result(data: dict[str, Any], tolerance: Fraction) -> EnclosureResult:
    status = data["status"]
    unit = data["unit"]
    budget = _number(data["tolerance"] if "tolerance" in data else data["absolute_width_budget"])
    replayed = data.get("kernel_replayed", False)
    if budget != tolerance or not isinstance(unit, str) or type(replayed) is not bool:
        raise ValueError("Inconsistent PhysicsLean result metadata")
    if status == "unresolved":
        if (
            any(data.get(k) is not None for k in ("interval", "lower", "upper", "width"))
            or replayed
        ):
            raise ValueError("Unresolved result contains endpoints or a concrete proof")
        return EnclosureResult(status, unit, budget, None, None, False, data)
    if status != "enclosed":
        raise ValueError(f"Unknown PhysicsLean result status {status!r}")
    interval = data.get("interval", data)
    lower, upper = _number(interval["lower"]), _number(interval["upper"])
    if not 0 <= upper - lower <= budget or _number(data["width"]) != upper - lower:
        raise ValueError("PhysicsLean result violates the requested width")
    return EnclosureResult("enclosed", unit, budget, lower, upper, replayed, data)


def _call(
    script: str,
    args: list[str],
    project: str | Path | None,
    budget: Fraction,
    unresolved_exit: int,
) -> EnclosureResult:
    if budget < 0:
        raise ValueError("tolerance must be nonnegative")
    root = (
        Path(project).expanduser()
        if project is not None
        else Path(__file__).resolve().parents[3] / "prototypes" / "physicslean"
    ).resolve()
    driver = root / "scripts" / script
    if not driver.is_file() or not (root / "lean-toolchain").is_file():
        raise FileNotFoundError(
            f"PhysicsLean driver or lean-toolchain missing in {root}. "
            "Pass project= pointing to the built prototypes/physicslean directory."
        )
    proc = subprocess.run(
        [sys.executable, str(driver), *args],
        cwd=root,
        text=True,
        capture_output=True,
    )
    if proc.returncode == 2 and not proc.stdout.strip():
        # Keep the CLI and the importable API on one set of domain and profile rules.
        for line in proc.stderr.splitlines():
            if line.startswith(f"{script}: error: "):
                raise ValueError(line.removeprefix(f"{script}: error: "))
    try:
        data = json.loads(proc.stdout)
        if isinstance(data, dict) and (
            (proc.returncode == 0 and data.get("status") == "enclosed")
            or (proc.returncode == unresolved_exit and data.get("status") == "unresolved")
        ):
            return _result(data, budget)
    except (KeyError, TypeError, ValueError) as exc:
        if proc.returncode in (0, unresolved_exit) and proc.stdout.strip():
            raise RuntimeError(f"Invalid response from {script}: {exc}") from exc
    raise RuntimeError(
        f"{script} failed ({proc.returncode}): {(proc.stderr or proc.stdout).strip()}"
    )


def evaluate_scalar(
    operation: str,
    params: Mapping[str, Exact],
    *,
    project: str | Path | None = None,
    profile: str | None = None,
    tolerance: Exact = "1e-8",
    bits: int = 96,
    degree: int = 32,
    fuel: int = 96,
    replay: bool = False,
) -> EnclosureResult:
    """Evaluate a named SASHIMI quantity with the scalar driver's units and profiles.

    The source checkout supplies the default Lake project. When using a wheel,
    set project to the matching built prototypes/physicslean directory and put
    Lake on PATH. replay checks the concrete output with Lean's kernel.
    """
    if type(replay) is not bool:
        raise ValueError("replay must be a boolean")
    if (
        not operation
        or not operation[0].isalnum()
        or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789-" for c in operation)
    ):
        raise ValueError("Invalid scalar operation name")
    if any(not key.isidentifier() for key in params):
        raise ValueError("Invalid scalar argument name")
    budget = _number(tolerance)
    args = [
        operation,
        *(f"{k}={_number(v)}" for k, v in params.items()),
        "--tol",
        str(budget),
        "--bits",
        _natural("bits", bits),
        "--degree",
        _natural("degree", degree),
        "--fuel",
        _natural("fuel", fuel),
    ]
    if profile is not None:
        args.extend(["--profile", profile])
    if replay:
        args.append("--replay")
    result = _call("sashimi.py", args, project, budget, unresolved_exit=1)
    if replay and result.status == "enclosed" and not result.kernel_replayed:
        raise RuntimeError("The scalar driver did not confirm the requested kernel replay")
    return result


def evaluate_table(
    table: str | Path,
    *,
    layout: Literal["canonical", "camb", "wmap7"],
    k: Exact,
    tolerance: Exact,
    interpolation: Literal["linear", "log-log"],
    transfer: Literal["cdm", "wdm-q5", "wdm-q10", "fdm", "gaussian"],
    placement: Literal["continuous", "at-knots"],
    project: str | Path | None = None,
    k_max: Exact | None = None,
    h: Exact | None = None,
    parameter: Exact | None = None,
    wdm_mass_kev: Exact | None = None,
    omega_wdm: Exact | None = None,
    thermal_h: Exact | None = None,
    radius: Exact = 0,
    order: int = 0,
    depth: int = 2,
    bits: int = 128,
    degree: int = 24,
    sigma8: Exact | None = None,
    anchor_radius: Exact | None = None,
    anchor_k_min: Exact | None = None,
    anchor_k_max: Exact | None = None,
    anchor_depth: int | None = None,
    replay: bool = False,
) -> EnclosureResult:
    """Evaluate power at k, or a finite-band spectral moment when k_max is supplied.

    k is in 1/Mpc and radius in Mpc. The report keeps the original table hash,
    interpolation, transfer convention and integration band. Set sigma8 to
    normalize using the finite anchor band. Canonical tables then require
    an explicit anchor_radius in Mpc. For a thermal WDM transfer, supply
    wdm_mass_kev and omega_wdm to derive alpha, with thermal_h for a canonical
    table. Alternatively, parameter supplies alpha in Mpc. replay checks the
    concrete endpoints.
    """
    if type(replay) is not bool:
        raise ValueError("replay must be a boolean")
    budget = _number(tolerance)
    args = [
        str(Path(table).expanduser().resolve()),
        "--format",
        layout,
        "--k",
        str(_number(k)),
        "--tolerance",
        str(budget),
        "--interpolation",
        interpolation,
        "--transfer",
        transfer,
        "--placement",
        placement,
        "--radius",
        str(_number(radius)),
        "--order",
        _natural("order", order),
        "--depth",
        _natural("depth", depth),
        "--bits",
        _natural("bits", bits),
        "--degree",
        _natural("degree", degree),
    ]
    for flag, val in [
        ("--k-max", k_max),
        ("--h", h),
        ("--parameter", parameter),
        ("--wdm-mass-kev", wdm_mass_kev),
        ("--omega-wdm", omega_wdm),
        ("--thermal-h", thermal_h),
        ("--sigma8", sigma8),
        ("--anchor-radius", anchor_radius),
        ("--anchor-k-min", anchor_k_min),
        ("--anchor-k-max", anchor_k_max),
    ]:
        if val is not None:
            args.extend([flag, str(_number(val))])
    if anchor_depth is not None:
        args.extend(["--anchor-depth", _natural("anchor_depth", anchor_depth)])
    if replay:
        args.append("--replay")
    result = _call("spectrum_table.py", args, project, budget, unresolved_exit=2)
    if replay and result.status == "enclosed" and not result.kernel_replayed:
        raise RuntimeError("The table driver did not confirm the requested kernel replay")
    return result


def evaluate_mass(
    table: str | Path,
    *,
    layout: Literal["canonical", "camb", "wmap7"],
    quantity: Literal[
        "sharp", "sharp-sigma", "sharp-derivative", "top-hat", "top-hat-sigma", "top-hat-derivative"
    ],
    interpolation: Literal["linear", "log-log"],
    placement: Literal["continuous", "at-knots"],
    transfer: Literal["cdm", "wdm-q5", "wdm-q10", "fdm"],
    density: Exact,
    mass: Exact,
    k_min: Exact,
    k_max: Exact,
    tolerance: Exact,
    growth: Exact | None = None,
    redshift: Exact | None = None,
    omega_m: Exact | None = None,
    omega_lambda: Exact | None = None,
    derivative_growth: Literal["squared", "legacy-f-linear"] = "squared",
    growth_tolerance: Exact | None = None,
    raw_tolerance: Exact | None = None,
    project: str | Path | None = None,
    h: Exact | None = None,
    parameter: Exact | None = None,
    wdm_mass_kev: Exact | None = None,
    omega_wdm: Exact | None = None,
    thermal_h: Exact | None = None,
    calibration: Exact = 1,
    depth: int = 2,
    bits: int = 80,
    degree: int = 16,
) -> EnclosureResult:
    """Enclose S(M), sigma(M), or dS/dM for a supplied finite-band spectrum.

    The mass is in Msun, density in Msun/Mpc³, and both wavenumbers in 1/Mpc.
    Supply an exact growth factor or a redshift and matter fraction for the
    normalized CPT growth fit. The historical F derivative convention uses
    one factor of growth. Sharp-k calibration accepts a positive rational c,
    "w-main" (c=2.5) or "w-family-legacy" (also called "sashimi-w").
    For thermal WDM supply wdm_mass_kev and omega_wdm, with thermal_h for
    canonical tables. A WDM parameter given alone is alpha in Mpc. The spectrum normalization
    is the supplied table's own. Results record the numerical policy used.
    """
    if (growth is None) == (redshift is None):
        raise ValueError("supply either growth or redshift")
    if redshift is not None and omega_m is None:
        raise ValueError("redshift growth requires omega_m")
    if growth is not None and (omega_m is not None or omega_lambda is not None):
        raise ValueError("omega_m and omega_lambda parameterize redshift growth")
    if parameter is not None and any(v is not None for v in (wdm_mass_kev, omega_wdm, thermal_h)):
        raise ValueError("WDM particle inputs replace parameter; use wdm_mass_kev for keV")
    budget = _number(tolerance)
    args = [
        str(Path(table).expanduser().resolve()),
        "--format",
        layout,
        "--quantity",
        quantity,
        "--interpolation",
        interpolation,
        "--placement",
        placement,
        "--transfer",
        transfer,
        "--density",
        str(_number(density)),
        "--mass",
        str(_number(mass)),
        "--k-min",
        str(_number(k_min)),
        "--k-max",
        str(_number(k_max)),
        "--tolerance",
        str(budget),
        "--calibration",
        str(calibration)
        if calibration in ("w-family-legacy", "sashimi-w", "w-main")
        else str(_number(calibration)),
        "--depth",
        _natural("depth", depth),
        "--bits",
        _natural("bits", bits),
        "--degree",
        _natural("degree", degree),
        "--derivative-growth",
        derivative_growth,
    ]
    for flag, val in [
        ("--h", h),
        ("--growth", growth),
        ("--redshift", redshift),
        ("--omega-m", omega_m),
        ("--omega-lambda", omega_lambda),
        ("--growth-tolerance", growth_tolerance),
        ("--raw-tolerance", raw_tolerance),
        ("--parameter", parameter),
        ("--wdm-mass-kev", wdm_mass_kev),
        ("--omega-wdm", omega_wdm),
        ("--thermal-h", thermal_h),
    ]:
        if val is not None:
            args.extend([flag, str(_number(val))])
    return _call("mass_variance.py", args, project, budget, unresolved_exit=2)
