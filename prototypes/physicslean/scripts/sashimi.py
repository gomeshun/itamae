"""Evaluate the PhysicsLean scalar models. Inputs are exact decimals or fractions.

From this directory run `python scripts/sashimi.py --list` for the operation
and unit catalogue, or `python scripts/sashimi.py growth z=2 --replay`.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

if __package__:
    from .table_replay import check_interval, literal
else:
    from table_replay import check_interval, literal

PROJECT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Operation:
    module: str
    name: str
    params: tuple[str, ...]
    unit: str


OPERATIONS = {}


def register(module, rows):
    for key, lean_name, params, unit in rows:
        OPERATIONS[key] = Operation(module, lean_name, tuple(params.split()), unit)


register(
    "SashimiBackground",
    [
        ("e2", "expansionSquared", "z", "1"),
        ("hubble", "hubble", "z", "km/s/Mpc"),
        ("omega-m", "omegaMatter", "z", "1"),
        ("rho-crit", "criticalDensity", "z", "Msun/Mpc^3"),
        ("rho-m", "matterDensity", "z", "Msun/Mpc^3"),
        ("delta-vir", "virialOverdensity", "z", "1"),
    ],
)
register(
    "SashimiGrowth",
    [
        ("growth", "growth", "z", "1"),
        ("growth-derivative", "growthDerivative", "z", "1/redshift"),
        ("growth-derivative-f-legacy", "legacyFDerivative", "z", "1/redshift"),
        ("collapse-threshold", "collapseThreshold", "z", "1"),
        ("growth-ratio", "ratio", "z1 z2", "1"),
        ("growth-secant", "secant", "z1 z2", "1/redshift"),
    ],
)
register(
    "SashimiTime",
    [
        ("age", "age", "z", "Gyr"),
        ("lookback", "lookback", "z", "Gyr"),
        ("elapsed", "elapsed", "z1 z2", "Gyr"),
    ],
)
register(
    "SashimiThermal",
    [
        ("g-rho-fit", "radiation", "T", "1"),
        ("g-s-fit", "entropy", "T", "1"),
    ],
)
register(
    "SashimiSpectrumScalar",
    [
        (n, n, p, u)
        for n, p, u in [
            ("cdm-cutoff", "kfs k", "1"),
            ("wdm-alpha", "mass_keV omega_c h", "Mpc/h"),
            ("wdm-q5", "alpha k", "1"),
            ("wdm-q10", "alpha k", "1"),
            ("wdm-half-mode", "alpha", "1/Mpc"),
            ("wdm-standard-power-half", "alpha", "1/Mpc"),
            ("wdm-amplitude", "alpha k", "1"),
            ("fdm-ratio", "m22 k", "1"),
            ("top-hat", "x", "1"),
            ("sharp-k", "x", "1"),
            ("delta2", "k power", "1"),
            ("sigma", "D S0", "1"),
            ("variance", "D S0", "1"),
            ("ludlow-sigma", "pivot mass", "1"),
            ("ludlow-slope", "pivot mass", "1"),
        ]
    ],
)
register(
    "SashimiHalo",
    [
        (n, n, p, u)
        for n, p, u in [
            ("nfw-f", "x", "1"),
            ("nfw-density", "r rs rho_s", "Msun/Mpc^3"),
            ("nfw-mass", "r rs rho_s", "Msun"),
            ("nfw-potential", "r rs rho_s", "(km/s)^2"),
            ("virial-overdensity", "omega_z", "1"),
            ("overdensity-radius", "mass delta rho_crit", "Mpc"),
            ("nfw-scale-density", "c delta rho_crit", "Msun/Mpc^3"),
            ("nfw-mean-density", "c rho_s", "Msun/Mpc^3"),
            ("c200-correa", "mass z", "1"),
            ("mvir-fit", "mass c delta", "Msun"),
            ("formation-redshift", "mass", "1"),
            ("formation-q", "mass", "1"),
            ("mass-history", "mass dz D dD dS", "Msun"),
            ("mass-rate", "mvir dz D dD dS", "Msun/redshift"),
            ("cold-m200", "omega_m h mass z zi", "Msun"),
            ("cold-virial-rate", "omega_m h mass z zi", "Msun/redshift"),
            ("wdm-z-eq", "omega_m h", "1"),
            ("wdm-jeans", "omega_m h mass_keV", "Msun"),
            ("wdm-collapse", "mass jeans_mass D", "1"),
        ]
    ],
)
register(
    "SashimiDerived",
    [
        (n, n, p, u)
        for n, p, u in [
            ("scale-factor", "z", "1"),
            ("omega-lambda-flat", "omega_m z", "1"),
            ("rho-cdm-flat", "omega_m omega_b h G z", "Msun/Mpc^3"),
            ("dt-dz-flat", "omega_m h z", "Gyr/redshift"),
            ("dynamical-time", "h E delta", "Gyr"),
            ("dynamical-time-flat", "omega_m h z", "Gyr"),
            ("stripping-a-csf", "mvir z", "1"),
            ("stripping-zeta-csf", "mvir z", "1"),
            ("stripping-a-w", "mvir z", "1"),
            ("stripping-zeta-w", "mvir z", "1"),
            ("stripping-phi-csf", "mvir z h E delta", "1/redshift"),
            ("stripping-phi-w", "mvir z h E delta", "1/redshift"),
        ]
    ],
)
register(
    "SashimiHaloWarm",
    [
        ("c200-w-ludlow", "concentration", "omega_m z dS_lo dS_hi", "1"),
    ],
)
register("SashimiErf", [("erf", "erf", "x", "1")])
register(
    "SashimiHaloRoots",
    [
        ("nfw-invert-mass", "enclosed", "target lo hi", "1"),
        ("nfw-invert-mean", "mean", "target lo hi", "1"),
        ("mvir-root", "virial", "mass c delta lo hi concentration_tol", "Msun"),
        ("cvir-root", "virial-concentration", "c delta lo hi", "1"),
    ],
)

# Shared units for inputs. The operation catalogue supplies the unit of its output.
UNITS = {
    "T": "GeV",
    "G": "Mpc (km/s)^2/Msun",
    "mass": "Msun",
    "pivot": "Msun",
    "mvir": "Msun",
    "jeans_mass": "Msun",
    "rs": "Mpc",
    "r": "Mpc",
    "rho_crit": "Msun/Mpc^3",
    "rho_s": "Msun/Mpc^3",
    "m22": "10^-22 eV",
    "mass_keV": "keV",
    "alpha": "Mpc",
    "k": "1/Mpc",
    "kfs": "1/Mpc",
    "power": "Mpc^3",
    "dD": "1/redshift",
}
# h multiplies 100 km/s/Mpc. The decimal gravity is copied from native.py.
FLAT_DEFAULTS = {"omega_m": "0.315", "h": "0.674"}
WMAP_DEFAULTS = {"omega_m": "0.27", "h": "0.70"}
PROFILES = ("native-flat", "w-consistent", "w-legacy", "flrw")
COSMOLOGY = {"SashimiBackground", "SashimiGrowth", "SashimiTime"}


def rational(value):
    if type(value) not in (str, int, Fraction, Decimal):
        raise ValueError(
            "Use an integer, exact decimal string, Decimal, or Fraction for each input"
        )
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError(f"Invalid finite rational input {value!r}") from exc


def profiles_for(op):
    if op.module not in COSMOLOGY:
        return ()
    if op.module == "SashimiTime":
        return PROFILES[:2]
    return PROFILES if op.module == "SashimiBackground" else PROFILES[:3]


def catalogue():
    return {
        name: {
            "inputs": {p: UNITS.get(p, "1") for p in op.params},
            "unit": op.unit,
            "profiles": list(profiles_for(op)),
        }
        for name, op in OPERATIONS.items()
    }


def prepare(op, supplied, profile):
    supported = profiles_for(op)
    if profile is None and supported:
        profile = "native-flat"
    if profile is not None and profile not in supported:
        raise ValueError(f"Available profiles for this operation: {', '.join(supported) or 'none'}")
    defaults = {}
    if profile:
        defaults.update(FLAT_DEFAULTS if profile in ("native-flat", "flrw") else WMAP_DEFAULTS)
        if op.module == "SashimiBackground":
            defaults["G"] = "4.30091e-9"
        if profile == "w-legacy":
            defaults.update(omega_b="0.0469", omega_r="0")
        if profile == "flrw":
            defaults.update(omega_r="0", omega_k="0")
    required = set(op.params)
    if profile == "flrw":
        required.add("omega_lambda")
    unknown = set(supplied) - required - defaults.keys()
    missing = required - supplied.keys()
    if unknown or missing:
        raise ValueError(f"Unexpected inputs {sorted(unknown)}, missing inputs {sorted(missing)}")
    values = {k: rational(v) for k, v in (defaults | supplied).items()}
    if profile:
        if values["h"] <= 0:
            raise ValueError("h must be strictly positive")
        if profile != "flrw" and not 0 <= values["omega_m"] <= 1:
            raise ValueError("omega_m must lie in [0,1] in a native or WMAP7 profile")
        if profile == "w-legacy" and not (
            0 <= values["omega_b"] <= values["omega_m"] and values["omega_r"] >= 0
        ):
            raise ValueError("W legacy requires 0 ≤ omega_b ≤ omega_m and omega_r ≥ 0")
    if op.module == "SashimiHaloRoots":
        if (
            values["lo"] < 0
            or values["lo"] > values["hi"]
            or (op.name != "enclosed" and values["lo"] == 0)
        ):
            raise ValueError("Concentration bracket must be ordered and in the root's domain")
        if op.name in ("virial", "virial-concentration") and (
            values["c"] <= 0 or values["delta"] <= 0
        ):
            raise ValueError("c and delta must be strictly positive")
    return values, profile


def lean_call(op, v, profile, tol, bits, degree, fuel):
    a = {k: literal(x) for k, x in v.items()}
    mod, name = op.module, op.name
    if profile:
        m, h = a["omega_m"], a["h"]
        hubble0 = f"(100 * {h})"
        if profile == "w-legacy":
            density = f"(SashimiBackground.wLegacy {m} {a['omega_b']} {a['omega_r']})"
            if mod == "SashimiBackground":
                func = {
                    "hubble": "wLegacyHz",
                    "criticalDensity": "wLegacyHz",
                    "omegaMatter": "wLegacyOmegaz",
                    "virialOverdensity": "wLegacyOmegaz",
                }.get(name, "wLegacy")
                density = f"(SashimiBackground.{func} {m} {a['omega_b']} {a['omega_r']})"
        elif profile == "flrw":
            density = f"(⟨{a['omega_lambda']},{m},{a['omega_r']},{a['omega_k']}⟩ : SashimiBackground.Densities)"
        else:
            density = f"(SashimiBackground.flat {m})"
        if mod == "SashimiBackground":
            return (
                f"{mod}.enclose .{name} {density} {hubble0} {a['G']} {a['z']} {tol} {bits} {degree}"
            )
        if mod == "SashimiGrowth":
            lam = f"{density}.lambda"
            if name in ("ratio", "secant"):
                return (
                    f"{mod}.twoPoint? .{name} {m} {lam} {a['z1']} {a['z2']} {tol} {bits} {degree}"
                )
            return f"{mod}.enclose .{name} {m} {lam} {h} {a['z']} {tol} {bits} {degree}"
        lo, hi = {
            "age": (a.get("z"), "none"),
            "lookback": ("0", f"(some {a.get('z')})"),
            "elapsed": (a.get("z1"), f"(some {a.get('z2')})"),
        }[name]
        return f"{mod}.enclose {hubble0} {m} {mod}.gyrFactor {lo} {hi} {tol} {bits} {degree}"
    if mod == "SashimiThermal":
        return (
            f"Cosmology.EffectiveDegreesOfFreedom.encloseFit .{name} {a['T']} {tol} {bits} {degree}"
        )
    if mod == "SashimiHaloWarm":
        return (
            f"SashimiHaloWarm.concentration? {a['omega_m']} {a['z']} "
            f"⟨{a['dS_lo']},{a['dS_hi']}⟩ {bits} {degree} {tol}"
        )
    if mod == "SashimiErf":
        return f"(SashimiErf.erf? (.point {a['x']}) {bits} {degree}).bind (SashimiNumerics.within? {tol})"
    if mod == "SashimiHaloRoots":
        bracket = f"⟨{a['lo']},{a['hi']}⟩"
        if name == "virial-concentration":
            target = f"(SashimiHalo.virialTarget {a['c']} {a['delta']})"
            return f"SashimiHalo.invert? .mean {target} {bracket} {bits} {degree} {fuel} {tol}"
        if name == "virial":
            return (
                f"(SashimiHalo.virial? {a['mass']} {a['c']} {a['delta']} {bracket} "
                f"{bits} {degree} {fuel} {a['concentration_tol']} {tol}).map (·.mass)"
            )
        return f"SashimiHalo.invert? .{name} (.const {a['target']}) {bracket} {bits} {degree} {fuel} {tol}"
    func = "SashimiSpectrum.scalar?" if mod == "SashimiSpectrumScalar" else f"{mod}.evaluate?"
    args = ", ".join(a[k] for k in op.params)
    return f'{func} "{name}" [{args}] {bits} {tol} {degree}'


def run(args):
    result = subprocess.run(
        args, cwd=PROJECT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    if result.returncode:
        raise RuntimeError(f"{args[0]} failed ({result.returncode})\n{result.stdout[-3500:]}")
    return result.stdout


def evaluate(
    operation, params=None, *, profile=None, tol="1e-8", bits=96, degree=32, fuel=96, replay=False
):
    """Return exact endpoints and run metadata for one scalar operation.

    params is a mapping; all physical values are exact. Set replay=True to have
    the Lean kernel reproduce the concrete returned endpoints.
    """
    if operation not in OPERATIONS:
        raise ValueError(f"Unknown operation {operation!r}")
    for name, value in (("bits", bits), ("degree", degree), ("fuel", fuel)):
        if type(value) is not int or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if type(replay) is not bool:
        raise ValueError("replay must be a boolean")
    op = OPERATIONS[operation]
    v, profile = prepare(op, dict(params or {}), profile)
    tol = rational(tol)
    if tol < 0:
        raise ValueError("tol must be nonnegative")
    call = lean_call(op, v, profile, literal(tol), bits, degree, fuel)
    modules = [op.module] + (["SashimiBackground"] if op.module == "SashimiGrowth" else [])
    run(["lake", "build", *("+" + mod for mod in modules)])
    source = (
        "".join(f"import {mod}\n" for mod in modules)
        + f"""
open FloatLib.Numerics
def result : Option RationalInterval := {call}
def main : IO Unit := match result with
  | none => IO.println "RESULT none"
  | some i => IO.println s!"RESULT {{i.lo}} {{i.hi}}"
"""
    )
    report = {
        "operation": operation,
        "profile": profile,
        "unit": op.unit,
        "inputs": {k: str(x) for k, x in v.items()},
        "tolerance": str(tol),
        "bits": bits,
        "degree": degree,
        "kernel_replayed": False,
    }
    if op.module == "SashimiHaloRoots":
        report["fuel"] = fuel
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "Scalar.lean"
        path.write_text(source, encoding="utf-8")
        output = run(["lake", "env", "lean", "--run", str(path)])
        matches = re.findall(r"^RESULT ([^\r\n]+)$", output, flags=re.MULTILINE)
        if len(matches) != 1:
            raise RuntimeError("Lean returned no unique result")
        if matches[0] == "none":
            return report | {"status": "unresolved", "interval": None}
        lo, hi = map(Fraction, matches[0].split())
        if not 0 <= hi - lo <= tol:
            raise RuntimeError("Returned bounds violate the requested width")
        if replay:
            check_interval(source, lo, hi, PROJECT)
            report["kernel_replayed"] = True
    return report | {
        "status": "enclosed",
        "interval": {"lower": str(lo), "upper": str(hi)},
        "width": str(hi - lo),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", nargs="?", choices=sorted(OPERATIONS))
    parser.add_argument("inputs", nargs="*", metavar="NAME=VALUE")
    parser.add_argument(
        "--list", action="store_true", help="show supported operations, units and profiles"
    )
    parser.add_argument(
        "--describe", action="store_true", help="show inputs and units for the selected operation"
    )
    parser.add_argument("--profile", choices=PROFILES)
    parser.add_argument("--tol", default="1e-8", help="absolute width budget in the output unit")
    parser.add_argument("--bits", type=int, default=96)
    parser.add_argument("--degree", type=int, default=32)
    parser.add_argument("--fuel", type=int, default=96, help="maximum root bisections")
    parser.add_argument(
        "--replay", action="store_true", help="check exact endpoints again with the Lean kernel"
    )
    args = parser.parse_args()
    if args.list:
        print(json.dumps(catalogue(), indent=2))
        return 0
    if not args.operation:
        parser.error("select an operation or --list")
    if args.describe:
        print(json.dumps(catalogue()[args.operation], indent=2))
        return 0
    params = {}
    for item in args.inputs:
        key, sep, value = item.partition("=")
        if not sep or not key or key in params:
            parser.error(f"Expected a unique NAME=VALUE, got {item!r}")
        params[key] = value
    try:
        report = evaluate(
            args.operation,
            params,
            profile=args.profile,
            tol=args.tol,
            bits=args.bits,
            degree=args.degree,
            fuel=args.fuel,
            replay=args.replay,
        )
    except ValueError as exc:
        parser.error(str(exc))
    except (RuntimeError, OSError) as exc:
        parser.exit(2, f"{exc}\n")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "enclosed" else 1


if __name__ == "__main__":
    sys.exit(main())
