"""Reference values from the original SASHIMI-C classes, with prompt cusps disabled.

The family snapshot and publicly downloadable author snapshot have identical
EPS calculations for models 1--3. The family snapshot adds input validation.
"""

import argparse
from fractions import Fraction
import hashlib
import importlib.util
from importlib.metadata import version
import json
from math import isfinite
from pathlib import Path
import sys
import warnings


_SNAPSHOTS = {
    "5433bcca11ca634df881d63da7180d53ef205081d6cac39752c1d1d6ad8797b9": (
        "gomeshun/sashimi-c",
        "9334d14622a0d66082375ea9b66ce21a6d8f1ad2",
    ),
    "d75e1b0765f08cc6a8d83a5ce9ac10f6ee091954b170ff7da59003f5b43a1bf6": (
        "shinichiroando/sashimi-c",
        "e07ef8fa13503bd80818995c0d366311104fb644",
    ),
}


def _load(source):
    path = Path(source).resolve()
    if path.is_dir():
        path = path / "sashimi_c.py"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest not in _SNAPSHOTS:
        raise ValueError(f"unrecognized SASHIMI-C source at {path} (SHA-256 {digest})")
    repository, revision = _SNAPSHOTS[digest]
    spec = importlib.util.spec_from_file_location("_eps_reference_sashimi_c", path)
    module = importlib.util.module_from_spec(spec)
    # The family snapshot imports a sibling module even with tidal stripping unused.
    sys.path.insert(0, str(path.parent))
    try:
        # Upstream otherwise installs its RuntimeWarning filter process-wide.
        with warnings.catch_warnings():
            spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    return module.subhalo_properties(prompt_cusps=False), {
        "repository": repository,
        "commit": revision,
        "sha256": digest,
        "url": f"https://github.com/{repository}/blob/{revision}/sashimi_c.py",
    }


def evaluate(source, grid):
    """Run the selected SASHIMI-C methods on the supplied physical grid."""
    model, provenance = _load(source)
    host = 1e12 * model.Msun
    z_values = grid["redshifts"]
    mass_values = grid["masses_msun"]
    redshifts = {s: float(Fraction(s)) for s in z_values}
    masses = {s: float(Fraction(s)) * model.Msun for s in mass_values}
    if any(not isfinite(z) or z < 0 for z in redshifts.values()):
        raise ValueError("redshift must be finite and nonnegative")
    if any(not isfinite(m) or m <= 0 for m in masses.values()):
        raise ValueError("mass must be finite and positive")

    background = []
    for key, z in redshifts.items():
        # dMdz uses the fitted virial mass times the logarithmic M200 rate.
        background.append(
            {
                "z": key,
                "H": float(model.Hubble(z) / (model.km / model.s / model.Mpc)),
                "D": float(model.growthD(z)),
                "dD_dz": float(model.dDdz(z)),
                "delta_c": float(model.deltac_func(z)),
                "M200": float(model.Mzzi(host, z, 0.0) / model.Msun),
                "dM_dz": float(model.dMdz(host, z, 0.0) / model.Msun),
            }
        )

    variance = []
    for z_key in grid.get("variance_redshifts", z_values):
        z = float(Fraction(z_key))
        for m_key, mass in masses.items():
            sigma = float(model.sigmaMz(mass, z))
            variance.append(
                {
                    "z": z_key,
                    "mass_msun": m_key,
                    "sigma": sigma,
                    "S": sigma**2,
                    "dS_dM": float(model.dsdm(mass, z) * model.Msun),
                }
            )

    barrier0 = model.deltac_func(0.0)
    barrier1 = model.deltac_func(1.0)
    s_host = model.s_func(host)
    conditional = [
        {
            "mass_msun": key,
            "F": float(model.Ffunc_Yang(barrier0, barrier1, s_host, model.s_func(mass))),
        }
        for key, mass in masses.items()
        if mass < host
    ]
    return {
        "source": provenance,
        "versions": {"numpy": version("numpy"), "scipy": version("scipy")},
        "parameters": {
            "omega_m": model.OmegaM,
            "omega_b": model.OmegaB,
            "h": model.h,
            "host_msun": host / model.Msun,
            "host_z": 0,
            "child_z": 1,
            "prompt_cusps": False,
        },
        "methods": {
            "H": "Hubble (km/s/Mpc)",
            "D": "growthD",
            "dD_dz": "dDdz",
            "delta_c": "deltac_func",
            "M200": "Mzzi (Msun)",
            "dM_dz": "signed dMdz prescription (Msun/redshift)",
            "sigma": "sigmaMz",
            "S": "sigmaMz squared",
            "dS_dM": "signed dsdm (1/Msun)",
            "F": "Ffunc_Yang using deltac_func and s_func, per unit S",
        },
        "background": background,
        "variance": variance,
        "conditional": conditional,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", required=True, type=Path, help="SASHIMI-C checkout or sashimi_c.py"
    )
    parser.add_argument("--grid", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = evaluate(args.source, json.loads(args.grid.read_text()))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
