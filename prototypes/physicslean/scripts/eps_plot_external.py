"""Floating-point comparators for the inputs to SASHIMI-C's EPS calculation.

Masses are physical solar masses. The Colossus spectrum uses Eisenstein & Hu
with baryons and a real-space spherical top hat. Astropy provides the expansion
history for the independent growing-mode integral in flat matter + Lambda.
"""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
from importlib.util import find_spec

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from colossus.cosmology import cosmology
from colossus.lss import peaks
from scipy.integrate import quad


def _grid(values, label, *, positive=False):
    converted = np.asarray(values, dtype=float)
    if converted.ndim != 1 or not len(converted) or not np.isfinite(converted).all():
        raise ValueError(f"{label} must be a nonempty finite one-dimensional grid")
    if np.any(converted <= 0 if positive else converted < 0):
        raise ValueError(f"out-of-domain {label}")
    return converted


def _growth_from_expansion(model, redshifts):
    """Normalized Heath growing mode and its analytic redshift derivative.

    For matter + Lambda, D(z) = E(z) I(a) / I(1), with
    I(a) = integral_0^a da' / (a'^3 E(a')^3). The change a'=a*x keeps
    integration on [0,1] and explicitly retains the early-time limit.
    """

    def integral(a):
        def integrand(x):
            if x == 0.0:
                return 0.0
            t = a * x
            return (1.0 / (t * model.efunc(1.0 / t - 1.0))) ** 3

        return a * quad(integrand, 0.0, 1.0, epsabs=0.0, epsrel=2e-11)[0]

    norm = integral(1.0)
    expansion = model.efunc(redshifts)
    growth = expansion * np.array([integral(1.0 / (1.0 + z)) for z in redshifts]) / norm
    derivative = 1.5 * model.Om(redshifts) * growth / (1.0 + redshifts) - (1.0 + redshifts) / (
        expansion**2 * norm
    )
    return growth, derivative


def external_curves(
    redshifts, masses_msun, *, h=0.674, omega_m=0.315, omega_b=0.049, sigma8=0.811, ns=0.965
):
    """Return redshift rows, mass rows at z=0, and the precise backend conventions."""
    z = _grid(redshifts, "redshift")
    mass = _grid(masses_msun, "mass", positive=True)
    if not np.isfinite([h, omega_m, omega_b, sigma8, ns]).all() or not (
        h > 0 and 0 < omega_b < omega_m <= 1 and sigma8 > 0 and ns > 0
    ):
        raise ValueError("invalid flat matter-plus-Lambda or spectrum parameters")
    astro = FlatLambdaCDM(H0=100 * h * u.km / u.s / u.Mpc, Om0=omega_m, Ob0=omega_b, Tcmb0=0 * u.K)
    col_params = dict(
        H0=100 * h,
        Om0=omega_m,
        Ob0=omega_b,
        flat=True,
        sigma8=sigma8,
        ns=ns,
        relspecies=False,
        persistence="",
    )
    col = cosmology.Cosmology(name="physicslean_c_eps", **col_params)
    # The peaks API reads the global cosmology. Restore even an unset state and
    # never replace a user's existing cosmology with a newly created instance.
    previous = cosmology.current_cosmo
    cosmology.setCurrent(col)
    try:
        radius = peaks.lagrangianR(mass * h)
        m8 = float(peaks.lagrangianM(8.0) / h)
        host_radius = peaks.lagrangianR(1e12 * h)
        collapse_eds = float(peaks.collapseOverdensity(corrections=False))
        collapse_lcdm = np.asarray(peaks.collapseOverdensity(corrections=True, z=z))
    finally:
        cosmology.setCurrent(previous)

    power = {"model": "eisenstein98"}
    sigma = np.asarray(col.sigma(radius, z=0.0, filt="tophat", ps_args=power))
    slope = np.asarray(col.sigma(radius, z=0.0, filt="tophat", derivative=True, ps_args=power))
    min_idx = int(np.argmin(mass))
    min_radius = float(radius[min_idx])
    direct = cosmology.Cosmology(name="physicslean_c_eps_direct", interpolation=False, **col_params)
    direct_sigma = float(direct.sigma(min_radius, z=0.0, filt="tophat", ps_args=power))
    variance = sigma**2
    log_slope = (2.0 / 3.0) * slope
    ds_dm = variance * log_slope / mass
    growth = np.array([col.growthFactor(float(redshift)) for redshift in z])
    d_growth = np.asarray(col.growthFactor(z, derivative=1))
    ast_growth, ast_d_growth = _growth_from_expansion(astro, z)
    ast_h = astro.H(z).to_value(u.km / u.s / u.Mpc)
    ast_age = astro.age(z).to_value(u.Gyr)
    col_h = col.Hz(z)

    col_background = [
        {
            "z": redshifts[i],
            "H": float(col_h[i]),
            "D": float(growth[i]),
            "dD_dz": float(d_growth[i]),
            "delta_c": float(1.686 / growth[i]),
            "delta_c_colossus_eds": float(collapse_eds / growth[i]),
            "delta_c_colossus_corrected": float(collapse_lcdm[i] / growth[i]),
        }
        for i in range(len(z))
    ]
    astro_background = [
        {
            "z": redshifts[i],
            "H": float(ast_h[i]),
            "D": float(ast_growth[i]),
            "dD_dz": float(ast_d_growth[i]),
            "delta_c": float(1.686 / ast_growth[i]),
            "age_gyr": float(ast_age[i]),
        }
        for i in range(len(z))
    ]
    masses = [
        {
            "mass_msun": masses_msun[i],
            "radius_mpc_h": float(radius[i]),
            "radius_mpc": float(radius[i] / h),
            "sigma": float(sigma[i]),
            "S": float(variance[i]),
            "dS_dM": float(ds_dm[i]),
            "dlnS_dlnM": float(log_slope[i]),
        }
        for i in range(len(mass))
    ]
    host_variance = float(col.sigma(host_radius, z=0.0, filt="tophat", ps_args=power) ** 2)
    delta_barrier = 1.686 * (
        1.0 / float(col.growthFactor(1.0)) - 1.0 / float(col.growthFactor(0.0))
    )
    children = np.flatnonzero(mass < 1.0e12)
    difference = variance[children] - host_variance
    if np.any(difference <= 0):
        raise ValueError("Eq. 14 requires positive child minus host variance")
    density = (
        delta_barrier
        * np.exp(-(delta_barrier**2) / (2 * difference))
        / (np.sqrt(2 * np.pi) * difference**1.5)
    )
    conditional = [{"mass_msun": masses_msun[i], "F": float(f)} for i, f in zip(children, density)]
    metadata = {
        "versions": {p: version(p) for p in ("astropy", "colossus", "numpy", "scipy")},
        "cosmology": {
            "H0_km_s_Mpc": 100 * h,
            "Om0": omega_m,
            "Ob0": omega_b,
            "Ode0": 1 - omega_m,
            "sigma8_colossus": sigma8,
            "ns_colossus": ns,
            "radiation_in_background": False,
        },
        "colossus": {
            "api": [
                "colossus.cosmology.cosmology.Cosmology.Hz",
                "colossus.cosmology.cosmology.Cosmology.growthFactor",
                "colossus.cosmology.cosmology.Cosmology.sigma",
                "colossus.lss.peaks.lagrangianR",
                "colossus.lss.peaks.lagrangianM",
                "colossus.lss.peaks.collapseOverdensity",
            ],
            "colossus_m_module_available": find_spec("colossus.m") is not None,
            "spectrum": "linear matter; eisenstein98; tophat; Colossus spline and k limits",
            "units": "lagrangianR input is h*M [Msun/h]; R is Mpc/h. dS/dM=2*S*(dlnsigma/dlnR)/(3*M) uses physical M [Msun].",
            "radius_interpolation": {
                "default_domain_mpc_h": [col.R_min_sigma, col.R_max_sigma],
                "requested_domain_mpc_h": [min_radius, float(np.max(radius))],
                "smallest_mass_msun": float(mass[min_idx]),
                "smallest_radius_physical_pc": min_radius * 1.0e6 / h,
                "reciprocal_smallest_radius_h_per_mpc": 1.0 / min_radius,
                "direct_integral_sigma_at_smallest_radius": direct_sigma,
                "spline_sigma_relative_to_direct_minus_one": float(
                    sigma[min_idx] / direct_sigma - 1.0
                ),
                "direct_api": "Cosmology(interpolation=False).sigma; same spectrum and spherical top-hat.",
            },
            "small_scale_assumptions": "The Eisenstein & Hu CDM+baryon transfer function is continued to small scales; particle free-streaming and a microhalo cutoff are absent. A reciprocal radius is a characteristic scale of the broad top-hat window. Agreement of the spline and direct integral tests numerical interpolation of this model.",
            "collapse_eds": collapse_eds,
            "collapse_corrected": "Colossus flat-universe delta_c,EdS * Om(z)^0.0055; barriers extrapolated to z=0 by dividing by D(z).",
            "background_cmb": "relspecies=False; Tcmb0=2.7255 K used by the transfer function",
            "M_at_8_Mpc_h_msun": m8,
            "spline_sigma_at_8_Mpc_h": float(col.sigma(8.0, z=0.0, filt="tophat", ps_args=power)),
        },
        "astropy": {
            "api": [
                "astropy.cosmology.FlatLambdaCDM.H",
                "astropy.cosmology.FlatLambdaCDM.efunc",
                "astropy.cosmology.FlatLambdaCDM.Om",
                "astropy.cosmology.FlatLambdaCDM.age",
                "scipy.integrate.quad",
            ],
            "growth": "Heath matter+Lambda growing mode from Astropy E(z) using SciPy quadrature, D(0)=1. Astropy supplies the expansion history.",
            "mass_variance": "No native linear density spectrum or sigma(M) supplied by FlatLambdaCDM.",
        },
        "conditional": {
            "source": "Yang et al. (2011), Eq. 14; analytic conditional first-crossing density per unit child variance.",
            "reference": "https://arxiv.org/abs/1104.1757",
            "host_mass_msun": 1e12,
            "host_z": 0,
            "child_z": 1,
            "host_variance": host_variance,
            "delta_barrier": delta_barrier,
            "construction": "Own Colossus top-hat S(M,0) and common 1.686/D with Colossus growth. This is an analytic post-processing of Colossus inputs.",
        },
        "source_links": {
            "colossus_cosmology": "https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html",
            "colossus_peaks": "https://bdiemer.bitbucket.io/colossus/lss_peaks.html",
            "astropy": "https://docs.astropy.org/en/stable/api/astropy.cosmology.FlatLambdaCDM.html",
        },
        "comparison": "The SASHIMI-C calculation has a fixed-amplitude Ludlow fit without independent sigma8/ns inputs. The common EPS barrier uses 1.686; separate Colossus thresholds report its own collapse constants.",
    }
    return {
        "astropy": {"background": astro_background},
        "colossus": {"background": col_background, "variance": masses, "conditional": conditional},
        "metadata": metadata,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    grid = json.loads(args.grid.read_text())
    curves = external_curves(grid["redshifts"], grid["masses_msun"])
    args.output.write_text(json.dumps(curves, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
