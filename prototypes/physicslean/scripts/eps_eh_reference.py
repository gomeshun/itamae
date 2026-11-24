"""Colossus 1.4.0 references for a matched Eisenstein & Hu EPS calculation.

The default spline and the interpolation-disabled numerical integration are
evaluated independently on the same physical-mass grid. Neither route is a
proof of the accuracy of the transfer model or of the floating-point result.
"""

import argparse
import hashlib
import inspect
import json
from importlib.metadata import version
from pathlib import Path

from colossus.cosmology import cosmology, power_spectrum
from colossus.lss import peaks
from colossus.utils import constants
import numpy as np


PS = {"model": "eisenstein98"}
PARAMS = dict(
    flat=True,
    H0=67.4,
    Om0=0.315,
    Ob0=0.049,
    sigma8=0.811,
    ns=0.965,
    Tcmb0=2.7255,
    relspecies=False,
    persistence="",
)
HOST_MASS = 1.0e12
# A centered difference in ln M at consecutive half steps with Richardson
# extrapolation. Save agreement diagnostics so a quadrature jump stays visible.
STEPS = (0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125)


BANDS = {"fixed_band": (1e-8, 1e10), "wide_band": (1e-10, 1e12)}


def _sigma(cosmo, radius, band=(None, None)):
    return float(cosmo.sigma(radius, z=0.0, filt="tophat", ps_args=PS, kmin=band[0], kmax=band[1]))


def _direct_derivative(cosmo, radius, mass, band=(None, None), scale=1.0):
    """Differentiate the public, interpolation-disabled sigma(M)^2 response."""
    central = []
    richardson = []
    for step in STEPS:
        upper = _sigma(cosmo, radius * np.exp(step / 3.0), band) ** 2
        lower = _sigma(cosmo, radius * np.exp(-step / 3.0), band) ** 2
        central.append(scale * (upper - lower) / (2 * step * mass))
        if len(central) >= 2:
            richardson.append((4 * central[-1] - central[-2]) / 3)
    # Give the finest successive comparison, even if the integration routine's
    # adaptive limits or its quadrature tolerance prevent convergence.
    last = richardson[-1]
    disagreement = abs(last - richardson[-2]) / max(abs(last), np.finfo(float).tiny)
    return {
        "dS_dM": last,
        "dS_dM_per_step": central,
        "richardson_per_half_step": richardson,
        "successive_relative_difference": disagreement,
        "stable_at_1e-4_relative": bool(disagreement <= 1.0e-4),
    }


def _automatic_scan_bounds(cosmo, radius):
    """Reconstruct the unlimited-model k bounds selected by 1.4.0 _sigmaExact."""
    # Natural-log step, first/last samples outside the relative-height cutoff.
    ln_k = np.arange(
        np.log(max(1.0e-7, cosmo.k_Pk[0] * 1.0001)),
        np.log(min(1.0e25, cosmo.k_Pk[-1] * 0.9999)),
        2.0,
    )
    k = np.exp(ln_k)
    w = np.asarray([cosmo.filterFunction("tophat", x, radius) for x in k])
    f = cosmo.matterPowerSpectrum(k, z=0.0, **PS) * k**3 * w**2
    threshold = 1.0e-6 * np.max(f)
    low, high = 0, len(k) - 1
    while f[low + 1] < threshold:
        low += 1
    while f[high - 1] < threshold:
        high -= 1
    return [float(k[low]), float(k[high])]


def _variance_row(mass_value, radius, spline, direct, tight, scales):
    mass = float(mass_value)
    sig, sig_direct = _sigma(spline, radius), _sigma(direct, radius)
    slope = float(spline.sigma(radius, z=0.0, filt="tophat", derivative=True, ps_args=PS))
    heuristic_band = _automatic_scan_bounds(direct, radius)
    stable_fd = _direct_derivative(tight, radius, mass, heuristic_band)
    row = {
        "mass_msun": mass_value,
        "radius_mpc_h": float(radius),
        "radius_mpc": float(radius / spline.h),
        "spline": {"sigma": sig, "S": sig**2, "dS_dM": 2 * sig**2 * slope / (3 * mass)},
        "direct": {
            "sigma": sig_direct,
            "S": sig_direct**2,
            **stable_fd,
            "loose_fd": _direct_derivative(direct, radius, mass),
        },
        "direct_k_range_h_mpc": heuristic_band,
    }
    for mode, band in {"tight": heuristic_band, **BANDS}.items():
        scale = scales[mode]
        variance = _sigma(tight, radius, band) ** 2 * scale
        fd = stable_fd if mode == "tight" else _direct_derivative(tight, radius, mass, band)
        row[mode] = {
            "S": variance,
            "sigma": float(np.sqrt(variance)),
            "dS_dM": fd["dS_dM"] * scale,
            "successive_relative_difference": fd["successive_relative_difference"],
            "stable_at_1e-4_relative": fd["stable_at_1e-4_relative"],
        }
    return row


def _background(cosmo, direct, redshift, raw_d0):
    z = float(redshift)
    d = float(cosmo.growthFactor(z))
    direct_d = float(direct.growthFactorUnnormalized(z) / raw_d0)
    # Colossus uses D_un(a)=5/2*Om0*E(a)*int_0^a da'/(a'^3 E(a')^3)
    # for flat matter + Lambda. Differentiate the same integral analytically.
    direct_der = 1.5 * direct.Om(z) * direct_d / (1 + z) - 2.5 * direct.Om0 * (1 + z) / (
        raw_d0 * direct.Ez(z) ** 2
    )
    h = float(cosmo.Hz(z))
    splined = {
        "H": h,
        "D": d,
        "dD_dz": float(cosmo.growthFactor(z, derivative=1)),
        "delta_c": 1.686 / d,
    }
    return {
        "z": redshift,
        **splined,
        "spline": splined,
        "direct": {
            "H": float(direct.Hz(z)),
            "D": direct_d,
            "dD_dz": float(direct_der),
            "delta_c": 1.686 / direct_d,
        },
    }


def reference(redshifts, masses_msun):
    if version("colossus") != "1.4.0":
        raise ValueError("This source-matched reference requires colossus==1.4.0")
    masses = np.array(masses_msun, dtype=float)
    z = np.array(redshifts, dtype=float)
    if (
        masses.ndim != 1
        or z.ndim != 1
        or not masses.size
        or not z.size
        or not (
            np.isfinite(masses).all()
            and (masses > 0).all()
            and np.isfinite(z).all()
            and (z >= 0).all()
        )
    ):
        raise ValueError("supply positive masses and nonnegative redshifts on finite grids")
    default = cosmology.Cosmology(name="eps_eh_spline", **PARAMS)
    direct = cosmology.Cosmology(name="eps_eh_direct", interpolation=False, **PARAMS)
    tight = cosmology.Cosmology(name="eps_eh_tight", interpolation=False, **PARAMS)
    # Cache the published default power amplitude before tightening the integrator.
    amplitude = float(tight.matterPowerSpectrumNorm(**PS))
    tight.accuracy_sigma = 2.0e-8
    r8_bounds = {"tight": _automatic_scan_bounds(direct, 8.0), **BANDS}
    scales = {
        mode: (tight.sigma8 / _sigma(tight, 8.0, band)) ** 2 for mode, band in r8_bounds.items()
    }
    log_host = np.log10(HOST_MASS)
    zf = -0.0064 * log_host**2 + 0.0237 * log_host + 1.8837
    q = 4.137 / zf**0.9476
    previous = cosmology.current_cosmo
    cosmology.setCurrent(default)
    try:
        radii = peaks.lagrangianR(masses * default.h)
        host_radius = peaks.lagrangianR(HOST_MASS * default.h)
        progenitor_radius = peaks.lagrangianR((HOST_MASS / q) * default.h)
    finally:
        cosmology.setCurrent(previous)
    variance = [
        _variance_row(m, r, default, direct, tight, scales) for m, r in zip(masses_msun, radii)
    ]
    host = next((row for row in variance if float(row["mass_msun"]) == HOST_MASS), None)
    if host is None:
        host = _variance_row(HOST_MASS, host_radius, default, direct, tight, scales)
    progenitor = _variance_row(HOST_MASS / q, progenitor_radius, default, direct, tight, scales)
    progenitor.update(host_mass_msun=HOST_MASS, zf=float(zf), q=float(q))
    raw_d0 = float(direct.growthFactorUnnormalized(0.0))
    background = [_background(default, direct, z_value, raw_d0) for z_value in redshifts]
    b0, b1 = (_background(default, direct, z, raw_d0) for z in (0.0, 1.0))
    conditional = []
    for row in variance:
        if float(row["mass_msun"]) >= HOST_MASS:
            continue
        value = {"mass_msun": row["mass_msun"]}
        for mode in ("spline", "direct", "tight", *BANDS):
            ds = row[mode]["S"] - host[mode]["S"]
            if ds <= 0:
                raise ValueError("EPS conditional density requires positive variance increment")
            growth_mode = "spline" if mode == "spline" else "direct"
            db = b1[growth_mode]["delta_c"] - b0[growth_mode]["delta_c"]
            value[f"F_{mode}"] = float(
                db * np.exp(-(db**2) / (2 * ds)) / (np.sqrt(2 * np.pi) * ds**1.5)
            )
        conditional.append(value)
    transfer = []
    for k in (1.0e-6, 1.0e-4, 1.0e-2, 0.1, 1.0, 10.0, 100.0, 1.0e4, 1.0e6, 1.0e8):
        t = float(
            power_spectrum.modelEisenstein98(k, direct.h, direct.Om0, direct.Ob0, direct.Tcmb0)
        )
        raw = float(power_spectrum.powerSpectrum(k, "eisenstein98", direct, output="ps"))
        p_direct = float(direct.matterPowerSpectrum(k, **PS))
        transfer.append(
            {
                "k_h_mpc": k,
                "T": t,
                "raw_P": raw,
                "P": p_direct,
                "P_direct": p_direct,
                "P_spline": float(default.matterPowerSpectrum(k, **PS)),
                "P_fixed_band": p_direct * scales["fixed_band"],
                "P_wide_band": p_direct * scales["wide_band"],
            }
        )
    sources = {
        Path(inspect.getfile(mod)).name: hashlib.sha256(
            Path(inspect.getfile(mod)).read_bytes()
        ).hexdigest()
        for mod in (cosmology, power_spectrum, peaks, constants)
    }
    metadata = {
        "versions": {p: version(p) for p in ("colossus", "scipy", "numpy")},
        "parameters": PARAMS,
        "source_files_sha256": sources,
        "reference_url": "https://pypi.org/project/colossus/1.4.0/",
        "transfer_reference": "https://ui.adsabs.harvard.edu/abs/1998ApJ...496..605E",
        "transfer_implementation": "colossus.cosmology.power_spectrum.modelEisenstein98, including CDM and baryon/acoustic terms. Input k in h/Mpc, internal k_phys = h*k; q=k_phys/(13.41*k_eq); y=(1+z_eq)/(1+z_d); total T = (Ob/Om)*T_b+(1-Ob/Om)*T_c. The source reports CAMB comparisons to about 100 h/Mpc; higher k and the lack of a particle free-streaming cutoff are model assumptions.",
        "mass_inputs": "physical Msun; peaks.lagrangianR(h*M) returns comoving Mpc/h. R_h^3=3 h M/(4 pi Om0*rho_crit_0).",
        "rho_crit_0_msun_h2_per_mpc3": constants.RHO_CRIT_0_MPC3,
        "rho_crit_0_msun_h2_per_kpc3": constants.RHO_CRIT_0_KPC3,
        "spectrum": "P = A k^ns T^2 in (Mpc/h)^3. A comes from the raw, non-interpolated top-hat integral at 8 Mpc/h.",
        "normalization": {
            "A": amplitude,
            "sigma8": direct.sigma8,
            "direct_sigma8": _sigma(direct, 8.0),
            "spline_sigma8": _sigma(default, 8.0),
            "direct_R8_k_range_h_mpc": _automatic_scan_bounds(direct, 8.0),
        },
        "quadrature": "Colossus _sigmaExact scans log(k) with spacing 2 between ln(1e-7) and ln(0.9999e20). It brackets cells above 1e-6 of the sampled maximum of k^3 P W^2, then uses scipy.integrate.quad over log(k), epsabs=0 and epsrel=accuracy_sigma. Output includes the selected bracket for each direct radius; the tail cut is a numerical heuristic. filterFunction uses W=1 when kR<1e-3.",
        "default_sigma_quadrature_relative_tolerance": direct.accuracy_sigma,
        "sigma_radius_range_mpc_h": [default.R_min_sigma, default.R_max_sigma],
        "default_spline": "log-sigma on a radius lookup grid. The default sigma table also interpolates P(k) during its construction. Direct uses interpolation=False and the same default quadrature tolerance.",
        "P_k_region_edges_h_mpc": default.k_Pk,
        "P_k_region_bins": default.k_Pk_Nbins,
        "derivative": "dS/dM per physical Msun. Spline uses 2S (dlnsigma/dlnR)/(3M). Direct uses centered differences of the interpolation-disabled S at M exp(+-step), followed by Richardson extrapolation. A successive-step check diagnoses the numerical procedure; it is not a verified error bound.",
        "direct_fd_lnM_half_steps": list(STEPS),
        "derivative_profiles": "direct.sigma and S retain public default integration; direct.dS_dM tightens the quadrature (epsrel=2e-8) while retaining its published amplitude and holding this R's k bracket fixed across differences, as the discrete automatic k cutoffs can jump between masses. direct.loose_fd repeats the default loose quadrature and its automatic bounds at every shifted R.",
        "tight_quadrature_relative_tolerance": tight.accuracy_sigma,
        "matched_bands": {
            mode: {
                "anchor_k_h_mpc": list(band),
                "amplitude_A": amplitude * scales[mode],
                "amplitude_ratio_to_published": scales[mode],
            }
            for mode, band in r8_bounds.items()
        },
        "matched_modes": "tight uses the published heuristic bracket for each R and high-accuracy normalization at R8. fixed_band, wide_band use the listed constant k interval both in variance/derivative and R8 normalization. All three use non-interpolated Colossus EH, its top-hat implementation and direct growth.",
        "growth": "D(0)=1. Spline evaluates growthFactor and growthFactor(z,derivative=1). Direct uses growthFactorUnnormalized(z)/growthFactorUnnormalized(0) and its analytic redshift derivative for flat matter + cosmological constant.",
        "direct_Dunnormalized_z0": raw_d0,
        "conditional": "Yang Eq. 14; host 1e12 physical Msun at z=0, child z=1; common 1.686/D, with branch-specific D and S.",
    }
    return {
        "metadata": metadata,
        "background": background,
        "variance": variance,
        "host_progenitor": progenitor,
        "conditional": conditional,
        "transfer": transfer,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("grid", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    values = json.loads(args.grid.read_text())
    report = reference(values["redshifts"], values["masses_msun"])
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
