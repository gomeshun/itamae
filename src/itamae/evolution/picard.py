"""Shared log-mass Picard integration with model-supplied coefficients.

Numerics derived from gomeshun/sashimi-c's standalone Picard implementation
(main 2ffefc8f02ff1db2f634aae73c1b2a394914a98e). Physical rates, host histories,
units and any validated domain are supplied by the caller, never by ITAMAE.
"""

from __future__ import annotations

from dataclasses import dataclass
import time
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .ode import solve_evolution
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline


class PicardFallbackWarning(UserWarning):
    """An explicitly counted high-accuracy ODE fallback was required.

    Notes
    -----
    A guarded Picard request required direct DOP853 integration. Inspect diagnostics() for
    reason and affected mass count; a completed fallback is not evidence that the table met
    the requested iteration criterion.
    """


def _real(value, name):
    """Validate real numerical input before its use in a physical calculation.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.

    Returns
    -------
    numpy.ndarray
        Numerical array with the original shape; the helper rejects the invalid values
        checked in its implementation.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    a = np.asarray(value)
    if a.dtype.kind not in "iuf" or np.any(~np.isfinite(a)):
        raise ValueError(f"{name} must contain finite real numbers.")
    return a.astype(float, copy=False)


def _integer(value, name, minimum):
    """Require a valid grid or iteration count before allocating a tidal table.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.
    minimum : int
        Smallest supported node/iteration count for this specific numerical control.

    Returns
    -------
    int
        Validated integer at least minimum.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    if (
        isinstance(value, (bool, np.bool_))
        or not np.isscalar(value)
        or not np.isfinite(value)
        or int(value) != value
        or value < minimum
    ):
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(value)


def physics_key(solver):
    """Read the complete physical-state identity that controls tidal-table reuse.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.

    Returns
    -------
    object
        Immutable equality-comparable callback identity for the current solver.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    return solver.state_key()


def _state(solver):
    """Invalidate all cached tidal coefficients and tables when the physical state changes.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.

    Returns
    -------
    object
        Current physical-state key after cache synchronization.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    key = physics_key(solver)
    if getattr(solver, "_picard_physics_key", None) != key:
        solver._picard_tables = {}
        solver._picard_coefficients = {}
        solver._picard_host_cache = {}
        solver._picard_events = []
        solver._picard_histories = {}
        solver._picard_physics_key = key
    return key


def _coefficients(solver, z_acc, z_final, count):
    """Sample the tidal rate, exponent and host history on one redshift path.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.
    z_acc : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    z_final : float
        Scalar final redshift, no greater than z_acc; count nodes sample the path from z_acc
        to this endpoint.
    count : int
        Number of nodes on a uniformly spaced redshift integration path, including both
        endpoints.

    Returns
    -------
    tuple of numpy.ndarray
        Redshift, nonnegative Phi, finite zeta and natural-log host mass arrays, each of
        length count.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    _state(solver)
    key = (float(z_acc), float(z_final), count)
    cached = solver._picard_coefficients.get(key)
    if cached is not None:
        return cached
    z = np.linspace(z_acc, z_final, count)
    phi = np.broadcast_to(np.asarray(solver.Phi(z), dtype=float), z.shape)
    zeta = np.broadcast_to(np.asarray(solver.zetaMz(z), dtype=float), z.shape)
    host = np.broadcast_to(np.asarray(solver.Mzvir(z), dtype=float), z.shape)
    if (
        not np.all(np.isfinite(phi))
        or not np.all(np.isfinite(zeta))
        or not np.all(np.isfinite(host))
        or np.any(host <= 0)
        or np.any(phi < 0)
    ):
        raise ValueError(
            "Tidal stripping requires finite, positive host mass and nonnegative rate."
        )
    result = z, phi, zeta, np.log(host)
    # Bound retained coefficient paths; SI catalogs may contain many slices.
    if len(solver._picard_coefficients) >= 512:
        solver._picard_coefficients.pop(next(iter(solver._picard_coefficients)))
    solver._picard_coefficients[key] = result
    return result


def _iterate(solver, log_ma, z_acc, z_final, n_integration, n_iterations, coefficients=None):
    """Solve for delta ln(m) with vectorized Picard updates on a fixed redshift path.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.
    log_ma : array_like
        Natural logarithm of subhalo mass in the same fixed mass unit as the host callback;
        only differences of these logarithms enter the rate.
    z_acc : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    z_final : float or array_like, optional
        Target redshift no greater than the accretion redshift. Explicit tables use one
        scalar target; None selects solver.z_min.
    n_integration : int
        Number of nodes on a uniformly spaced redshift integration path, including both
        endpoints.
    n_iterations : int
        Number of nonlinear Picard updates after the initial zero-exponent integral. One
        extra update estimates the remaining iteration error.
    coefficients : tuple, optional
        Precomputed redshift, Phi, zeta and ln(host mass) arrays on this integration path.

    Returns
    -------
    tuple
        Integration grid (n_integration,), delta ln(m) (n_integration,n_mass), and per-mass
        remaining-iteration estimates (n_mass,).

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    z, phi, zeta, log_host = (
        _coefficients(solver, z_acc, z_final, n_integration)
        if coefficients is None
        else coefficients
    )
    log_ma = np.asarray(log_ma)
    initial = cumulative_trapezoid(phi, x=z, initial=0.0)
    delta = np.broadcast_to(initial[:, None], (len(z), log_ma.size)).copy()
    log_ratio = log_ma[None, :] - log_host[:, None]

    def update(current):
        """Integrate the nonlinear tidal rate along the current log-mass trajectory.

        Parameters
        ----------
        current : numpy.ndarray
            Current estimate of delta ln(m), with shape (n_redshift, n_mass).

        Returns
        -------
        numpy.ndarray
            Updated delta ln(m) with shape (n_redshift,n_mass), zero at accretion.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        rhs = phi[:, None] * np.exp(zeta[:, None] * (log_ratio + current))
        return cumulative_trapezoid(rhs, x=z, axis=0, initial=0.0)

    for _ in range(n_iterations):
        delta = update(delta)
    next_delta = update(delta)
    residual = np.max(np.abs(next_delta - delta), axis=0)
    # A conservative local contraction estimate from the integrated rate.
    contraction = np.max(np.abs(zeta)) * np.max(np.abs(next_delta), axis=0)
    estimate = np.full(residual.shape, np.inf)
    valid = contraction < 1.0
    estimate[valid] = residual[valid] / (1 - contraction[valid])
    return z, delta, estimate


def direct_log_mass(solver, ma, za, z, *, rtol=1e-11, atol=1e-12):
    """DOP853 integration of log mass, preserving scalar or full-history shape.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.
    ma : float or array_like
        Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
        queries broadcast over mass and redshift.
    za : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    z : float or array_like
        Dimensionless redshift. The cosmological background requires z > -1; population
        nodes use nonnegative redshifts.
    rtol : float or array_like, optional
        Positive relative or absolute local ODE error tolerance. Absolute tolerance uses the
        integrated state units; it is not a physical-model uncertainty.
    atol : float or array_like, optional
        Positive relative or absolute local ODE error tolerance. Absolute tolerance uses the
        integrated state units; it is not a physical-model uncertainty.

    Returns
    -------
    float or numpy.ndarray
        Positive evolved masses; scalar target retains ma.shape, while a history prepends
        the redshift axis. Repeated target epochs retain their input order.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    mass = _real(ma, "ma")
    redshifts = _real(z, "z")
    za = float(_real(za, "za"))
    if np.any(mass <= 0) or za <= -1 or np.any(redshifts <= -1) or np.any(redshifts > za):
        raise ValueError("Require ma>0 and -1<z<=za.")
    scalar_z = redshifts.ndim == 0
    times = np.atleast_1d(redshifts)
    if times.ndim != 1 or np.any(np.diff(times) > 0):
        raise ValueError("History redshifts must be one-dimensional and non-increasing.")
    if mass.size == 0 or times.size == 0:
        return np.empty(mass.shape if scalar_z else times.shape + mass.shape)
    if times[-1] == za:
        return mass.copy() if scalar_z else np.broadcast_to(mass, times.shape + mass.shape).copy()

    def rhs(redshift, log_mass):
        """Evaluate d ln(m)/dz using the model rate and contemporaneous host mass.

        Parameters
        ----------
        redshift : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.
        log_mass : array_like
            Natural logarithm of subhalo mass in the same fixed mass unit as the host callback;
            only differences of these logarithms enter the rate.

        Returns
        -------
        numpy.ndarray
            Log-mass derivative with the same flattened mass-state shape.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        return solver.Phi(redshift) * np.exp(
            solver.zetaMz(redshift) * (log_mass - np.log(solver.Mzvir(redshift)))
        )

    # Integrate from the declared accretion time even if it is not in t_eval.
    unique, reverse = np.unique(times, return_inverse=True)
    evaluation = unique[::-1]
    include_start = evaluation[0] != za
    grid = np.concatenate(([za], evaluation)) if include_start else evaluation
    values = solve_evolution(rhs, np.log(mass.ravel()), grid, method="DOP853", rtol=rtol, atol=atol)
    if include_start:
        values = values[1:]
    answer = np.exp(values[::-1][reverse]).reshape(times.shape + mass.shape)
    answer = answer[0] if scalar_z else answer
    return float(answer) if answer.ndim == 0 else answer


def _fallback(solver, ma, za, z, reason):
    """Record and warn about a Picard fallback, then integrate the same tidal equation
    directly.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.
    ma : float or array_like
        Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
        queries broadcast over mass and redshift.
    za : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    z : float or array_like
        Dimensionless redshift. The cosmological background requires z > -1; population
        nodes use nonnegative redshifts.
    reason : str
        Physical-domain or iteration-convergence reason recorded for a direct-ODE fallback.

    Returns
    -------
    float or numpy.ndarray
        DOP853 mass solution with the same endpoint/history shape as the request.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """
    _state(solver)
    event = dict(
        reason=reason,
        mass_count=int(np.size(ma)),
        z_acc=float(za),
        z_final=float(np.min(np.atleast_1d(z))),
    )
    solver._picard_events.append(event)
    warnings.warn(
        f"Picard fallback to log-mass DOP853: {reason}; {event['mass_count']} masses.",
        PicardFallbackWarning,
        stacklevel=3,
    )
    return direct_log_mass(solver, ma, za, z)


class PicardMassLossTable:
    """Tabulate tidal delta ln(m) over accretion epoch and initial host-mass ratio.

    The convergence estimate covers the nonlinear iterations. Interpolation
    and quadrature errors are separately validated by the saved ODE campaign.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.
    z_final : float or array_like, optional
        Target redshift no greater than the accretion redshift. Explicit tables use one
        scalar target; None selects solver.z_min.
    z_acc_max : float, optional
        Largest tabulated accretion redshift; None selects solver.z_max.
    n_z_acc : int
        Number of endpoint-table accretion-redshift nodes, at least two (four for cubic
        interpolation).
    n_log_ratio : int
        Number of endpoint-table log-mass-ratio nodes, at least two (four for cubic
        interpolation).
    log10_ratio_min : float
        Lower or upper endpoint-table bound in log10[ma/M_host(z_acc)].
    log10_ratio_max : float
        Lower or upper endpoint-table bound in log10[ma/M_host(z_acc)].
    n_integration : int
        Number of nodes on a uniformly spaced redshift integration path, including both
        endpoints.
    n_iterations : int
        Number of nonlinear Picard updates after the initial zero-exponent integral. One
        extra update estimates the remaining iteration error.
    convergence_tolerance : float
        Positive threshold for the estimated remaining change in ln(m). It does not bound
        quadrature or interpolation error.
    interpolation : str, optional
        Interpolation rule supported by the owning table; see Notes for its coordinates and
        permitted choices.

    Notes
    -----
    Call delta_log_mass or mass only within the stored coordinate range and check converged
    for the iteration estimate. These methods do not automatically fall back to an ODE. Use
    solver.endpoint for guarded reuse and fallback. A changed physical-state key makes an
    existing table stale. Cubic interpolation uses at least four nodes per axis; otherwise
    linear is available.
    """

    def __init__(
        self,
        solver,
        z_final=None,
        z_acc_max=None,
        n_z_acc=48,
        n_log_ratio=32,
        log10_ratio_min=-24.0,
        log10_ratio_max=3.0,
        n_integration=129,
        n_iterations=3,
        convergence_tolerance=8e-4,
        interpolation=None,
    ):
        """Precompute endpoint mass loss and an iteration estimate on the requested grid.

        Parameters
        ----------
        solver : PicardMassLossSolver
            Model-supplied rate, exponent and host-mass callbacks with an identity that changes
            when their physical inputs change.
        z_final : float or array_like, optional
            Target redshift no greater than the accretion redshift. Explicit tables use one
            scalar target; None selects solver.z_min.
        z_acc_max : float, optional
            Largest tabulated accretion redshift; None selects solver.z_max.
        n_z_acc : int
            Number of endpoint-table accretion-redshift nodes, at least two (four for cubic
            interpolation).
        n_log_ratio : int
            Number of endpoint-table log-mass-ratio nodes, at least two (four for cubic
            interpolation).
        log10_ratio_min : float
            Lower or upper endpoint-table bound in log10[ma/M_host(z_acc)].
        log10_ratio_max : float
            Lower or upper endpoint-table bound in log10[ma/M_host(z_acc)].
        n_integration : int
            Number of nodes on a uniformly spaced redshift integration path, including both
            endpoints.
        n_iterations : int
            Number of nonlinear Picard updates after the initial zero-exponent integral. One
            extra update estimates the remaining iteration error.
        convergence_tolerance : float
            Positive threshold for the estimated remaining change in ln(m). It does not bound
            quadrature or interpolation error.
        interpolation : str, optional
            Interpolation rule supported by the owning table; see Notes for its coordinates and
            permitted choices.

        Notes
        -----
        None interpolation selects cubic when both axes have at least four nodes, otherwise
        linear. Bounds describe table coordinates; they are not evidence that the physical model
        was validated there.
        """
        self.solver = solver
        self.z_final = float(_real(solver.z_min if z_final is None else z_final, "z_final"))
        self.z_acc_max = float(_real(solver.z_max if z_acc_max is None else z_acc_max, "z_acc_max"))
        if self.z_final <= -1 or self.z_acc_max < self.z_final:
            raise ValueError("Require -1 < z_final <= z_acc_max.")
        self.n_z_acc = _integer(n_z_acc, "n_z_acc", 2)
        self.n_log_ratio = _integer(n_log_ratio, "n_log_ratio", 2)
        self.n_integration = _integer(n_integration, "n_integration", 2)
        self.n_iterations = _integer(n_iterations, "n_iterations", 1)
        self.log10_ratio_min = float(_real(log10_ratio_min, "log10_ratio_min"))
        self.log10_ratio_max = float(_real(log10_ratio_max, "log10_ratio_max"))
        if self.log10_ratio_max <= self.log10_ratio_min:
            raise ValueError("log10_ratio_max must exceed log10_ratio_min.")
        self.convergence_tolerance = float(_real(convergence_tolerance, "convergence_tolerance"))
        if self.convergence_tolerance <= 0:
            raise ValueError("convergence_tolerance must be positive.")
        if interpolation is None:
            interpolation = "cubic" if min(self.n_z_acc, self.n_log_ratio) >= 4 else "linear"
        if interpolation not in ("linear", "cubic"):
            raise ValueError("interpolation must be linear or cubic.")
        if interpolation == "cubic" and min(self.n_z_acc, self.n_log_ratio) < 4:
            raise ValueError("Cubic tables need at least four nodes on each axis.")
        self.interpolation = interpolation
        self.z_acc_grid = np.linspace(self.z_final, self.z_acc_max, self.n_z_acc)
        self.log_ratio_grid = np.linspace(
            self.log10_ratio_min, self.log10_ratio_max, self.n_log_ratio
        )
        self._ln_ratio_grid = self.log_ratio_grid * np.log(10.0)
        self.delta_ln_mass = np.zeros((self.n_z_acc, self.n_log_ratio))
        self.iteration_error = np.zeros_like(self.delta_ln_mass)
        self._physics_key = _state(solver)
        start = time.perf_counter()
        # Evaluate each host/background array once for the complete table.
        # Flattening retains the one-dimensional provider interface.
        paths = np.linspace(self.z_acc_grid, self.z_final, self.n_integration, axis=-1)
        flat = paths.ravel()
        phi = np.broadcast_to(np.asarray(solver.Phi(flat)), flat.shape).reshape(paths.shape)
        zeta = np.broadcast_to(np.asarray(solver.zetaMz(flat)), flat.shape).reshape(paths.shape)
        host = np.broadcast_to(np.asarray(solver.Mzvir(flat)), flat.shape).reshape(paths.shape)
        if (
            np.any(~np.isfinite(host))
            or np.any(host <= 0)
            or np.any(~np.isfinite(phi))
            or np.any(phi < 0)
            or np.any(~np.isfinite(zeta))
        ):
            raise ValueError(
                "Tidal table requires finite, positive host mass and nonnegative rate."
            )
        log_host = np.log(host)
        for i, za in enumerate(self.z_acc_grid):
            if za == self.z_final:
                continue
            log_ma = log_host[i, 0] + self._ln_ratio_grid
            _, history, error = _iterate(
                solver,
                log_ma,
                za,
                self.z_final,
                self.n_integration,
                self.n_iterations,
                (paths[i], phi[i], zeta[i], log_host[i]),
            )
            self.delta_ln_mass[i] = history[-1]
            self.iteration_error[i] = error
        self.build_seconds = time.perf_counter() - start
        if self.z_final == self.z_acc_max:
            self._interpolator = None
        elif interpolation == "linear":
            self._interpolator = RegularGridInterpolator(
                (self.z_acc_grid, self._ln_ratio_grid), self.delta_ln_mass, bounds_error=True
            )
        else:
            self._interpolator = RectBivariateSpline(
                self.z_acc_grid, self._ln_ratio_grid, self.delta_ln_mass, kx=3, ky=3, s=0
            )

    def _points(self, ma, z_acc):
        """Validate an endpoint query and map it to (z_acc, ln[ma/M_host]) coordinates.

        Parameters
        ----------
        ma : float or array_like
            Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
            queries broadcast over mass and redshift.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.

        Returns
        -------
        tuple
            Broadcast accretion masses and a flattened coordinate array of shape (mass.size,2).
            Stale or out-of-range tables raise ValueError.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        if physics_key(self.solver) != getattr(self, "_physics_key", None):
            raise ValueError("Picard table is stale after a host or background change.")
        mass, za = np.broadcast_arrays(_real(ma, "ma"), _real(z_acc, "z_acc"))
        if np.any(mass <= 0):
            raise ValueError("Accretion masses must be positive.")
        if np.any(za < self.z_final) or np.any(za > self.z_acc_max):
            raise ValueError("z_acc lies outside the precomputed redshift range.")
        if mass.size == 0:
            return mass, np.empty((0, 2))
        host = np.asarray(self.solver.Mzvir(za))
        ratio = np.log(mass / host)
        if (
            np.any(~np.isfinite(ratio))
            or np.any(ratio < self._ln_ratio_grid[0])
            or np.any(ratio > self._ln_ratio_grid[-1])
        ):
            raise ValueError("ma/Mvir(z_acc) lies outside the precomputed mass-ratio range.")
        return mass, np.column_stack((za.ravel(), ratio.ravel()))

    def delta_log_mass(self, ma, z_acc):
        """Interpolate ln[m(z_final)/ma] without extrapolating the endpoint table.

        Parameters
        ----------
        ma : float or array_like
            Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
            queries broadcast over mass and redshift.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.

        Returns
        -------
        float or numpy.ndarray
            Dimensionless log-mass change, with broadcast mass/redshift shape.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        mass, points = self._points(ma, z_acc)
        if self._interpolator is None or mass.size == 0:
            values = np.zeros_like(mass)
        else:
            if self.interpolation == "linear":
                values = self._interpolator(points).reshape(mass.shape)
            else:
                values = self._interpolator.ev(points[:, 0], points[:, 1]).reshape(mass.shape)
        return float(values) if values.ndim == 0 else values

    def mass(self, ma, z_acc):
        """Convert the interpolated log-mass change to the final bound mass.

        Parameters
        ----------
        ma : float or array_like
            Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
            queries broadcast over mass and redshift.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.

        Returns
        -------
        float or numpy.ndarray
            ma*exp(delta ln(m)), with broadcast shape and the input mass unit.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        mass = np.asarray(ma)
        values = mass * np.exp(self.delta_log_mass(ma, z_acc))
        return float(values) if values.ndim == 0 else values

    def converged(self, ma, z_acc):
        """Check the largest stored iteration estimate over the interpolation stencil.

        Parameters
        ----------
        ma : float or array_like
            Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
            queries broadcast over mass and redshift.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.

        Returns
        -------
        numpy.ndarray of bool
            True where the stencil estimate is at most convergence_tolerance; this is not a
            total-accuracy certificate.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        mass, points = self._points(ma, z_acc)
        if mass.size == 0 or self._interpolator is None:
            return np.ones(mass.shape, dtype=bool)
        i = np.clip(np.searchsorted(self.z_acc_grid, points[:, 0]) - 1, 0, self.n_z_acc - 2)
        j = np.clip(np.searchsorted(self._ln_ratio_grid, points[:, 1]) - 1, 0, self.n_log_ratio - 2)
        offsets = (-1, 0, 1, 2) if self.interpolation == "cubic" else (0, 1)
        error = np.maximum.reduce(
            [
                self.iteration_error[
                    np.clip(i + di, 0, self.n_z_acc - 1), np.clip(j + dj, 0, self.n_log_ratio - 1)
                ]
                for di in offsets
                for dj in offsets
            ]
        )
        return (error <= self.convergence_tolerance).reshape(mass.shape)


def endpoint_mass(
    solver,
    ma,
    za,
    z,
    *,
    n_z_acc=48,
    n_log_ratio=32,
    log10_ratio_min=-24.0,
    log10_ratio_max=3.0,
    n_integration=129,
    n_iterations=3,
    convergence_tolerance=8e-4,
    interpolation=None,
):
    """Evaluate guarded tidal endpoints, reusing tables with identical physical and numerical
    settings.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.
    ma : float or array_like
        Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
        queries broadcast over mass and redshift.
    za : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    z : float or array_like
        Dimensionless redshift. The cosmological background requires z > -1; population
        nodes use nonnegative redshifts.
    n_z_acc : int
        Number of endpoint-table accretion-redshift nodes, at least two (four for cubic
        interpolation).
    n_log_ratio : int
        Number of endpoint-table log-mass-ratio nodes, at least two (four for cubic
        interpolation).
    log10_ratio_min : float
        Lower or upper endpoint-table bound in log10[ma/M_host(z_acc)].
    log10_ratio_max : float
        Lower or upper endpoint-table bound in log10[ma/M_host(z_acc)].
    n_integration : int
        Number of nodes on a uniformly spaced redshift integration path, including both
        endpoints.
    n_iterations : int
        Number of nonlinear Picard updates after the initial zero-exponent integral. One
        extra update estimates the remaining iteration error.
    convergence_tolerance : float
        Positive threshold for the estimated remaining change in ln(m). It does not bound
        quadrature or interpolation error.
    interpolation : str, optional
        Interpolation rule supported by the owning table; see Notes for its coordinates and
        permitted choices.

    Returns
    -------
    float or numpy.ndarray
        Final masses with the broadcast shape of ma, za and z; no evolution is applied when
        z equals za.

    Notes
    -----
    Out-of-domain or unconverged points use counted, warned DOP853 fallbacks. Requested
    masses can explicitly expand a cached table range; the actual settings are exposed by
    diagnostics().
    """
    n_z_acc = _integer(n_z_acc, "n_z_acc", 2)
    n_log_ratio = _integer(n_log_ratio, "n_log_ratio", 2)
    n_integration = _integer(n_integration, "n_integration", 2)
    n_iterations = _integer(n_iterations, "n_iterations", 1)
    log10_ratio_min = float(_real(log10_ratio_min, "log10_ratio_min"))
    log10_ratio_max = float(_real(log10_ratio_max, "log10_ratio_max"))
    if log10_ratio_max <= log10_ratio_min:
        raise ValueError("log10_ratio_max must exceed log10_ratio_min.")
    convergence_tolerance = float(_real(convergence_tolerance, "convergence_tolerance"))
    if convergence_tolerance <= 0:
        raise ValueError("convergence_tolerance must be positive.")
    if interpolation is None:
        interpolation = "cubic" if min(n_z_acc, n_log_ratio) >= 4 else "linear"
    if interpolation not in ("linear", "cubic"):
        raise ValueError("interpolation must be linear or cubic.")
    if interpolation == "cubic" and min(n_z_acc, n_log_ratio) < 4:
        raise ValueError("Cubic tables need at least four nodes on each axis.")
    mass, acc, final = np.broadcast_arrays(_real(ma, "ma"), _real(za, "za"), _real(z, "z"))
    if np.any(mass <= 0) or np.any(acc <= -1) or np.any(final <= -1) or np.any(final > acc):
        raise ValueError("Require ma>0 and -1<z<=za.")
    _state(solver)
    output = mass.copy()
    if mass.size == 0:
        return output
    # Group distinct final redshifts while preserving the usual scalar-target
    # table reuse. SI uses history_mass instead of creating endpoint tables.
    for target in np.unique(final):
        selected = (final == target) & (acc > target)
        if not np.any(selected):
            continue
        m, a = mass[selected], acc[selected]
        host = np.asarray(solver.Mzvir(a))
        if np.any(~np.isfinite(host)) or np.any(host <= 0):
            raise ValueError("Tidal stripping requires finite positive host mass.")
        ratio = np.log10(m / host)
        outside = solver.outside_domain(a, target, ratio)
        selected_output = np.empty_like(m)
        for start in np.unique(a[outside]):
            use = outside & (a == start)
            selected_output[use] = _fallback(
                solver, m[use], start, target, "outside validated redshift or mass-ratio domain"
            )
        inside = ~outside
        if np.any(inside):
            # Required coordinates are included even when custom table limits
            # are narrower. Expansion is explicit in the cached table metadata.
            lo = min(float(log10_ratio_min), float(np.floor(np.min(ratio[inside])))) - 1e-10
            hi = max(float(log10_ratio_max), float(np.ceil(np.max(ratio[inside])))) + 1e-10
            upper = (
                solver.z_max
                if solver.domain is None
                else min(solver.z_max, solver.domain.redshift[1])
            )
            zmax = max(float(upper) + 32 * np.finfo(float).eps, float(np.max(a[inside])))
            key = (
                float(target),
                zmax,
                n_z_acc,
                n_log_ratio,
                lo,
                hi,
                n_integration,
                n_iterations,
                convergence_tolerance,
                interpolation,
            )
            table = solver._picard_tables.get(key)
            if table is None:
                table = PicardMassLossTable(
                    solver,
                    target,
                    zmax,
                    n_z_acc,
                    n_log_ratio,
                    lo,
                    hi,
                    n_integration,
                    n_iterations,
                    convergence_tolerance,
                    interpolation,
                )
                solver._picard_tables[key] = table
            ok = table.converged(m[inside], a[inside])
            values = table.mass(m[inside], a[inside])
            for start in np.unique(a[inside][~ok]):
                use = (~ok) & (a[inside] == start)
                values[use] = _fallback(
                    solver,
                    m[inside][use],
                    start,
                    target,
                    "Picard iteration convergence estimate exceeded tolerance",
                )
            selected_output[inside] = values
        output[selected] = selected_output
    return float(output) if output.ndim == 0 else output


def history_mass(
    solver,
    ma,
    za,
    z,
    *,
    n_integration=100,
    n_iterations=3,
    convergence_tolerance=8e-4,
    n_mass_nodes=32,
):
    """Evaluate complete tidal histories from one accretion epoch.

    Parameters
    ----------
    solver : PicardMassLossSolver
        Model-supplied rate, exponent and host-mass callbacks with an identity that changes
        when their physical inputs change.
    ma : float or array_like
        Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
        queries broadcast over mass and redshift.
    za : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    z : float or array_like
        Dimensionless redshift. The cosmological background requires z > -1; population
        nodes use nonnegative redshifts.
    n_integration : int
        Number of nodes on a uniformly spaced redshift integration path, including both
        endpoints.
    n_iterations : int
        Number of nonlinear Picard updates after the initial zero-exponent integral. One
        extra update estimates the remaining iteration error.
    convergence_tolerance : float
        Positive threshold for the estimated remaining change in ln(m). It does not bound
        quadrature or interpolation error.
    n_mass_nodes : int
        Minimum history-interpolation mass-node count. The implementation adds nodes when
        necessary to keep spacing at most 0.25 dex.

    Returns
    -------
    float or numpy.ndarray
        Shape (n_redshift,*ma.shape) for a history or ma.shape for a scalar target, in the
        input mass unit.

    Notes
    -----
    Redshifts must be non-increasing; repeated epochs and omitted accretion endpoints are
    supported. Linear redshift/mass interpolation is used when needed. Domain violations or
    failed iteration estimates trigger warned DOP853 fallback; refine the grid separately to
    assess interpolation and quadrature error.
    """
    mass = _real(ma, "ma")
    times = _real(z, "z")
    za = float(_real(za, "za"))
    scalar_z = times.ndim == 0
    times = np.atleast_1d(times)
    count = _integer(n_integration, "n_integration", 2)
    iterations = _integer(n_iterations, "n_iterations", 1)
    mass_nodes = _integer(n_mass_nodes, "n_mass_nodes", 2)
    if np.any(mass <= 0) or za <= -1 or np.any(times <= -1) or np.any(times > za):
        raise ValueError("Require ma>0 and -1<z<=za.")
    if times.ndim != 1 or np.any(np.diff(times) > 0):
        raise ValueError("History redshifts must be one-dimensional and non-increasing.")
    if not np.isfinite(convergence_tolerance) or convergence_tolerance <= 0:
        raise ValueError("convergence_tolerance must be finite and positive.")
    _state(solver)
    if mass.size == 0 or times.size == 0:
        return np.empty(mass.shape if scalar_z else times.shape + mass.shape)
    if times[-1] == za:
        return mass.copy() if scalar_z else np.broadcast_to(mass, times.shape + mass.shape).copy()
    ratio = np.log10(mass / solver.Mzvir(za))
    if np.any(solver.outside_domain(za, times[-1], ratio)):
        return _fallback(solver, mass, za, z, "outside validated redshift or mass-ratio domain")
    log_mass = np.log(mass.ravel())
    # Wider requests receive enough mass nodes to keep the interpolation
    # spacing within the validated quarter-dex bound.
    mass_nodes = max(mass_nodes, int(np.ceil(np.ptp(log_mass) / (np.log(10.0) * 0.25))) + 1)
    use_mass_grid = log_mass.size > mass_nodes and np.ptp(log_mass) > 0
    nodes = (
        np.linspace(np.min(log_mass), np.max(log_mass), mass_nodes) if use_mass_grid else log_mass
    )
    history_key = (
        count,
        iterations,
        float(convergence_tolerance),
        int(nodes.size) if use_mass_grid else 0,
    )
    solver._picard_histories[history_key] = solver._picard_histories.get(history_key, 0) + 1
    grid, delta, error = _iterate(solver, nodes, za, float(times[-1]), count, iterations)
    # Vectorized linear interpolation of the entire requested trajectory.
    if np.array_equal(grid, times):
        requested = delta
    else:
        positions = (za - times) / (za - times[-1]) * (count - 1)
        left = np.minimum(positions.astype(int), count - 2)
        fraction = positions - left
        requested = delta[left] * (1 - fraction[:, None]) + delta[left + 1] * fraction[:, None]
    if use_mass_grid:
        coordinate = (log_mass - nodes[0]) / (nodes[-1] - nodes[0]) * (mass_nodes - 1)
        left = np.minimum(coordinate.astype(int), mass_nodes - 2)
        fraction = coordinate - left
        requested = (
            requested[:, left] * (1 - fraction[None, :])
            + requested[:, left + 1] * fraction[None, :]
        )
        error = np.maximum(error[left], error[left + 1])
    output = (mass.ravel()[None, :] * np.exp(requested)).reshape(times.shape + mass.shape)
    bad = error > convergence_tolerance
    if np.any(bad):
        flattened = output.reshape(len(times), -1)
        flattened[:, bad] = _fallback(
            solver,
            mass.ravel()[bad],
            za,
            times,
            "Picard iteration convergence estimate exceeded tolerance",
        )
    output = output[0] if scalar_z else output
    return float(output) if output.ndim == 0 else output


@dataclass(frozen=True)
class PicardDomain:
    """Caller-validated domain; outside points use direct log-mass integration.

    ``redshift`` and ``log10_mass_ratio`` are inclusive lower/upper bounds.
    Domain validity is a scientific responsibility of the supplying model.

    Parameters
    ----------
    redshift : tuple of float
        Inclusive lower and upper validated redshift bounds, with -1 < lower < upper.
    log10_mass_ratio : tuple of float
        Inclusive lower and upper bounds for log10[ma/M_host(z_acc)] in the caller-validated
        domain.

    Notes
    -----
    The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
    integration from accretion to a lower redshift decreases mass. All masses must share the
    caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
    errors.
    """

    redshift: tuple[float, float]
    log10_mass_ratio: tuple[float, float]

    def __post_init__(self):
        """Require finite increasing physical-domain bounds before accepting a validation claim.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        for name in ("redshift", "log10_mass_ratio"):
            value = _real(getattr(self, name), name)
            if value.shape != (2,) or value[0] >= value[1]:
                raise ValueError(f"{name} must contain two increasing finite bounds.")
            object.__setattr__(self, name, tuple(float(x) for x in value))
        if self.redshift[0] <= -1:
            raise ValueError("redshift bounds must exceed -1.")


class PicardMassLossSolver:
    """Integrate d ln(m)/dz = rate(z) * (m/reference_mass(z))**exponent(z).

    Callbacks accept scalar or one-dimensional redshifts and return compatible
    finite values; rate is nonnegative and reference mass positive. Mass units
    must agree, but ITAMAE does not choose them or the model's coefficients.
    ``state_key`` must return an immutable equality-comparable value covering
    all callback inputs (including cosmology, units and mutable arrays). A
    changed key invalidates caches and rejects previously retained tables.

    Endpoint and full-history queries share the same nonlinear Picard updates.
    The iteration check does not bound quadrature/interpolation error; callers
    must verify those against a direct reference and grid refinement.

    Parameters
    ----------
    rate : callable
        Nonnegative Phi(z) in d ln(m)/dz; accepts scalar or one-dimensional redshifts.
    exponent : callable
        Dimensionless zeta(z) multiplying ln(m/M_host) in the tidal-loss law.
    reference_mass : callable
        Positive host/reference mass M_host(z), in the same units as accretion masses.
    z_min : float
        Lower or upper redshift extent declared by the caller; require -1 < z_min <= z_max.
    z_max : float
        Lower or upper redshift extent declared by the caller; require -1 < z_min <= z_max.
    state_key : callable
        Zero-argument callback returning an immutable equality-comparable identity covering
        cosmology, units and every mutable coefficient input.
    domain : PicardDomain, optional
        Caller-validated redshift and initial mass-ratio bounds. None declares no such
        restriction; it does not establish validation.

    Notes
    -----
    The model supplies Phi, zeta and M_host in
    dln(m)/dz=Phi(z)*exp(zeta(z)*ln(m/M_host(z))). Integrating toward lower redshift
    therefore decreases mass when Phi is nonnegative. Use endpoint or history for guarded
    Picard evaluation, direct for a DOP853 reference, and diagnostics for actual settings
    and fallback counts. The caller state_key must change whenever any physical input
    changes. Iteration estimates alone do not bound interpolation or quadrature error.
    """

    def __init__(self, rate, exponent, reference_mass, *, z_min, z_max, state_key, domain=None):
        """Bind a model-owned tidal-loss equation to reusable numerical integration state.

        Parameters
        ----------
        rate : callable
            Nonnegative Phi(z) in d ln(m)/dz; accepts scalar or one-dimensional redshifts.
        exponent : callable
            Dimensionless zeta(z) multiplying ln(m/M_host) in the tidal-loss law.
        reference_mass : callable
            Positive host/reference mass M_host(z), in the same units as accretion masses.
        z_min : float
            Lower or upper redshift extent declared by the caller; require -1 < z_min <= z_max.
        z_max : float
            Lower or upper redshift extent declared by the caller; require -1 < z_min <= z_max.
        state_key : callable
            Zero-argument callback returning an immutable equality-comparable identity covering
            cosmology, units and every mutable coefficient input.
        domain : PicardDomain, optional
            Caller-validated redshift and initial mass-ratio bounds. None declares no such
            restriction; it does not establish validation.

        Raises
        ------
        ValueError
            If redshift bounds are non-finite, reversed, or at/below z=-1.
        TypeError
            If the physical callbacks or state_key are not callable.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        for callback in (rate, exponent, reference_mass, state_key):
            if not callable(callback):
                raise TypeError("Coefficients and state_key must be callable.")
        self.z_min = float(_real(z_min, "z_min"))
        self.z_max = float(_real(z_max, "z_max"))
        if self.z_min <= -1 or self.z_max < self.z_min:
            raise ValueError("Require -1 < z_min <= z_max.")
        if domain is not None and not isinstance(domain, PicardDomain):
            raise TypeError("domain must be a PicardDomain or None.")
        self.domain = domain
        self._rate = rate
        self._exponent = exponent
        self._reference_mass = reference_mass
        self._input_key = state_key
        self._picard_events = []
        _state(self)

    def state_key(self):
        """Combine callback identity, declared redshift extent and validated domain for cache
        invalidation.

        Returns
        -------
        tuple
            State identity including callback objects and the caller physical-input key; not a
            durable source hash.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        return (
            self._input_key(),
            self.z_min,
            self.z_max,
            self.domain,
            id(self._rate),
            id(self._exponent),
            id(self._reference_mass),
        )

    def Phi(self, z):
        """Evaluate the nonnegative tidal rate per unit redshift.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Phi(z), broadcast to the redshift shape.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        value = np.broadcast_to(_real(self._rate(z), "rate"), np.shape(z))
        if np.any(value < 0):
            raise ValueError("rate must be nonnegative.")
        return value

    def zetaMz(self, z):
        """Evaluate the dimensionless host-mass-ratio exponent.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            zeta(z), broadcast to the redshift shape.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        return np.broadcast_to(_real(self._exponent(z), "exponent"), np.shape(z))

    def Mzvir(self, z):
        """Evaluate and cache a copy of the positive reference host mass.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Reference mass M_host(z), in the caller mass unit and redshift shape.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        _state(self)
        redshift = np.asarray(z, dtype=float)
        key = (redshift.shape, redshift.tobytes())
        cache = self._picard_host_cache
        if key not in cache:
            value = np.broadcast_to(
                _real(self._reference_mass(z), "reference_mass"), redshift.shape
            )
            if np.any(value <= 0):
                raise ValueError("reference_mass must be positive.")
            if len(cache) >= 4096:
                cache.pop(next(iter(cache)))
            cache[key] = value.copy()
        result = cache[key].copy()
        return float(result) if result.ndim == 0 else result

    def outside_domain(self, z_acc, z_final, log10_ratio):
        """Mark queries outside the model-supplied validation rectangle.

        Parameters
        ----------
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.
        z_final : float or array_like, optional
            Target redshift no greater than the accretion redshift. Explicit tables use one
            scalar target; None selects solver.z_min.
        log10_ratio : array_like
            Initial log10[ma/M_host(z_acc)] used for the caller-supplied domain decision.

        Returns
        -------
        numpy.ndarray of bool
            Broadcast mask for unsupported accretion/target redshifts or initial mass ratios;
            all false if no domain was supplied.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        a, z, ratio = np.broadcast_arrays(z_acc, z_final, log10_ratio)
        if self.domain is None:
            return np.zeros(a.shape, dtype=bool)
        lo_z, hi_z = self.domain.redshift
        lo_m, hi_m = self.domain.log10_mass_ratio
        epsilon = 32 * np.finfo(float).eps
        return (a > hi_z + epsilon) | (z < lo_z) | (ratio < lo_m) | (ratio > hi_m + epsilon)

    def endpoint(self, ma, z_acc, z_final, **options):
        """Return final masses, with broadcasting across mass and redshift.

        Parameters
        ----------
        ma : float or array_like
            Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
            queries broadcast over mass and redshift.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.
        z_final : float or array_like, optional
            Target redshift no greater than the accretion redshift. Explicit tables use one
            scalar target; None selects solver.z_min.
        **options : dict
            Keyword numerical settings accepted by endpoint_mass, including grid counts,
            iteration count, interpolation and convergence_tolerance.

        Returns
        -------
        float or numpy.ndarray
            Guarded endpoint mass with broadcast mass/redshift shape.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        return endpoint_mass(self, ma, z_acc, z_final, **options)

    def history(self, ma, z_acc, redshifts, **options):
        """Return shape (n_redshifts, *mass_shape); scalar targets omit that axis.

        Parameters
        ----------
        ma : float or array_like
            Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
            queries broadcast over mass and redshift.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.
        redshifts : float or array_like
            Scalar target or one-dimensional non-increasing redshift history, at or below the
            scalar accretion redshift.
        **options : dict
            Keyword settings accepted by history_mass: n_integration, n_iterations, n_mass_nodes
            and convergence_tolerance.

        Returns
        -------
        float or numpy.ndarray
            History mass with a leading redshift axis, or ma.shape for a scalar target.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        return history_mass(self, ma, z_acc, redshifts, **options)

    def direct(self, ma, z_acc, redshifts, **options):
        """High-accuracy DOP853 reference using ITAMAE's common ODE controller.

        Parameters
        ----------
        ma : float or array_like
            Positive subhalo mass at accretion, in the same units as reference_mass; endpoint
            queries broadcast over mass and redshift.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.
        redshifts : float or array_like
            Scalar target or one-dimensional non-increasing redshift history, at or below the
            scalar accretion redshift.
        **options : dict
            Optional rtol and atol for direct_log_mass; these control local error in ln(m).

        Returns
        -------
        float or numpy.ndarray
            Direct log-mass ODE solution with the requested endpoint/history shape.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        return direct_log_mass(self, ma, z_acc, redshifts, **options)

    def table(self, **options):
        """Build an explicit table; use endpoint() for guarded cached queries.

        Parameters
        ----------
        **options : dict
            Keyword constructor settings for PicardMassLossTable, including target redshift and
            table coordinates.

        Returns
        -------
        PicardMassLossTable
            Explicit endpoint table bound to the present physical state.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        return PicardMassLossTable(self, **options)

    def diagnostics(self):
        """Return JSON-compatible settings and actual fallback counts.

        Returns
        -------
        dict
            JSON-compatible backend, actual cached-table/history settings, domain and fallback
            events/counts for the current physical state.

        Notes
        -----
        The equation is d ln(m)/dz = Phi(z) * [m/M_host(z)]**zeta(z). With nonnegative Phi,
        integration from accretion to a lower redshift decreases mass. All masses must share the
        caller unit. The iteration estimate excludes trapezoidal quadrature and interpolation
        errors.
        """
        _state(self)
        names = (
            "z_final",
            "z_acc_max",
            "n_z_acc",
            "n_log_ratio",
            "log10_ratio_min",
            "log10_ratio_max",
            "n_integration",
            "n_iterations",
            "convergence_tolerance",
            "interpolation",
        )
        tables = [
            {name: getattr(table, name) for name in names} for table in self._picard_tables.values()
        ]
        return {
            "backend": "itamae:picard-log-mass:v1",
            "validated_domain": None
            if self.domain is None
            else {
                "redshift": list(self.domain.redshift),
                "log10_mass_ratio": list(self.domain.log10_mass_ratio),
            },
            "tables": tables,
            "histories": [
                dict(
                    n_integration=key[0],
                    n_iterations=key[1],
                    convergence_tolerance=key[2],
                    n_mass_nodes_used=key[3],
                    calls=count,
                )
                for key, count in self._picard_histories.items()
            ],
            "fallback_calls": len(self._picard_events),
            "fallback_mass_count": sum(event["mass_count"] for event in self._picard_events),
            "fallback_events": [dict(event) for event in self._picard_events],
            "fallback_method": "DOP853",
            "fallback_rtol": 1e-11,
            "fallback_atol": 1e-12,
        }
