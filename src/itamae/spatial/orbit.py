"""Generic spherical-orbit kernels for model-supplied potentials.

The functions in this module implement numerical mechanisms only. SASHIMI
variants remain responsible for selecting a host potential, infall
distribution, dynamical-friction law, and disruption prescription.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from .radial import RadialMeasure


def _validate_orbit(energy: float, angular_momentum: float) -> tuple[float, float]:
    """Return finite scalar orbit invariants.

    Parameters
    ----------
    energy : float or array_like
        Specific orbital energy, in the same velocity-squared units as the supplied
        potential.
    angular_momentum : float or array_like
        Nonnegative specific angular momentum, in length times velocity units consistent
        with the potential.

    Returns
    -------
    tuple of float
        Finite specific energy and nonnegative specific angular momentum.

    Notes
    -----
    These kernels assume a time-independent spherical potential and a bound orbit with two
    resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
    km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
    Gyr.
    """
    energy = float(energy)
    angular_momentum = float(angular_momentum)
    if not np.isfinite(energy):
        raise ValueError("Orbital energy must be finite.")
    if not np.isfinite(angular_momentum) or angular_momentum < 0.0:
        raise ValueError("Angular momentum must be finite and nonnegative.")
    return energy, angular_momentum


def radial_velocity_squared(
    radius: Any,
    potential: Callable[[Any], Any],
    energy: float,
    angular_momentum: float,
) -> np.ndarray:
    r"""Return :math:`v_r^2=2(E-\Phi)-L^2/r^2` in consistent units.

    Parameters
    ----------
    radius : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    potential : callable
        Spherical potential Phi(r) with the same energy zero and velocity-squared units as
        energy; accepts the supplied radius shape.
    energy : float or array_like
        Specific orbital energy, in the same velocity-squared units as the supplied
        potential.
    angular_momentum : float or array_like
        Nonnegative specific angular momentum, in length times velocity units consistent
        with the potential.

    Returns
    -------
    numpy.ndarray
        v_r**2 in the potential velocity-squared unit; negative values mark inaccessible
        radii.

    Notes
    -----
    These kernels assume a time-independent spherical potential and a bound orbit with two
    resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
    km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
    Gyr.
    """
    energy, angular_momentum = _validate_orbit(energy, angular_momentum)
    radius = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("Radii must be finite and positive.")
    potential_value = np.asarray(potential(radius), dtype=float)
    try:
        potential_value = np.broadcast_to(potential_value, radius.shape)
    except ValueError as exc:
        raise ValueError("Potential output must broadcast to the radius shape.") from exc
    if not np.all(np.isfinite(potential_value)):
        raise ValueError("Potential returned non-finite values.")
    return 2.0 * (energy - potential_value) - angular_momentum**2 / radius**2


def turning_points(
    potential: Callable[[Any], Any],
    energy: float,
    angular_momentum: float,
    r_min: float,
    r_max: float,
    *,
    n_scan: int = 2049,
) -> tuple[float, float]:
    """Bracket and solve the pericenter and apocenter of a bound orbit.

    Parameters
    ----------
    potential : callable
        Spherical potential Phi(r) with the same energy zero and velocity-squared units as
        energy; accepts the supplied radius shape.
    energy : float or array_like
        Specific orbital energy, in the same velocity-squared units as the supplied
        potential.
    angular_momentum : float or array_like
        Nonnegative specific angular momentum, in length times velocity units consistent
        with the potential.
    r_min : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    r_max : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    n_scan : int
        Number of logarithmic scan points used to bracket turning points. Very narrow
        allowed intervals require a refined scan.

    Returns
    -------
    tuple of float
        Pericenter and apocenter in the input radius unit; unresolved two-root intervals
        raise ValueError.

    Notes
    -----
    These kernels assume a time-independent spherical potential and a bound orbit with two
    resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
    km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
    Gyr.
    """
    energy, angular_momentum = _validate_orbit(energy, angular_momentum)
    r_min = float(r_min)
    r_max = float(r_max)
    if not np.isfinite(r_min) or not np.isfinite(r_max) or r_min <= 0.0 or r_min >= r_max:
        raise ValueError("Turning-point bounds must satisfy 0 < r_min < r_max.")
    if isinstance(n_scan, bool) or not isinstance(n_scan, (int, np.integer)):
        raise TypeError("n_scan must be an integer.")
    if n_scan < 3:
        raise ValueError("n_scan must be at least three.")

    grid = np.geomspace(r_min, r_max, n_scan)
    velocity_squared = radial_velocity_squared(grid, potential, energy, angular_momentum)
    roots: list[float] = []
    for left, right, f_left, f_right in zip(
        grid[:-1],
        grid[1:],
        velocity_squared[:-1],
        velocity_squared[1:],
        strict=True,
    ):
        if f_left == 0.0:
            roots.append(float(left))
        elif f_left * f_right < 0.0:
            roots.append(
                brentq(
                    lambda radius: float(
                        radial_velocity_squared(radius, potential, energy, angular_momentum)
                    ),
                    float(left),
                    float(right),
                )
            )
    if velocity_squared[-1] == 0.0:
        roots.append(float(grid[-1]))
    unique = np.unique(np.asarray(roots, dtype=float))
    if unique.size < 2:
        raise ValueError("A bound orbit requires two turning points in the scan interval.")
    return float(unique[0]), float(unique[-1])


def radial_period(
    potential: Callable[[Any], Any],
    energy: float,
    angular_momentum: float,
    pericenter: float,
    apocenter: float,
) -> float:
    """Return the complete radial period between supplied turning points.

    Parameters
    ----------
    potential : callable
        Spherical potential Phi(r) with the same energy zero and velocity-squared units as
        energy; accepts the supplied radius shape.
    energy : float or array_like
        Specific orbital energy, in the same velocity-squared units as the supplied
        potential.
    angular_momentum : float or array_like
        Nonnegative specific angular momentum, in length times velocity units consistent
        with the potential.
    pericenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    apocenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.

    Returns
    -------
    float
        Full radial oscillation period, twice the integral dr/|v_r|, in consistent
        length/velocity units.

    Notes
    -----
    These kernels assume a time-independent spherical potential and a bound orbit with two
    resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
    km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
    Gyr.
    """
    energy, angular_momentum = _validate_orbit(energy, angular_momentum)
    pericenter = float(pericenter)
    apocenter = float(apocenter)
    if (
        not np.isfinite(pericenter)
        or not np.isfinite(apocenter)
        or pericenter <= 0.0
        or pericenter >= apocenter
    ):
        raise ValueError("Turning points must satisfy 0 < pericenter < apocenter.")
    midpoint = 0.5 * (pericenter + apocenter)
    half_width = 0.5 * (apocenter - pericenter)

    def transformed_integrand(theta: float) -> float:
        """Remove the simple turning-point endpoint singularity from the radial time integral.

        Parameters
        ----------
        theta : float
            Dimensionless integration angle in [0, pi] mapping the two radial turning points to
            regular endpoints.

        Returns
        -------
        float
            One-way time element after r=midpoint-half_width*cos(theta).

        Notes
        -----
        These kernels assume a time-independent spherical potential and a bound orbit with two
        resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
        km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
        Gyr.
        """
        radius = midpoint - half_width * np.cos(theta)
        velocity_squared = float(
            radial_velocity_squared(radius, potential, energy, angular_momentum)
        )
        if velocity_squared <= 0.0:
            # Both numerator and radial speed vanish at a simple turning point.
            # Quad never relies on the exact endpoint value, so zero is a
            # stable limiting placeholder there.
            return 0.0
        return half_width * np.sin(theta) / np.sqrt(velocity_squared)

    one_way = quad(
        transformed_integrand,
        0.0,
        np.pi,
        epsabs=1.0e-10,
        epsrel=1.0e-10,
        limit=300,
    )[0]
    period = 2.0 * one_way
    if not np.isfinite(period) or period <= 0.0:
        raise ValueError("Orbit integration produced an invalid radial period.")
    return float(period)


def radial_shell_pdf(
    radius: Any,
    potential: Callable[[Any], Any],
    energy: float,
    angular_momentum: float,
    pericenter: float,
    apocenter: float,
    *,
    period: float | None = None,
) -> np.ndarray:
    """Return the time-averaged probability density per unit radius.

    Parameters
    ----------
    radius : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    potential : callable
        Spherical potential Phi(r) with the same energy zero and velocity-squared units as
        energy; accepts the supplied radius shape.
    energy : float or array_like
        Specific orbital energy, in the same velocity-squared units as the supplied
        potential.
    angular_momentum : float or array_like
        Nonnegative specific angular momentum, in length times velocity units consistent
        with the potential.
    pericenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    apocenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    period : float, optional
        Full radial period, in length/velocity units consistent with the potential. If
        omitted it is integrated from the supplied turning points.

    Returns
    -------
    numpy.ndarray
        Probability per unit radius, 2/(T_r*|v_r|) strictly inside the turning points and
        zero elsewhere.

    Notes
    -----
    These kernels assume a time-independent spherical potential and a bound orbit with two
    resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
    km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
    Gyr.
    """
    radius = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("Radii must be finite and positive.")
    radial_period_value = (
        radial_period(potential, energy, angular_momentum, pericenter, apocenter)
        if period is None
        else float(period)
    )
    if not np.isfinite(radial_period_value) or radial_period_value <= 0.0:
        raise ValueError("period must be finite and positive.")
    velocity_squared = radial_velocity_squared(radius, potential, energy, angular_momentum)
    inside = (radius > pericenter) & (radius < apocenter) & (velocity_squared > 0.0)
    result = np.zeros(radius.shape, dtype=float)
    result[inside] = 2.0 / (radial_period_value * np.sqrt(velocity_squared[inside]))
    return result


def orbit_radial_measure(
    radius_edges: Any,
    potential: Callable[[Any], Any],
    energy: float,
    angular_momentum: float,
    pericenter: float,
    apocenter: float,
) -> RadialMeasure:
    """Integrate normalized time fractions over radial shell edges.

    Parameters
    ----------
    radius_edges : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    potential : callable
        Spherical potential Phi(r) with the same energy zero and velocity-squared units as
        energy; accepts the supplied radius shape.
    energy : float or array_like
        Specific orbital energy, in the same velocity-squared units as the supplied
        potential.
    angular_momentum : float or array_like
        Nonnegative specific angular momentum, in length times velocity units consistent
        with the potential.
    pericenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    apocenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.

    Returns
    -------
    RadialMeasure
        Shell-center nodes and integrated, normalized orbital time fractions.

    Notes
    -----
    These kernels assume a time-independent spherical potential and a bound orbit with two
    resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
    km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
    Gyr.
    """
    edges = np.asarray(radius_edges, dtype=float)
    if (
        edges.ndim != 1
        or edges.size < 2
        or not np.all(np.isfinite(edges))
        or np.any(np.diff(edges) <= 0.0)
        or not np.isclose(edges[0], pericenter, rtol=1.0e-12, atol=0.0)
        or not np.isclose(edges[-1], apocenter, rtol=1.0e-12, atol=0.0)
    ):
        raise ValueError("Shell edges must increase from pericenter through apocenter.")
    period = radial_period(potential, energy, angular_momentum, pericenter, apocenter)
    midpoint = 0.5 * (pericenter + apocenter)
    half_width = 0.5 * (apocenter - pericenter)

    def theta(radius: float) -> float:
        """Map a shell boundary to the angle used for regularized orbit integration.

        Parameters
        ----------
        radius : float or array_like
            Positive radius or radial bound in one consistent length unit. Pericenter and
            apocenter bound the allowed radial motion.

        Returns
        -------
        float
            Angle in [0,pi] corresponding to the supplied radius.

        Notes
        -----
        These kernels assume a time-independent spherical potential and a bound orbit with two
        resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
        km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
        Gyr.
        """
        argument = np.clip((midpoint - radius) / half_width, -1.0, 1.0)
        return float(np.arccos(argument))

    def transformed_density(angle: float) -> float:
        """Integrate orbital time fraction with the transformed radial Jacobian.

        Parameters
        ----------
        angle : float
            Dimensionless integration angle in [0, pi] mapping the two radial turning points to
            regular endpoints.

        Returns
        -------
        float
            Orbital probability per unit angle, including the radial Jacobian.

        Notes
        -----
        These kernels assume a time-independent spherical potential and a bound orbit with two
        resolved turning points. Units remain caller-consistent: if radius is Mpc and velocity
        km/s, a returned period is Mpc/(km/s), requiring conversion before interpretation in
        Gyr.
        """
        radius = midpoint - half_width * np.cos(angle)
        velocity_squared = float(
            radial_velocity_squared(radius, potential, energy, angular_momentum)
        )
        if velocity_squared <= 0.0:
            return 0.0
        return 2.0 * half_width * np.sin(angle) / (period * np.sqrt(velocity_squared))

    weights = np.asarray(
        [
            quad(
                transformed_density,
                theta(float(left)),
                theta(float(right)),
                epsabs=1.0e-10,
                epsrel=1.0e-10,
                limit=200,
            )[0]
            for left, right in zip(edges[:-1], edges[1:], strict=True)
        ]
    )
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("Orbit shells have invalid total probability.")
    centers = 0.5 * (edges[:-1] + edges[1:])
    return RadialMeasure(radius=centers, weight=weights / total)


__all__ = [
    "orbit_radial_measure",
    "radial_period",
    "radial_shell_pdf",
    "radial_velocity_squared",
    "turning_points",
]
