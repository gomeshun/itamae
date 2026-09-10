"""Numerically stable Navarro-Frenk-White profile utilities."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

_G_MPC_KMS2_MSUN = 4.30091e-9


def nfw_mass_function(x):
    """Return ``ln(1+x)-x/(1+x)`` for finite nonnegative ``x``.

    A small-radius series avoids subtracting two nearly equal terms. Values
    smaller than the floating-point subnormal range can still underflow.
    """
    x = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x)) or np.any(x < 0.0):
        raise ValueError("NFW radius ratio must be finite and nonnegative.")
    result = np.empty_like(x)
    small = x < 0.01
    # f(x) = sum_{n=2}^infinity (-1)^n (n-1)/n x^n.
    coefficient = [(-1.0) ** n * (n - 1.0) / n for n in range(2, 11)]
    result[small] = x[small] ** 2 * np.polynomial.polynomial.polyval(x[small], coefficient)
    result[~small] = np.log1p(x[~small]) - x[~small] / (1.0 + x[~small])
    return result


def invert_nfw_mass_function(y):
    """Invert the NFW enclosed-mass function over representable positive radii.

    Solve in log-radius so an absolute tolerance in radius cannot erase a
    small, positive solution. The output retains the input shape.
    """
    y = np.asarray(y, dtype=float)
    if not np.all(np.isfinite(y)) or np.any(y < 0.0):
        raise ValueError("Enclosed-mass function values must be finite and nonnegative.")
    max_radius = np.finfo(float).max
    max_log_radius = np.log(max_radius)
    max_enclosed = float(nfw_mass_function(max_radius))
    if np.any(y > max_enclosed):
        raise ValueError(
            "The inverse NFW radius is outside the representable floating-point range."
        )

    def one(value: float) -> float:
        if value == 0.0:
            return 0.0
        if value == max_enclosed:
            return max_radius
        # f(x) <= x^2/2 and f(exp(y+1)) >= y bound the positive root.
        lower = 0.5 * (np.log(2.0) + np.log(value)) - np.log(2.0)
        upper = min(value + 1.0, max_log_radius)

        def residual(log_radius):
            radius = max_radius if log_radius == max_log_radius else np.exp(log_radius)
            return float(nfw_mass_function(radius) / value - 1.0)

        root = brentq(residual, lower, upper, xtol=5.0e-14, rtol=4.0 * np.finfo(float).eps)
        return float(np.exp(root))

    out = np.vectorize(one, otypes=[float])(y)
    return float(out) if out.ndim == 0 else out


@dataclass(frozen=True, slots=True)
class NFWProfile:
    """Spherical NFW profile parameterized by scale radius and density."""

    r_s: float
    rho_s: float

    def __post_init__(self):
        """Require a physical finite positive scale radius and density."""
        for name in ("r_s", "rho_s"):
            value = getattr(self, name)
            if np.ndim(value) != 0 or not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")

    def enclosed_mass(self, r):
        """Return enclosed mass in the profile's mass unit."""
        r = np.asarray(r, dtype=float)
        if not np.all(np.isfinite(r)) or np.any(r < 0.0):
            raise ValueError("Radius must be finite and nonnegative.")
        return 4.0 * np.pi * self.rho_s * self.r_s**3 * nfw_mass_function(r / self.r_s)

    def density(self, r):
        """Return density at positive radius."""
        r = np.asarray(r, dtype=float)
        if not np.all(np.isfinite(r)) or np.any(r <= 0.0):
            raise ValueError("Density requires finite positive radii; it is singular at zero.")
        x = r / self.r_s
        return self.rho_s / (x * (1.0 + x) ** 2)

    def potential(self, r):
        """Return gravitational potential with zero at infinity."""
        r = np.asarray(r, dtype=float)
        if not np.all(np.isfinite(r)) or np.any(r < 0.0):
            raise ValueError("Radius must be finite and nonnegative.")
        x = r / self.r_s
        ratio = np.divide(np.log1p(x), x, out=np.ones_like(x), where=x != 0.0)
        return -4.0 * np.pi * _G_MPC_KMS2_MSUN * self.rho_s * self.r_s**2 * ratio
