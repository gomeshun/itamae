"""Numerically stable Navarro-Frenk-White profile utilities."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

_G_MPC_KMS2_MSUN = 4.30091e-9


def nfw_mass_function(x):
    """Return ``ln(1+x)-x/(1+x)`` for nonnegative ``x``."""
    x = np.asarray(x, dtype=float)
    if np.any(x < 0.0):
        raise ValueError("NFW radius ratio must be nonnegative.")
    return np.log1p(x) - x / (1.0 + x)


def invert_nfw_mass_function(y):
    """Invert the monotonic NFW enclosed-mass function."""
    y = np.asarray(y, dtype=float)
    if np.any(y < 0.0):
        raise ValueError("Enclosed-mass function values must be nonnegative.")

    def one(value: float) -> float:
        if value == 0.0:
            return 0.0
        upper = max(1.0, np.exp(min(value + 1.0, 700.0)))
        while nfw_mass_function(upper) < value:
            upper *= 2.0
        return brentq(lambda x: float(nfw_mass_function(x) - value), 0.0, upper)

    out = np.vectorize(one, otypes=[float])(y)
    return float(out) if out.ndim == 0 else out


@dataclass(frozen=True, slots=True)
class NFWProfile:
    """Spherical NFW profile parameterized by scale radius and density."""

    r_s: float
    rho_s: float

    def enclosed_mass(self, r):
        """Return enclosed mass in the profile's mass unit."""
        r = np.asarray(r, dtype=float)
        if np.any(r < 0.0):
            raise ValueError("Radius must be nonnegative.")
        return 4.0 * np.pi * self.rho_s * self.r_s**3 * nfw_mass_function(r / self.r_s)

    def density(self, r):
        """Return density at positive radius."""
        r = np.asarray(r, dtype=float)
        if np.any(r <= 0.0):
            raise ValueError("Density is singular at nonpositive radius.")
        x = r / self.r_s
        return self.rho_s / (x * (1.0 + x) ** 2)

    def potential(self, r):
        """Return gravitational potential with zero at infinity."""
        r = np.asarray(r, dtype=float)
        x = r / self.r_s
        ratio = np.where(x == 0.0, 1.0, np.log1p(x) / x)
        return -4.0 * np.pi * _G_MPC_KMS2_MSUN * self.rho_s * self.r_s**2 * ratio
