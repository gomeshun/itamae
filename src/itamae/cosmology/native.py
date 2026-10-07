"""Native flat-LCDM cosmology implementation."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

_G_MPC_KMS2_MSUN = 4.30091e-9
_KMS_MPC_TO_GYR_INV = 1.0227121650537077e-3


def _redshift(value):
    array = np.asarray(value)
    if array.dtype.kind not in "biuf":
        raise ValueError("Redshift must be real, finite and greater than -1.")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array <= -1.0):
        raise ValueError("Redshift must be real, finite and greater than -1.")
    return array


def _finite_background(evaluate, quantity):
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            result = evaluate()
    except FloatingPointError as error:
        raise ValueError(f"{quantity} exceeds the finite numerical background domain.") from error
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{quantity} must be finite.")
    return result


@dataclass(frozen=True, slots=True)
class NativeFlatLCDM:
    """Flat-LCDM background with NumPy/SciPy only.

    Parameters
    ----------
    omega_m0
        Present-day matter density fraction.
    h
        Reduced Hubble constant.
    """

    omega_m0: float = 0.315
    h: float = 0.674

    def __post_init__(self):
        """Validate the matter-plus-nonnegative-Lambda background domain."""
        for name in ("omega_m0", "h"):
            value = np.asarray(getattr(self, name))
            if value.ndim != 0 or value.dtype.kind not in "iuf" or not np.isfinite(value):
                raise ValueError(f"{name} must be a finite real scalar.")
            object.__setattr__(self, name, float(value))
        if not 0.0 < self.omega_m0 <= 1.0:
            raise ValueError("omega_m0 must satisfy 0 < omega_m0 <= 1.")
        if self.h <= 0.0:
            raise ValueError("h must be positive.")

    @property
    def identifier(self) -> str:
        """Return a stable identifier including cosmological parameters."""
        return f"native-flatlcdm:Om={self.omega_m0}:h={self.h}"

    @property
    def omega_lambda0(self) -> float:
        """Return the present-day dark-energy density fraction."""
        return 1.0 - self.omega_m0

    def e2(self, z):
        """Return the dimensionless squared expansion rate ``E(z)^2``."""
        z = _redshift(z)
        try:
            with np.errstate(over="raise", invalid="raise"):
                result = self.omega_m0 * (1.0 + z) ** 3 + self.omega_lambda0
        except FloatingPointError as error:
            raise ValueError(
                "Redshift exceeds the finite numerical expansion-rate domain."
            ) from error
        return result

    def H(self, z):
        """Return the Hubble rate in km s^-1 Mpc^-1."""
        return _finite_background(lambda: 100.0 * self.h * np.sqrt(self.e2(z)), "Hubble rate")

    def rho_crit(self, z):
        """Return critical density in Msun Mpc^-3."""
        return _finite_background(
            lambda: 3.0 * self.H(z) ** 2 / (8.0 * np.pi * _G_MPC_KMS2_MSUN),
            "Critical density",
        )

    def rho_m(self, z):
        """Return physical matter density in Msun Mpc^-3."""
        return _finite_background(
            lambda: self.omega_m0 * self.rho_crit(0.0) * (1.0 + _redshift(z)) ** 3,
            "Matter density",
        )

    def omega_m(self, z):
        """Return the redshift-dependent matter density fraction."""
        z = _redshift(z)
        return self.omega_m0 * (1.0 + z) ** 3 / self.e2(z)

    def growth_factor(self, z):
        """Return the Carroll-Press-Turner growth approximation normalized at z=0."""
        z = _redshift(z)
        om = self.omega_m(z)
        ol = 1.0 - om
        g = 2.5 * om / (om ** (4.0 / 7.0) - ol + (1.0 + om / 2.0) * (1.0 + ol / 70.0))
        om0 = self.omega_m0
        ol0 = self.omega_lambda0
        g0 = 2.5 * om0 / (om0 ** (4.0 / 7.0) - ol0 + (1.0 + om0 / 2.0) * (1.0 + ol0 / 70.0))
        return g / (g0 * (1.0 + z))

    def collapse_threshold(self, z):
        """Return the spherical-collapse threshold scaled by the growth factor."""
        return 1.686 / self.growth_factor(z)

    def cosmic_time(self, z):
        """Return cosmic age in Gyr.

        Notes
        -----
        The integral is evaluated independently for scalar inputs. Array inputs
        are vectorized to preserve a simple public API.
        """

        def one(zi: float) -> float:
            def integrand(zp: float) -> float:
                return 1.0 / ((1.0 + zp) * self.H(zp) * _KMS_MPC_TO_GYR_INV)

            return quad(integrand, zi, np.inf, epsabs=1e-9, epsrel=1e-9)[0]

        arr = _redshift(z)
        out = np.vectorize(one, otypes=[float])(arr)
        return float(out) if out.ndim == 0 else out

    def lookback_time(self, z):
        """Return lookback time in Gyr."""
        return self.cosmic_time(0.0) - self.cosmic_time(z)
