"""Power-spectrum integration behind the common variance protocol."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import simpson

from itamae.power.windows import SharpKWindow
from itamae.protocols import PowerSpectrum, WindowFunction


@dataclass(frozen=True, slots=True)
class IntegratedVarianceModel:
    r"""Integrate a supplied power spectrum and smoothing window.

    Parameters
    ----------
    power
        Model-supplied power spectrum. ITAMAE does not select its cosmology,
        normalization, or dark-matter transfer function.
    window
        Dimensionless Fourier-space smoothing window.
    rho_mean
        Mean density used in
        :math:`R=(3M/(4\pi\rho_\mathrm{mean}))^{1/3}`. Mass and density must use
        mutually consistent units.
    k_min, k_max, n_k
        Positive integration bounds and logarithmic grid size.
    filter_scale
        Dimensionless mapping in ``x = k R / filter_scale``. For a top-hat it
        is normally one. A sharp-k SASHIMI variant may explicitly supply its
        calibrated mass-assignment coefficient.
    growth_function
        Optional callable returning a growth factor normalized consistently
        with the supplied redshift-zero spectrum.
    growth_identifier
        Required stable identifier when ``growth_function`` is supplied.
    derivative_step
        Symmetric logarithmic step used for ``dvariance_dmass`` with smooth
        windows. Sharp-k windows use their exact moving-boundary derivative.
    chunk_size
        Maximum number of masses integrated in one vectorized allocation.

    Notes
    -----
    The implemented convention is

    .. math::

       S(M)=\int d\ln k\,\frac{k^3P(k)}{2\pi^2}W^2(kR/c).

    WDM/FDM transfer functions and their physical parameters remain owned by
    the calling SASHIMI package.
    """

    power: PowerSpectrum
    window: WindowFunction
    rho_mean: float
    k_min: float
    k_max: float
    n_k: int = 2049
    filter_scale: float = 1.0
    growth_function: Callable[[Any], Any] | None = None
    growth_identifier: str | None = None
    derivative_step: float = 1.0e-4
    chunk_size: int = 256

    def __post_init__(self) -> None:
        """Validate physical and numerical configuration."""
        if not isinstance(self.power, PowerSpectrum):
            raise TypeError("power must implement the ITAMAE power-spectrum protocol.")
        if not isinstance(self.window, WindowFunction):
            raise TypeError("window must implement the ITAMAE window protocol.")
        for name in ("rho_mean", "k_min", "k_max", "filter_scale", "derivative_step"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.k_min >= self.k_max:
            raise ValueError("k_min must be smaller than k_max.")
        for name in ("n_k", "chunk_size"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"{name} must be an integer.")
            if int(value) < (3 if name == "n_k" else 1):
                raise ValueError(f"{name} is too small.")
        if self.growth_function is None and self.growth_identifier is not None:
            raise ValueError("growth_identifier requires growth_function.")
        if self.growth_function is not None and (
            not callable(self.growth_function)
            or not isinstance(self.growth_identifier, str)
            or not self.growth_identifier.strip()
        ):
            raise ValueError("A callable growth_function requires a non-empty growth_identifier.")

    @property
    def identifier(self) -> str:
        """Return a complete numerical and physical-component identifier."""
        growth = self.growth_identifier or "growth:redshift-independent"
        return (
            "integrated-variance:v2;"
            f"power=({self.power.identifier});window=({self.window.identifier});"
            f"rho_mean={self.rho_mean:.17g};k=[{self.k_min:.17g},{self.k_max:.17g}];"
            f"n_k={self.n_k};filter_scale={self.filter_scale:.17g};{growth}"
        )

    def _validated_mass(self, mass: Any) -> np.ndarray:
        """Return finite positive mass values as a floating-point array."""
        values = np.asarray(mass, dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("Masses must be finite and positive.")
        return values

    def _sharp_k_cutoff(self, mass: np.ndarray) -> np.ndarray:
        """Return the mass-dependent sharp-k integration boundary."""
        radius = (3.0 * mass / (4.0 * np.pi * self.rho_mean)) ** (1.0 / 3.0)
        return self.filter_scale / radius

    def _sharp_k_variance_z0(self, mass: Any) -> np.ndarray:
        """Integrate a sharp-k spectrum with the cutoff as an exact endpoint.

        A discontinuous window sampled on one fixed wavenumber grid makes the
        numerical variance piecewise constant as mass changes. Integrating each
        mass to its physical cutoff instead keeps the result continuous and
        makes the moving-boundary derivative well defined.
        """
        values = self._validated_mass(mass)
        shape = values.shape
        flattened = values.reshape(-1)
        cutoff = np.minimum(self._sharp_k_cutoff(flattened), self.k_max)
        result = np.zeros(flattened.shape, dtype=float)
        active = cutoff > self.k_min
        fractions = np.linspace(0.0, 1.0, self.n_k)
        log_k_min = np.log(self.k_min)

        active_indices = np.flatnonzero(active)
        for start in range(0, active_indices.size, self.chunk_size):
            indices = active_indices[start : start + self.chunk_size]
            upper = cutoff[indices]
            log_upper = np.log(upper)
            log_k = log_k_min + (log_upper[:, None] - log_k_min) * fractions
            k = np.exp(log_k)
            # Assigning the exact endpoints avoids one-ulp excursions beyond a
            # strict tabulated-spectrum domain.
            k[:, 0] = self.k_min
            k[:, -1] = upper
            dimensionless_power = k**3 * np.asarray(self.power(k), dtype=float) / (2.0 * np.pi**2)
            if (
                dimensionless_power.shape != k.shape
                or not np.all(np.isfinite(dimensionless_power))
                or np.any(dimensionless_power < 0.0)
            ):
                raise ValueError("Power spectrum must return aligned finite nonnegative values.")
            result[indices] = simpson(dimensionless_power, x=log_k, axis=1)

        if not np.all(np.isfinite(result)) or np.any(result < 0.0):
            raise ValueError("Variance integration produced invalid values.")
        return result.reshape(shape)

    def _variance_z0(self, mass: Any) -> np.ndarray:
        """Integrate the redshift-zero variance in bounded memory chunks."""
        if isinstance(self.window, SharpKWindow):
            return self._sharp_k_variance_z0(mass)

        values = self._validated_mass(mass)
        shape = values.shape
        flattened = values.reshape(-1)
        result = np.empty(flattened.shape, dtype=float)
        # geomspace assigns the requested endpoints exactly. Re-exponentiating
        # a logarithmic linspace can cross a strict tabulated-spectrum boundary
        # by one floating-point ulp.
        k = np.geomspace(self.k_min, self.k_max, self.n_k)
        log_k = np.log(k)
        dimensionless_power = k**3 * np.asarray(self.power(k), dtype=float) / (2.0 * np.pi**2)
        if (
            dimensionless_power.shape != k.shape
            or not np.all(np.isfinite(dimensionless_power))
            or np.any(dimensionless_power < 0.0)
        ):
            raise ValueError("Power spectrum must return aligned finite nonnegative values.")

        for start in range(0, flattened.size, self.chunk_size):
            stop = min(start + self.chunk_size, flattened.size)
            radius = (3.0 * flattened[start:stop] / (4.0 * np.pi * self.rho_mean)) ** (1.0 / 3.0)
            argument = radius[:, None] * k[None, :] / self.filter_scale
            window = np.asarray(self.window(argument), dtype=float)
            integrand = dimensionless_power[None, :] * window * window
            result[start:stop] = simpson(integrand, x=log_k, axis=1)
        if not np.all(np.isfinite(result)) or np.any(result < 0.0):
            raise ValueError("Variance integration produced invalid values.")
        return result.reshape(shape)

    def _growth(self, redshift: Any) -> np.ndarray:
        """Evaluate and validate the optional growth factor."""
        if self.growth_function is None:
            redshift_array = np.asarray(redshift, dtype=float)
            if not np.all(np.isfinite(redshift_array)):
                raise ValueError("Redshift must be finite.")
            return np.ones(redshift_array.shape, dtype=float)
        growth = np.asarray(self.growth_function(redshift), dtype=float)
        if not np.all(np.isfinite(growth)) or np.any(growth < 0.0):
            raise ValueError("Growth function must return finite nonnegative values.")
        return growth

    def variance(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        r"""Return :math:`S(M,z)` with mass/redshift broadcasting."""
        mass_array, redshift = np.broadcast_arrays(
            np.asarray(mass, dtype=float), np.asarray(z, dtype=float)
        )
        return self._variance_z0(mass_array) * self._growth(redshift) ** 2

    def sigma(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        r"""Return :math:`\sigma(M,z)=\sqrt{S(M,z)}`."""
        return np.sqrt(self.variance(mass, z))

    def dvariance_dmass(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        r"""Return :math:`dS/dM` under the configured window convention.

        Smooth windows use a symmetric logarithmic difference. For a sharp-k
        window, ITAMAE evaluates the exact moving-boundary term

        .. math::

           \frac{dS}{dM}
           = -\frac{k_c^3 P(k_c)}{6\pi^2 M}D(z)^2,

        when the physical cutoff lies inside the configured power-spectrum
        interval. Outside that interval the truncated integral is constant.
        """
        mass_array, redshift = np.broadcast_arrays(
            np.asarray(mass, dtype=float), np.asarray(z, dtype=float)
        )
        self._validated_mass(mass_array)
        if isinstance(self.window, SharpKWindow):
            cutoff = self._sharp_k_cutoff(mass_array)
            growth = np.broadcast_to(self._growth(redshift), mass_array.shape)
            result = np.zeros(mass_array.shape, dtype=float)
            active = (cutoff > self.k_min) & (cutoff < self.k_max)
            if np.any(active):
                active_cutoff = cutoff[active]
                boundary_power = np.asarray(self.power(active_cutoff), dtype=float)
                if (
                    boundary_power.shape != active_cutoff.shape
                    or not np.all(np.isfinite(boundary_power))
                    or np.any(boundary_power < 0.0)
                ):
                    raise ValueError(
                        "Power spectrum must return aligned finite nonnegative values."
                    )
                result[active] = (
                    -(active_cutoff**3)
                    * boundary_power
                    * growth[active] ** 2
                    / (6.0 * np.pi**2 * mass_array[active])
                )
            return result

        upper_mass = mass_array * np.exp(self.derivative_step)
        lower_mass = mass_array * np.exp(-self.derivative_step)
        upper = self.variance(upper_mass, redshift)
        lower = self.variance(lower_mass, redshift)
        return (upper - lower) / (upper_mass - lower_mass)


__all__ = ["IntegratedVarianceModel"]
