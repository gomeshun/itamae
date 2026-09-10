"""Adapter from existing callables to the ITAMAE variance protocol."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np

from ._validation import aligned, coordinates, numerical_errors


@dataclass(frozen=True, slots=True)
class CallableVarianceModel:
    r"""Wrap model-supplied variance functions behind a stable interface.

    Parameters
    ----------
    identifier
        Stable description of the physical and numerical prescription.
    sigma_function
        Callable with signature ``sigma_function(mass, z)`` returning
        :math:`\sigma(M,z)`.
    derivative_function
        Callable with signature ``derivative_function(mass, z)`` returning
        :math:`\mathrm{d}S/\mathrm{d}M`, where :math:`S=\sigma^2`.

    Notes
    -----
    This class is primarily a migration tool. It lets existing SASHIMI formulae
    participate in ITAMAE pipelines while retaining explicit ownership of their
    transfer, filter, concentration, and normalization assumptions.
    """

    identifier: str
    sigma_function: Callable[[Any, Any], Any]
    derivative_function: Callable[[Any, Any], Any]

    def __post_init__(self):
        """Require a meaningful identifier and callable numerical implementations."""
        if not isinstance(self.identifier, str) or not self.identifier.strip():
            raise ValueError("identifier must be a non-empty string.")
        if not callable(self.sigma_function) or not callable(self.derivative_function):
            raise TypeError("sigma_function and derivative_function must be callable.")

    def sigma(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        """Return the rms fluctuation supplied by the wrapped implementation."""
        mass_array, _ = coordinates(mass, z)
        with numerical_errors("Sigma callback"):
            return aligned(
                self.sigma_function(mass, z),
                mass_array.shape,
                "Sigma callback",
                nonnegative=True,
                scalar_ok=True,
            )

    def variance(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        """Return the square of the wrapped rms fluctuation."""
        sigma = self.sigma(mass, z)
        with numerical_errors("Variance"):
            return sigma * sigma

    def dvariance_dmass(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        """Return the wrapped derivative of variance with respect to mass."""
        mass_array, _ = coordinates(mass, z)
        with numerical_errors("Variance derivative callback"):
            return aligned(
                self.derivative_function(mass, z),
                mass_array.shape,
                "Variance derivative callback",
                scalar_ok=True,
            )
