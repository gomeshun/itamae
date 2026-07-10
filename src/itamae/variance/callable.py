"""Adapter from existing callables to the ITAMAE variance protocol."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True, slots=True)
class CallableVarianceModel:
    """Wrap model-supplied variance functions behind a stable interface.

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

    def sigma(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        """Return the rms fluctuation supplied by the wrapped implementation."""
        return np.asarray(self.sigma_function(mass, z), dtype=float)

    def variance(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        """Return the square of the wrapped rms fluctuation."""
        sigma = self.sigma(mass, z)
        return sigma * sigma

    def dvariance_dmass(self, mass: Any, z: Any = 0.0) -> np.ndarray:
        """Return the wrapped derivative of variance with respect to mass."""
        return np.asarray(self.derivative_function(mass, z), dtype=float)
