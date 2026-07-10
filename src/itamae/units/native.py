"""Lightweight unit backend using canonical floating-point units."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class NativeUnits:
    """Interpret numeric inputs in the canonical ITAMAE unit system."""

    identifier: str = "native-v1"

    def to_internal(self, value, physical_type: str):
        """Convert a numeric value to a floating NumPy array.

        Parameters
        ----------
        value
            Scalar or array expressed in the documented canonical unit.
        physical_type
            Name of the physical quantity. It is retained for API symmetry.

        Returns
        -------
        numpy.ndarray
            Floating representation of ``value``.
        """
        del physical_type
        return np.asarray(value, dtype=float)

    def from_internal(self, value, unit=None):
        """Return an internal value without attaching units."""
        del unit
        return np.asarray(value, dtype=float)

    def validate(self, value, physical_type: str) -> None:
        """Validate that a value can be represented numerically."""
        del physical_type
        arr = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(arr)):
            raise ValueError("Unit input contains non-finite values.")
