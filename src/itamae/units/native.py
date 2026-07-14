"""Lightweight unit backend using canonical floating-point units."""

from dataclasses import dataclass

import numpy as np

from .schema import CANONICAL_UNIT_SCHEMA_VERSION, canonical_unit


@dataclass(frozen=True, slots=True)
class NativeUnits:
    """Interpret numeric inputs in the canonical ITAMAE unit system."""

    identifier: str = f"native-units:{CANONICAL_UNIT_SCHEMA_VERSION}"

    def to_internal(self, value, physical_type: str):
        """Convert a numeric value to a floating NumPy array.

        Parameters
        ----------
        value
            Scalar or array expressed in the documented canonical unit.
        physical_type
            Name of the physical quantity in the canonical unit schema.

        Returns
        -------
        numpy.ndarray
            Floating representation of ``value``.
        """
        canonical_unit(physical_type)
        array = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError("Unit input contains non-finite values.")
        return array

    def from_internal(self, value, unit=None):
        """Return a finite internal value without attaching units.

        Notes
        -----
        The native backend cannot convert to a non-canonical display unit.
        Callers requesting explicit output units should use ``AstropyUnits``.
        """
        if unit not in (None, ""):
            raise ValueError("NativeUnits cannot attach or convert explicit output units.")
        array = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError("Unit output contains non-finite values.")
        return array

    def validate(self, value, physical_type: str) -> None:
        """Validate that a value can be represented numerically."""
        self.to_internal(value, physical_type)
