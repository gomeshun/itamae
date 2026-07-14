"""Astropy Quantity adapter for ITAMAE public interfaces."""

from dataclasses import dataclass

import astropy.units as u
import numpy as np

from .schema import CANONICAL_UNITS, CANONICAL_UNIT_SCHEMA_VERSION, canonical_unit


def _astropy_unit(physical_type: str) -> u.UnitBase:
    """Return an Astropy unit for a canonical physical type."""
    unit = canonical_unit(physical_type)
    return u.dimensionless_unscaled if unit == "" else u.Unit(unit)


@dataclass(frozen=True, slots=True)
class AstropyUnits:
    """Validate Quantity inputs and convert them to canonical units."""

    identifier: str = f"astropy-units:{CANONICAL_UNIT_SCHEMA_VERSION}"

    def to_internal(self, value, physical_type: str):
        """Convert a Quantity to a plain floating array in canonical units."""
        quantity = u.Quantity(value)
        array = np.asarray(quantity.to_value(_astropy_unit(physical_type)), dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError("Unit input contains non-finite values.")
        return array

    def from_internal(self, value, unit):
        """Convert a canonical floating value to an equivalent Astropy unit."""
        target = u.Unit(unit)
        array = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError("Unit output contains non-finite values.")
        for physical_type in CANONICAL_UNITS:
            source = _astropy_unit(physical_type)
            if source.is_equivalent(target):
                return (array * source).to(target)
        raise u.UnitConversionError(f"No canonical ITAMAE unit is equivalent to {target}.")

    def validate(self, value, physical_type: str) -> None:
        """Raise when a value is dimensionally incompatible."""
        self.to_internal(value, physical_type)
