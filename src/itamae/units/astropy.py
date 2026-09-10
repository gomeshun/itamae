"""Astropy Quantity adapter for ITAMAE public interfaces."""

from dataclasses import dataclass

import astropy.units as u

from .native import _finite_real_array
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
        return _finite_real_array(quantity.to_value(_astropy_unit(physical_type)))

    def from_internal(self, value, unit):
        """Convert a canonical floating value to an equivalent Astropy unit."""
        target = u.Unit(unit)
        array = _finite_real_array(value)
        for physical_type in CANONICAL_UNITS:
            source = _astropy_unit(physical_type)
            if source.is_equivalent(target):
                return (array * source).to(target)
        raise u.UnitConversionError(f"No canonical ITAMAE unit is equivalent to {target}.")

    def validate(self, value, physical_type: str) -> None:
        """Raise when a value is dimensionally incompatible."""
        self.to_internal(value, physical_type)
