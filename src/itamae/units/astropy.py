"""Astropy Quantity adapter for ITAMAE public interfaces."""

from dataclasses import dataclass

import astropy.units as u
import numpy as np

_UNITS = {
    "mass": u.Msun,
    "length": u.Mpc,
    "velocity": u.km / u.s,
    "time": u.Gyr,
    "density": u.Msun / u.Mpc**3,
    "cross_section_per_mass": u.cm**2 / u.g,
}


@dataclass(frozen=True, slots=True)
class AstropyUnits:
    """Validate Quantity inputs and convert them to canonical units."""

    identifier: str = "astropy-v1"

    def to_internal(self, value, physical_type: str):
        """Convert a Quantity to a plain floating array in canonical units."""
        if physical_type not in _UNITS:
            raise KeyError(f"Unknown physical type: {physical_type}")
        quantity = u.Quantity(value)
        return np.asarray(quantity.to_value(_UNITS[physical_type]), dtype=float)

    def from_internal(self, value, unit):
        """Attach an Astropy unit to a canonical floating value."""
        return np.asarray(value, dtype=float) * u.Unit(unit)

    def validate(self, value, physical_type: str) -> None:
        """Raise when a value is dimensionally incompatible."""
        self.to_internal(value, physical_type)
