"""Astropy Quantity adapter for ITAMAE public interfaces."""

from dataclasses import dataclass

import astropy.units as u

from .native import _finite_real_array
from .schema import CANONICAL_UNITS, CANONICAL_UNIT_SCHEMA_VERSION, canonical_unit


def _astropy_unit(physical_type: str) -> u.UnitBase:
    """Return an Astropy unit for a canonical physical type.

    Parameters
    ----------
    physical_type : str
        Quantity name in the canonical unit schema, such as mass, length, time or density.

    Returns
    -------
    astropy.units.UnitBase
        Canonical physical unit, or the dimensionless unit for an empty schema value.

    Notes
    -----
    Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
    are assumed to already use that schema; only the quantity-aware backend performs
    dimensional conversion.
    """
    unit = canonical_unit(physical_type)
    return u.dimensionless_unscaled if unit == "" else u.Unit(unit)


@dataclass(frozen=True, slots=True)
class AstropyUnits:
    """Validate Quantity inputs and convert them to canonical units.

    Parameters
    ----------
    identifier : str
        Stable identifier of the physical and numerical prescription, used to distinguish
        calculations and cached data.

    Notes
    -----
    Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
    are assumed to already use that schema; only the quantity-aware backend performs
    dimensional conversion.
    """

    identifier: str = f"astropy-units:{CANONICAL_UNIT_SCHEMA_VERSION}"

    def to_internal(self, value, physical_type: str):
        """Convert a Quantity to a plain floating array in canonical units.

        Parameters
        ----------
        value : array_like
            Numerical input to be checked before use; shape and physical interpretation are set
            by the caller and retained by the validation helper.
        physical_type : str
            Quantity name in the canonical unit schema, such as mass, length, time or density.

        Returns
        -------
        numpy.ndarray
            Finite plain numerical array in the canonical unit for physical_type.

        Notes
        -----
        Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
        are assumed to already use that schema; only the quantity-aware backend performs
        dimensional conversion.
        """
        quantity = u.Quantity(value)
        return _finite_real_array(quantity.to_value(_astropy_unit(physical_type)))

    def from_internal(self, value, unit):
        """Convert a canonical floating value to an equivalent Astropy unit.

        Parameters
        ----------
        value : array_like
            Numerical input to be checked before use; shape and physical interpretation are set
            by the caller and retained by the validation helper.
        unit : object, optional
            Requested display unit. NativeUnits only accepts None or an empty unit; AstropyUnits
            requires a dimensionally equivalent unit.

        Returns
        -------
        numpy.ndarray or astropy.units.Quantity
            Values in the requested unit when supported; native output remains a plain canonical
            array.

        Notes
        -----
        Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
        are assumed to already use that schema; only the quantity-aware backend performs
        dimensional conversion.
        """
        target = u.Unit(unit)
        array = _finite_real_array(value)
        for physical_type in CANONICAL_UNITS:
            source = _astropy_unit(physical_type)
            if source.is_equivalent(target):
                return (array * source).to(target)
        raise u.UnitConversionError(f"No canonical ITAMAE unit is equivalent to {target}.")

    def validate(self, value, physical_type: str) -> None:
        """Raise when a value is dimensionally incompatible.

        Parameters
        ----------
        value : array_like
            Numerical input to be checked before use; shape and physical interpretation are set
            by the caller and retained by the validation helper.
        physical_type : str
            Quantity name in the canonical unit schema, such as mass, length, time or density.

        Notes
        -----
        Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
        are assumed to already use that schema; only the quantity-aware backend performs
        dimensional conversion.
        """
        self.to_internal(value, physical_type)
