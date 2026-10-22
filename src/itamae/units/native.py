"""Lightweight unit backend using canonical floating-point units."""

from dataclasses import dataclass

import numpy as np

from .schema import CANONICAL_UNIT_SCHEMA_VERSION, canonical_unit


def _finite_real_array(value):
    """Reject non-real or non-finite quantities before a unit-boundary conversion.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.

    Returns
    -------
    numpy.ndarray
        Finite floating values retaining the input shape.

    Notes
    -----
    Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
    are assumed to already use that schema; only the quantity-aware backend performs
    dimensional conversion.
    """
    array = np.asarray(value)
    if array.dtype.kind not in "biuf":
        raise ValueError("Unit input/output must contain real numeric values.")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError("Unit input/output contains non-finite values.")
    return array


@dataclass(frozen=True, slots=True)
class NativeUnits:
    """Interpret numeric inputs in the canonical ITAMAE unit system.

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

        Notes
        -----
        Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
        are assumed to already use that schema; only the quantity-aware backend performs
        dimensional conversion.
        """
        canonical_unit(physical_type)
        return _finite_real_array(value)

    def from_internal(self, value, unit=None):
        """Return a finite internal value without attaching units.

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
        The native backend cannot convert to a non-canonical display unit.
        Callers requesting explicit output units should use ``AstropyUnits``.

        Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
        are assumed to already use that schema; only the quantity-aware backend performs
        dimensional conversion.
        """
        if unit not in (None, ""):
            raise ValueError("NativeUnits cannot attach or convert explicit output units.")
        return _finite_real_array(value)

    def validate(self, value, physical_type: str) -> None:
        """Validate that a value can be represented numerically.

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
