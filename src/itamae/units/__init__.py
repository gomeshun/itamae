"""Unit conversion backends."""

from typing import TYPE_CHECKING, Any

from .native import NativeUnits
from .schema import CANONICAL_UNITS, CANONICAL_UNIT_SCHEMA_VERSION, canonical_unit

if TYPE_CHECKING:
    from .astropy import AstropyUnits


def __getattr__(name: str) -> Any:
    """Load a requested optional unit converter by its public name.

    Parameters
    ----------
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.

    Returns
    -------
    type
        Requested unit-backend class; unknown names raise AttributeError.

    Notes
    -----
    Canonical values use Msun, Mpc, km/s, Gyr, Msun/Mpc^3 and cm^2/g. Native numeric values
    are assumed to already use that schema; only the quantity-aware backend performs
    dimensional conversion.
    """
    if name == "AstropyUnits":
        from .astropy import AstropyUnits

        return AstropyUnits
    raise AttributeError(name)


__all__ = [
    "AstropyUnits",
    "CANONICAL_UNITS",
    "CANONICAL_UNIT_SCHEMA_VERSION",
    "NativeUnits",
    "canonical_unit",
]
