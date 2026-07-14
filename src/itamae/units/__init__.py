"""Unit conversion backends."""

from typing import TYPE_CHECKING, Any

from .native import NativeUnits
from .schema import CANONICAL_UNITS, CANONICAL_UNIT_SCHEMA_VERSION, canonical_unit

if TYPE_CHECKING:
    from .astropy import AstropyUnits


def __getattr__(name: str) -> Any:
    """Load optional unit backends lazily."""
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
