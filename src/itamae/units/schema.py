"""Canonical internal unit schema used by all ITAMAE backends.

The schema is intentionally represented by dependency-free strings so that the
native backend does not import Astropy. Optional unit backends translate these
definitions into their own unit objects.
"""

from types import MappingProxyType

CANONICAL_UNIT_SCHEMA_VERSION = "1.0"

CANONICAL_UNITS = MappingProxyType(
    {
        "dimensionless": "",
        "mass": "Msun",
        "length": "Mpc",
        "velocity": "km / s",
        "time": "Gyr",
        "density": "Msun / Mpc3",
        "cross_section_per_mass": "cm2 / g",
    }
)


def canonical_unit(physical_type: str) -> str:
    """Return the canonical unit string for a physical quantity.

    Parameters
    ----------
    physical_type
        Key in :data:`CANONICAL_UNITS`.

    Returns
    -------
    str
        Backend-independent unit string. The empty string denotes a
        dimensionless value.

    Raises
    ------
    KeyError
        If ``physical_type`` is not part of the public schema.
    """
    try:
        return CANONICAL_UNITS[physical_type]
    except KeyError as exc:
        supported = ", ".join(sorted(CANONICAL_UNITS))
        raise KeyError(
            f"Unknown physical type {physical_type!r}; expected one of {supported}."
        ) from exc


__all__ = ["CANONICAL_UNITS", "CANONICAL_UNIT_SCHEMA_VERSION", "canonical_unit"]
