"""Structural interfaces for cosmology and unit backends."""

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class CosmologyBackend(Protocol):
    """Provide background-cosmology quantities in canonical ITAMAE units."""

    @property
    def identifier(self) -> str:
        """Return a stable identifier suitable for metadata and cache keys."""
        ...

    def H(self, z: Any) -> Any:
        """Return the Hubble rate in km s^-1 Mpc^-1."""
        ...

    def rho_crit(self, z: Any) -> Any:
        """Return critical density in Msun Mpc^-3."""
        ...

    def rho_m(self, z: Any) -> Any:
        """Return physical matter density in Msun Mpc^-3."""
        ...

    def omega_m(self, z: Any) -> Any:
        """Return the matter density fraction."""
        ...

    def growth_factor(self, z: Any) -> Any:
        """Return the linear growth factor normalized at redshift zero."""
        ...

    def collapse_threshold(self, z: Any) -> Any:
        """Return the spherical-collapse threshold at redshift ``z``."""
        ...

    def cosmic_time(self, z: Any) -> Any:
        """Return cosmic age in Gyr."""
        ...

    def lookback_time(self, z: Any) -> Any:
        """Return lookback time in Gyr."""
        ...


@runtime_checkable
class UnitBackend(Protocol):
    """Convert public values to and from the canonical internal unit system."""

    @property
    def identifier(self) -> str:
        """Return a stable identifier including the unit-schema version."""
        ...

    def to_internal(self, value: Any, physical_type: str) -> Any:
        """Convert a value to a plain array in canonical units."""
        ...

    def from_internal(self, value: Any, unit: Any = None) -> Any:
        """Convert a canonical value to a requested public unit."""
        ...

    def validate(self, value: Any, physical_type: str) -> None:
        """Validate dimensional compatibility and finite values."""
        ...


__all__ = ["CosmologyBackend", "UnitBackend"]
