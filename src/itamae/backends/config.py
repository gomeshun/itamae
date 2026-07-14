"""Immutable backend selection for ITAMAE calculations."""

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Mapping

from itamae.protocols import CosmologyBackend, UnitBackend


@dataclass(frozen=True, slots=True)
class BackendConfig:
    """Collect backend implementations used by one calculation.

    Parameters
    ----------
    cosmology
        Object implementing the cosmology backend protocol.
    units
        Object implementing the unit backend protocol.
    array
        Name of the numerical array backend. NumPy is the initial backend.
    """

    cosmology: CosmologyBackend
    units: UnitBackend
    array: str = "numpy"

    def __post_init__(self) -> None:
        """Validate backend contracts when a configuration is created."""
        if self.array != "numpy":
            raise ValueError("NumPy is the only supported array backend in ITAMAE 0.1.")
        if not isinstance(self.cosmology, CosmologyBackend):
            raise TypeError("cosmology must implement the CosmologyBackend protocol.")
        if not isinstance(self.units, UnitBackend):
            raise TypeError("units must implement the UnitBackend protocol.")

    @property
    def identifier(self) -> str:
        """Return a stable human-readable backend identifier."""
        return (
            f"array={self.array};cosmology={self.cosmology.identifier};"
            f"units={self.units.identifier}"
        )

    def metadata(self) -> Mapping[str, Any]:
        """Return immutable backend metadata for catalogs and cache keys."""
        return MappingProxyType(
            {
                "backend_identifier": self.identifier,
                "array_backend": self.array,
                "cosmology_backend": self.cosmology.identifier,
                "unit_backend": self.units.identifier,
            }
        )
