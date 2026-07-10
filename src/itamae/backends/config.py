"""Immutable backend selection for ITAMAE calculations."""

from dataclasses import dataclass
from typing import Any


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

    cosmology: Any
    units: Any
    array: str = "numpy"

    @property
    def identifier(self) -> str:
        """Return a stable human-readable backend identifier."""
        return f"{self.array}:{self.cosmology.identifier}:{self.units.identifier}"
