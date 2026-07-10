"""Cosmology backends."""

from typing import TYPE_CHECKING, Any

from .native import NativeFlatLCDM

if TYPE_CHECKING:
    from .colossus import ColossusCosmology


def __getattr__(name: str) -> Any:
    """Load optional cosmology backends lazily."""
    if name == "ColossusCosmology":
        from .colossus import ColossusCosmology

        return ColossusCosmology
    raise AttributeError(name)


__all__ = ["NativeFlatLCDM", "ColossusCosmology"]
