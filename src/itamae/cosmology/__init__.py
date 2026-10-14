"""Cosmology backends."""

from typing import TYPE_CHECKING, Any

from .native import NativeFlatLCDM

if TYPE_CHECKING:
    from .colossus import ColossusCosmology


def __getattr__(name: str) -> Any:
    """Load a requested optional cosmology backend by its public name.

    Parameters
    ----------
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.

    Returns
    -------
    type
        Requested backend class; unknown names raise AttributeError.

    Notes
    -----
    The optional Colossus backend is loaded only when requested, so a native-background
    calculation does not depend on importing Colossus.
    """
    if name == "ColossusCosmology":
        from .colossus import ColossusCosmology

        return ColossusCosmology
    raise AttributeError(name)


__all__ = ["NativeFlatLCDM", "ColossusCosmology"]
