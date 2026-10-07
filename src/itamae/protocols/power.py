"""Structural interfaces for spectra and smoothing windows.

ITAMAE supplies numerical integration machinery while SASHIMI variants retain
ownership of their CDM, WDM, and FDM transfer functions and normalizations.
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class PowerSpectrum(Protocol):
    """Describe a positive power spectrum sampled as a function of wavenumber."""

    @property
    def identifier(self) -> str:
        """Return a stable identifier suitable for metadata and cache keys."""
        ...

    def __call__(self, wavenumber: Any) -> Any:
        """Evaluate the power spectrum at positive wavenumber."""
        ...


@runtime_checkable
class WindowFunction(Protocol):
    """Describe a dimensionless Fourier-space smoothing window."""

    @property
    def identifier(self) -> str:
        """Return a stable identifier suitable for metadata and cache keys."""
        ...

    def __call__(self, argument: Any) -> Any:
        """Evaluate the window at a nonnegative dimensionless argument."""
        ...


__all__ = ["PowerSpectrum", "WindowFunction"]
