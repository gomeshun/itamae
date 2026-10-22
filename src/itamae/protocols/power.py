"""Structural interfaces for spectra and smoothing windows.

ITAMAE supplies numerical integration machinery while SASHIMI variants retain
ownership of their CDM, WDM, and FDM transfer functions and normalizations.
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class PowerSpectrum(Protocol):
    """Describe a positive power spectrum sampled as a function of wavenumber.

    Notes
    -----
    Power is a density-fluctuation spectrum P(k); a transfer modification supplies a power
    ratio, not an unsquared amplitude. The model owns normalization, k units, factors of h
    and any physical cutoff. Content identities distinguish changed tables even when labels
    are reused.
    """

    @property
    def identifier(self) -> str:
        """Return a stable identifier suitable for metadata and cache keys.

        Returns
        -------
        str
            Stable prescription identity for cache separation and catalogue provenance.

        Notes
        -----
        Power is a density-fluctuation spectrum P(k); a transfer modification supplies a power
        ratio, not an unsquared amplitude. The model owns normalization, k units, factors of h
        and any physical cutoff. Content identities distinguish changed tables even when labels
        are reused.
        """
        ...

    def __call__(self, wavenumber: Any) -> Any:
        """Evaluate the power spectrum at positive wavenumber.

        Parameters
        ----------
        wavenumber : float or array_like
            Positive Fourier wavenumber, in the inverse-length convention of the supplied
            spectrum; factors of h are owned by that spectrum.

        Returns
        -------
        array_like
            Nonnegative P(k) with the query shape and the provider-declared volume units.

        Notes
        -----
        Power is a density-fluctuation spectrum P(k); a transfer modification supplies a power
        ratio, not an unsquared amplitude. The model owns normalization, k units, factors of h
        and any physical cutoff. Content identities distinguish changed tables even when labels
        are reused.
        """
        ...


@runtime_checkable
class WindowFunction(Protocol):
    """Describe a dimensionless Fourier-space smoothing window.

    Notes
    -----
    Power is a density-fluctuation spectrum P(k); a transfer modification supplies a power
    ratio, not an unsquared amplitude. The model owns normalization, k units, factors of h
    and any physical cutoff. Content identities distinguish changed tables even when labels
    are reused.
    """

    @property
    def identifier(self) -> str:
        """Return a stable identifier suitable for metadata and cache keys.

        Returns
        -------
        str
            Stable prescription identity for cache separation and catalogue provenance.

        Notes
        -----
        Power is a density-fluctuation spectrum P(k); a transfer modification supplies a power
        ratio, not an unsquared amplitude. The model owns normalization, k units, factors of h
        and any physical cutoff. Content identities distinguish changed tables even when labels
        are reused.
        """
        ...

    def __call__(self, argument: Any) -> Any:
        """Evaluate the window at a nonnegative dimensionless argument.

        Parameters
        ----------
        argument : float or array_like
            Nonnegative dimensionless Fourier-window argument x=kR/c; output retains its shape.

        Returns
        -------
        array_like
            Dimensionless W(kR/c), with the query shape.

        Notes
        -----
        Power is a density-fluctuation spectrum P(k); a transfer modification supplies a power
        ratio, not an unsquared amplitude. The model owns normalization, k units, factors of h
        and any physical cutoff. Content identities distinguish changed tables even when labels
        are reused.
        """
        ...


__all__ = ["PowerSpectrum", "WindowFunction"]
