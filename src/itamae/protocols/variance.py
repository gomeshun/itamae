"""Variance-model interface for halo population calculations.

A variance model supplies the smoothed linear-density variance used by EPS-like
accretion prescriptions. The protocol does not prescribe a power spectrum,
transfer function, window function, or growth convention; those choices remain
explicit properties of the concrete SASHIMI model.
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class VarianceModel(Protocol):
    r"""Describe a mass-variance implementation consumed by ITAMAE pipelines.

    Attributes
    ----------
    identifier : str
        Stable description of the physical and numerical variance prescription.
        It should include enough information to distinguish transfer functions,
        window functions, normalization choices, and relevant backend settings.

    Notes
    -----
    Mass inputs use the active ITAMAE unit convention. Implementations may accept
    scalar or array inputs, but must document broadcasting and factors of ``h``.
    ``dvariance_dmass`` denotes the derivative of :math:`S=\sigma^2`, not the
    derivative of :math:`\sigma`.

    EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
    mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
    dS/dM. Validation rejects malformed inputs without clipping a physical variance or
    abundance.
    """

    @property
    def identifier(self) -> str:
        """Return a stable model identifier.

        Returns
        -------
        str
            Stable prescription identity for cache separation and catalogue provenance.

        Notes
        -----
        EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
        mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
        dS/dM. Validation rejects malformed inputs without clipping a physical variance or
        abundance.
        """
        ...

    def sigma(self, mass: Any, z: Any = 0.0) -> Any:
        r"""Return the linear rms fluctuation :math:`\sigma(M,z)`.

        Parameters
        ----------
        mass : float or array_like
            Positive halo mass. Use solar masses with the canonical backend; a supplied variance
            prescription must use the same mass convention and factors of h.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        numpy.ndarray
            Dimensionless rms density fluctuation sigma(M,z), with broadcast mass/redshift
            shape.

        Notes
        -----
        EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
        mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
        dS/dM. Validation rejects malformed inputs without clipping a physical variance or
        abundance.
        """
        ...

    def variance(self, mass: Any, z: Any = 0.0) -> Any:
        r"""Return the linear variance :math:`S(M,z)=\sigma(M,z)^2`.

        Parameters
        ----------
        mass : float or array_like
            Positive halo mass. Use solar masses with the canonical backend; a supplied variance
            prescription must use the same mass convention and factors of h.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        numpy.ndarray
            Dimensionless linear variance S(M,z)=sigma(M,z)**2, with broadcast mass/redshift
            shape.

        Notes
        -----
        EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
        mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
        dS/dM. Validation rejects malformed inputs without clipping a physical variance or
        abundance.
        """
        ...

    def dvariance_dmass(self, mass: Any, z: Any = 0.0) -> Any:
        r"""Return :math:`\mathrm{d}S/\mathrm{d}M`.

        Parameters
        ----------
        mass : float or array_like
            Positive halo mass. Use solar masses with the canonical backend; a supplied variance
            prescription must use the same mass convention and factors of h.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        numpy.ndarray
            Derivative dS/dM, in inverse mass units, with broadcast mass/redshift shape.

        Notes
        -----
        EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
        mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
        dS/dM. Validation rejects malformed inputs without clipping a physical variance or
        abundance.
        """
        ...
