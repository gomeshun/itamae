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
    """

    @property
    def identifier(self) -> str:
        """Return a stable model identifier."""
        ...

    def sigma(self, mass: Any, z: Any = 0.0) -> Any:
        r"""Return the linear rms fluctuation :math:`\sigma(M,z)`."""
        ...

    def variance(self, mass: Any, z: Any = 0.0) -> Any:
        r"""Return the linear variance :math:`S(M,z)=\sigma(M,z)^2`."""
        ...

    def dvariance_dmass(self, mass: Any, z: Any = 0.0) -> Any:
        r"""Return :math:`\mathrm{d}S/\mathrm{d}M`."""
        ...
