"""Tabulated and transfer-modified power-spectrum mechanisms."""

from __future__ import annotations

from collections.abc import Callable
from hashlib import sha256
from typing import Any

import numpy as np
from scipy.interpolate import interp1d

from itamae.protocols import PowerSpectrum


def _array_digest(*arrays: np.ndarray) -> str:
    """Return a platform-independent digest for floating-point table data."""
    digest = sha256()
    for array in arrays:
        canonical = np.ascontiguousarray(array, dtype="<f8")
        digest.update(np.asarray(canonical.shape, dtype="<i8").tobytes())
        digest.update(canonical.tobytes())
    return digest.hexdigest()


class TabulatedPowerSpectrum:
    """Interpolate a strictly positive spectrum in logarithmic coordinates.

    Parameters
    ----------
    wavenumber, power
        Aligned one-dimensional arrays. Wavenumbers must be finite, positive,
        and strictly increasing. Power values must be finite and nonnegative.
    identifier
        Optional provenance identifier. When omitted, a SHA-256 digest of the
        canonicalized table is used.
    interpolation
        ``"log-log"`` for a strictly positive smooth spectrum or ``"linear"``
        for a nonnegative spectrum that may contain physical zeros, such as an
        oscillatory FDM transfer table.
    extrapolate
        If false, evaluation outside the tabulated domain raises
        :class:`ValueError`. If true, the two endpoint slopes are extended in
        the configured interpolation coordinates.

    Notes
    -----
    ITAMAE does not attach units or factors of ``h`` to the table. The SASHIMI
    variant supplying it must document those conventions in ``identifier``.
    """

    def __init__(
        self,
        wavenumber: Any,
        power: Any,
        *,
        identifier: str | None = None,
        interpolation: str = "log-log",
        extrapolate: bool = False,
    ) -> None:
        k = np.asarray(wavenumber, dtype=float)
        values = np.asarray(power, dtype=float)
        if k.ndim != 1 or values.ndim != 1 or k.shape != values.shape or k.size < 2:
            raise ValueError("Power-spectrum arrays must be aligned one-dimensional tables.")
        if not np.all(np.isfinite(k)) or np.any(k <= 0.0) or np.any(np.diff(k) <= 0.0):
            raise ValueError("Wavenumbers must be finite, positive, and strictly increasing.")
        if interpolation not in {"log-log", "linear"}:
            raise ValueError("interpolation must be 'log-log' or 'linear'.")
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("Power-spectrum values must be finite and nonnegative.")
        if interpolation == "log-log" and np.any(values == 0.0):
            raise ValueError("Log-log interpolation requires strictly positive power.")
        if identifier is not None and (not isinstance(identifier, str) or not identifier.strip()):
            raise ValueError("identifier must be None or a non-empty string.")

        self._wavenumber = k.copy()
        self._power = values.copy()
        source_identifier = identifier or f"sha256={_array_digest(k, values)}"
        self._identifier = (
            f"tabulated-power:interpolation={interpolation};source=({source_identifier})"
        )
        self._interpolation = interpolation
        self._extrapolate = bool(extrapolate)
        interpolation_x = np.log(k) if interpolation == "log-log" else k
        interpolation_y = np.log(values) if interpolation == "log-log" else values
        self._interpolator = interp1d(
            interpolation_x,
            interpolation_y,
            kind="linear",
            bounds_error=not self._extrapolate,
            fill_value="extrapolate" if self._extrapolate else np.nan,
            assume_sorted=True,
        )

    @property
    def identifier(self) -> str:
        """Return the user-supplied provenance or content-derived identifier."""
        return self._identifier

    @property
    def domain(self) -> tuple[float, float]:
        """Return the inclusive tabulated wavenumber interval."""
        return float(self._wavenumber[0]), float(self._wavenumber[-1])

    def __call__(self, wavenumber: Any) -> np.ndarray:
        """Evaluate the configured interpolant."""
        k = np.asarray(wavenumber, dtype=float)
        if not np.all(np.isfinite(k)) or np.any(k <= 0.0):
            raise ValueError("Evaluation wavenumbers must be finite and positive.")
        if not self._extrapolate and (
            np.any(k < self._wavenumber[0]) or np.any(k > self._wavenumber[-1])
        ):
            raise ValueError(f"Wavenumber lies outside the tabulated domain {self.domain}.")
        if self._interpolation == "log-log":
            return np.exp(self._interpolator(np.log(k)))
        return np.asarray(self._interpolator(k), dtype=float)


class TransferModifiedPowerSpectrum:
    """Multiply a base spectrum by a model-supplied power ratio.

    Parameters
    ----------
    base
        Spectrum satisfying :class:`itamae.protocols.PowerSpectrum`.
    power_ratio
        Callable returning ``P_modified(k) / P_base(k)``. This is explicitly a
        power ratio, not a transfer amplitude, so ITAMAE never guesses whether
        it should be squared.
    ratio_identifier
        Stable identifier owned by the SASHIMI variant.
    """

    def __init__(
        self,
        base: PowerSpectrum,
        power_ratio: Callable[[Any], Any],
        *,
        ratio_identifier: str,
    ) -> None:
        if not isinstance(base, PowerSpectrum):
            raise TypeError("base must implement the ITAMAE power-spectrum protocol.")
        if not callable(power_ratio):
            raise TypeError("power_ratio must be callable.")
        if not isinstance(ratio_identifier, str) or not ratio_identifier.strip():
            raise ValueError("ratio_identifier must be a non-empty string.")
        self._base = base
        self._power_ratio = power_ratio
        self._identifier = f"transfer-modified:base=({base.identifier});ratio=({ratio_identifier})"

    @property
    def identifier(self) -> str:
        """Return the composed base-spectrum and power-ratio identifier."""
        return self._identifier

    def __call__(self, wavenumber: Any) -> np.ndarray:
        """Return the base spectrum multiplied by the supplied power ratio."""
        ratio = np.asarray(self._power_ratio(wavenumber), dtype=float)
        if not np.all(np.isfinite(ratio)) or np.any(ratio < 0.0):
            raise ValueError("The model-supplied power ratio must be finite and nonnegative.")
        return np.asarray(self._base(wavenumber), dtype=float) * ratio


__all__ = ["TabulatedPowerSpectrum", "TransferModifiedPowerSpectrum"]
