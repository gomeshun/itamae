"""Dimensionless Fourier-space smoothing windows."""

from dataclasses import dataclass
from typing import Any

import numpy as np


def _validated_argument(argument: Any) -> np.ndarray:
    """Return a finite nonnegative dimensionless window argument."""
    value = np.asarray(argument, dtype=float)
    if not np.all(np.isfinite(value)) or np.any(value < 0.0):
        raise ValueError("Window arguments must be finite and nonnegative.")
    return value


@dataclass(frozen=True, slots=True)
class SphericalTopHatWindow:
    r"""Evaluate the real-space spherical top-hat Fourier transform.

    Notes
    -----
    The implementation uses a small-argument series for
    :math:`3(\sin x-x\cos x)/x^3` to avoid cancellation.
    """

    @property
    def identifier(self) -> str:
        """Return the stable window identifier."""
        return "window:spherical-top-hat:v1"

    def __call__(self, argument: Any) -> np.ndarray:
        """Evaluate the window with scalar/array broadcasting."""
        x = _validated_argument(argument)
        result = np.empty_like(x, dtype=float)
        small = x < 1.0e-3
        x2 = x[small] * x[small]
        result[small] = 1.0 - x2 / 10.0 + x2 * x2 / 280.0
        regular = ~small
        xr = x[regular]
        result[regular] = 3.0 * (np.sin(xr) - xr * np.cos(xr)) / xr**3
        return result


@dataclass(frozen=True, slots=True)
class SharpKWindow:
    """Evaluate an ideal sharp cutoff in Fourier space."""

    @property
    def identifier(self) -> str:
        """Return the stable window identifier."""
        return "window:sharp-k:v1"

    def __call__(self, argument: Any) -> np.ndarray:
        """Return one at and below the cutoff and zero above it."""
        return (_validated_argument(argument) <= 1.0).astype(float)


__all__ = ["SharpKWindow", "SphericalTopHatWindow"]
