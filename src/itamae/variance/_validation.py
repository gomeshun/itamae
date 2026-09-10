"""Internal validation shared by the two variance adapters."""

from contextlib import contextmanager

import numpy as np


def real_array(value, name):
    """Reject non-real/non-numeric values before conversion, then require finiteness."""
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values.")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values.")
    return array


def coordinates(mass, redshift):
    """Validate coordinates before callbacks and obtain their common shape."""
    mass = real_array(mass, "Mass")
    redshift = real_array(redshift, "Redshift")
    if np.any(mass <= 0.0):
        raise ValueError("Masses must be finite and positive.")
    return np.broadcast_arrays(mass, redshift)


def aligned(value, shape, name, *, nonnegative=False, scalar_ok=False):
    """Validate a callback result without allowing an extra output dimension."""
    array = real_array(value, name)
    if scalar_ok and array.ndim == 0:
        array = np.broadcast_to(array, shape)
    if array.shape != shape:
        raise ValueError(f"{name} must return shape {shape}; got {array.shape}.")
    if nonnegative and np.any(array < 0.0):
        raise ValueError(f"{name} must return nonnegative values.")
    return array


@contextmanager
def numerical_errors(name):
    """Turn arithmetic overflow/invalid operations into an explicit API failure."""
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            yield
    except (FloatingPointError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} exceeded the finite numerical domain: {exc}") from exc
