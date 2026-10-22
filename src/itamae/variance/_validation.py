"""Internal validation shared by the two variance adapters."""

from contextlib import contextmanager

import numpy as np


def real_array(value, name):
    """Validate real numerical input before its use in a physical calculation.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.

    Returns
    -------
    numpy.ndarray
        Numerical array with the original shape; the helper rejects the invalid values
        checked in its implementation.

    Notes
    -----
    EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
    mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
    dS/dM. Validation rejects malformed inputs without clipping a physical variance or
    abundance.
    """
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name} must contain real numeric values.")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values.")
    return array


def coordinates(mass, redshift):
    """Validate coordinates before callbacks and obtain their common shape.

    Parameters
    ----------
    mass : float or array_like
        Positive halo mass. Use solar masses with the canonical backend; a supplied variance
        prescription must use the same mass convention and factors of h.
    redshift : float or array_like
        Dimensionless redshift. The cosmological background requires z > -1; population
        nodes use nonnegative redshifts.

    Returns
    -------
    tuple of numpy.ndarray
        Mass and redshift arrays broadcast to their common query shape.

    Notes
    -----
    EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
    mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
    dS/dM. Validation rejects malformed inputs without clipping a physical variance or
    abundance.
    """
    mass = real_array(mass, "Mass")
    redshift = real_array(redshift, "Redshift")
    if np.any(mass <= 0.0):
        raise ValueError("Masses must be finite and positive.")
    return np.broadcast_arrays(mass, redshift)


def aligned(value, shape, name, *, nonnegative=False, scalar_ok=False):
    """Validate a callback result without allowing an extra output dimension.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.
    shape : tuple of int
        Required population or coordinate shape. Extra axes cannot be silently treated as
        additional objects.
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.
    nonnegative : bool
        Require nonnegative values, as needed for power, variance or statistical weights.
    scalar_ok : bool
        Allow a constant scalar provider result to broadcast to the requested coordinate
        shape.

    Returns
    -------
    numpy.ndarray
        Finite callback result with exactly the requested shape, optionally broadcasting a
        scalar.

    Notes
    -----
    EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
    mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
    dS/dM. Validation rejects malformed inputs without clipping a physical variance or
    abundance.
    """
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
    """Turn arithmetic overflow/invalid operations into an explicit API failure.

    Parameters
    ----------
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.

    Yields
    ------
    None
        Context in which arithmetic overflow, invalid operations and division by zero become
        explicit ValueError failures.

    Notes
    -----
    EPS accretion uses S=sigma**2 and dS/dM. Keep the same transfer, growth, smoothing and
    mass-unit conventions in both quantities; a derivative of sigma cannot substitute for
    dS/dM. Validation rejects malformed inputs without clipping a physical variance or
    abundance.
    """
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            yield
    except (FloatingPointError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} exceeded the finite numerical domain: {exc}") from exc
