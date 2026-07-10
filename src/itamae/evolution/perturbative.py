"""Perturbative sequence acceleration utilities."""

import numpy as np


def shanks_transform(s0, s1, s2, *, tolerance: float = 1e-14):
    """Apply the three-term Shanks transformation.

    Parameters
    ----------
    s0, s1, s2
        Three consecutive partial approximations.
    tolerance
        Minimum absolute denominator used to accept the transformed value.

    Returns
    -------
    numpy.ndarray
        Accelerated estimate. ``s2`` is returned where the denominator is too
        small for a stable transformation.
    """
    s0, s1, s2 = np.broadcast_arrays(
        np.asarray(s0, dtype=float), np.asarray(s1, dtype=float), np.asarray(s2, dtype=float)
    )
    denominator = s2 - 2.0 * s1 + s0
    safe = np.abs(denominator) > tolerance
    return np.where(safe, s2 - (s2 - s1) ** 2 / denominator, s2)
