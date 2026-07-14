"""Validated one-dimensional grids used by population calculations."""

import numpy as np


def log_grid(minimum: float, maximum: float, size: int) -> np.ndarray:
    """Return a strictly increasing logarithmic grid.

    Parameters
    ----------
    minimum, maximum
        Positive finite endpoints with ``maximum > minimum``.
    size
        Number of nodes, including both endpoints.
    """
    if not np.isfinite(minimum) or not np.isfinite(maximum):
        raise ValueError("Logarithmic-grid endpoints must be finite.")
    if minimum <= 0.0 or maximum <= minimum:
        raise ValueError("Require 0 < minimum < maximum for a logarithmic grid.")
    if isinstance(size, bool) or not isinstance(size, (int, np.integer)) or size < 2:
        raise ValueError("Logarithmic-grid size must be an integer of at least two.")
    return np.geomspace(minimum, maximum, int(size))


def redshift_grid(
    target: float,
    maximum: float,
    step: float,
    *,
    include_maximum: bool = True,
) -> np.ndarray:
    """Return increasing accretion-redshift nodes above a target redshift.

    Parameters
    ----------
    target
        Target redshift, excluded from the returned grid.
    maximum
        Maximum accretion redshift.
    step
        Positive nominal spacing.
    include_maximum
        Append ``maximum`` when it is not already on the regular grid.

    Notes
    -----
    The final interval may be shorter than ``step`` when
    ``include_maximum=True``. This behavior makes the integration domain
    explicit and avoids silently extending beyond ``maximum``.
    """
    if not all(np.isfinite(value) for value in (target, maximum, step)):
        raise ValueError("Redshift-grid parameters must be finite.")
    if target < 0.0 or maximum <= target or step <= 0.0:
        raise ValueError("Require 0 <= target < maximum and step > 0.")
    count = int(np.floor((maximum - target) / step))
    grid = target + step * np.arange(1, count + 1, dtype=float)
    tolerance = 16.0 * np.finfo(float).eps * max(1.0, abs(maximum))
    if include_maximum and (grid.size == 0 or maximum - grid[-1] > tolerance):
        grid = np.append(grid, maximum)
    elif grid.size and grid[-1] > maximum + tolerance:
        grid = grid[grid <= maximum + tolerance]
    return grid


__all__ = ["log_grid", "redshift_grid"]
