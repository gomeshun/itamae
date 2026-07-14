"""Perturbative Picard evolution and sequence acceleration utilities."""

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import cumulative_trapezoid


def shanks_transform(s0: Any, s1: Any, s2: Any, *, tolerance: float = 1e-14) -> np.ndarray:
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
    if tolerance < 0.0 or not np.isfinite(tolerance):
        raise ValueError("Shanks tolerance must be finite and nonnegative.")
    s0, s1, s2 = np.broadcast_arrays(
        np.asarray(s0, dtype=float), np.asarray(s1, dtype=float), np.asarray(s2, dtype=float)
    )
    denominator = s2 - 2.0 * s1 + s0
    safe = np.abs(denominator) > tolerance
    with np.errstate(divide="ignore", invalid="ignore"):
        accelerated = s2 - (s2 - s1) ** 2 / denominator
    return np.where(safe, accelerated, s2)


@dataclass(frozen=True, slots=True)
class PerturbativeEvolutionSolver:
    """Solve a scalar evolution equation by successive Picard iteration.

    Parameters
    ----------
    order
        Highest perturbative correction order. ``order=0`` performs one
        integration of the right-hand side evaluated at the initial state.
    shanks
        Apply elementwise Shanks acceleration to the final three iterates.

    Notes
    -----
    The supplied right-hand side must accept the full time grid and current
    trajectory, returning values broadcastable to the grid shape. Calibrated
    mass-loss equations remain owned by individual SASHIMI packages.
    """

    order: int = 2
    shanks: bool = False
    tolerance: float = 1e-14

    def __post_init__(self) -> None:
        """Validate solver configuration."""
        if isinstance(self.order, bool) or not isinstance(self.order, (int, np.integer)):
            raise TypeError("Perturbative order must be an integer.")
        if self.order < 0:
            raise ValueError("Perturbative order must be nonnegative.")
        if self.shanks and self.order < 1:
            raise ValueError("Shanks acceleration requires order >= 1.")
        if self.tolerance < 0.0 or not np.isfinite(self.tolerance):
            raise ValueError("Solver tolerance must be finite and nonnegative.")

    def solve(
        self,
        rhs: Callable[..., Any],
        y0: float,
        t_grid: Any,
        *,
        args: tuple[Any, ...] = (),
    ) -> np.ndarray:
        """Return the perturbative trajectory on an explicit time grid.

        Parameters
        ----------
        rhs
            Vectorized callable with signature ``rhs(t_grid, trajectory,
            *args)``.
        y0
            Scalar initial state.
        t_grid
            Strictly increasing or decreasing one-dimensional grid.
        args
            Additional positional arguments forwarded to ``rhs``.
        """
        time = np.asarray(t_grid, dtype=float)
        if time.ndim != 1 or time.size < 2:
            raise ValueError("t_grid must be one-dimensional with at least two points.")
        if not np.all(np.isfinite(time)):
            raise ValueError("t_grid must contain only finite values.")
        delta = np.diff(time)
        if not (np.all(delta > 0.0) or np.all(delta < 0.0)):
            raise ValueError("t_grid must be strictly monotonic.")
        initial = float(y0)
        if not np.isfinite(initial):
            raise ValueError("y0 must be finite.")

        approximations = [np.full(time.shape, initial, dtype=float)]
        for _ in range(self.order + 1):
            derivative = np.asarray(rhs(time, approximations[-1], *args), dtype=float)
            try:
                derivative = np.broadcast_to(derivative, time.shape)
            except ValueError as exc:
                raise ValueError("rhs output must be broadcastable to t_grid.") from exc
            if not np.all(np.isfinite(derivative)):
                raise ValueError("rhs returned non-finite derivatives.")
            correction = cumulative_trapezoid(derivative, time, initial=0.0)
            approximations.append(initial + correction)

        if self.shanks and len(approximations) >= 3:
            return shanks_transform(*approximations[-3:], tolerance=self.tolerance)
        return approximations[-1]


__all__ = ["PerturbativeEvolutionSolver", "shanks_transform"]
