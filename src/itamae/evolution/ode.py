"""Generic ODE integration for model-supplied evolution equations."""

from collections.abc import Callable

import numpy as np
from scipy.integrate import solve_ivp


def solve_evolution(rhs: Callable, y0, t_eval, *, args=(), rtol=1e-8, atol=1e-10):
    """Integrate an evolution equation on an explicit output grid.

    Parameters
    ----------
    rhs
        Callable with signature ``rhs(t, y, *args)``.
    y0
        Initial state vector.
    t_eval
        Strictly increasing or decreasing output grid.
    args
        Additional positional arguments forwarded to ``rhs``.
    rtol, atol
        Relative and absolute error tolerances.

    Returns
    -------
    numpy.ndarray
        State values with the time axis first.
    """
    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or t_eval.size < 2:
        raise ValueError("t_eval must be a one-dimensional grid with at least two points.")
    delta = np.diff(t_eval)
    if not (np.all(delta > 0.0) or np.all(delta < 0.0)):
        raise ValueError("t_eval must be strictly monotonic.")
    result = solve_ivp(
        lambda t, y: rhs(t, y, *args),
        (float(t_eval[0]), float(t_eval[-1])),
        np.atleast_1d(np.asarray(y0, dtype=float)),
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.y.T
