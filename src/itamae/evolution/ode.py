"""Generic ODE integration for model-supplied evolution equations."""

import warnings
from collections.abc import Callable

import numpy as np
from scipy.integrate import ODEintWarning, odeint, solve_ivp


def _real_finite(value, name):
    array = np.asarray(value)
    if np.iscomplexobj(array):
        raise ValueError(f"{name} must be real and finite.")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be real and finite.")
    return array


def solve_evolution(rhs: Callable, y0, t_eval, *, args=(), method="RK45", rtol=None, atol=None):
    """Integrate a real evolution equation on a strictly monotonic output grid.

    ``rhs(t, y, *args)`` returns one derivative for each initial state. Output
    has shape ``(len(t_eval), len(y0))``. Scalar initial states are accepted.

    ``method='RK45'`` preserves ITAMAE's existing solve_ivp defaults
    (rtol=1e-8, atol=1e-10). ``method='odeint'`` explicitly selects ODEPACK
    LSODA with SciPy's original odeint tolerance defaults when None is passed.
    Other solve_ivp method names are forwarded explicitly. Tolerances must
    be finite positive scalars or arrays matching y0. No failure or non-finite
    callback/output is silently returned as a successful evolution.
    """
    grid = _real_finite(t_eval, "t_eval")
    initial = np.atleast_1d(_real_finite(y0, "y0"))
    if initial.ndim != 1 or initial.size == 0:
        raise ValueError("y0 must be a nonempty one-dimensional state vector.")
    if grid.ndim != 1 or grid.size < 2:
        raise ValueError("t_eval must be a one-dimensional grid with at least two points.")
    delta = np.diff(grid)
    if not (np.all(delta > 0.0) or np.all(delta < 0.0)):
        raise ValueError("t_eval must be strictly monotonic.")
    for name, value in (("rtol", rtol), ("atol", atol)):
        if value is not None:
            tolerance = _real_finite(value, name)
            if tolerance.shape not in ((), initial.shape) or np.any(tolerance <= 0):
                raise ValueError(f"{name} must be positive and scalar or match y0.")

    def evaluate(time, state):
        derivative = _real_finite(rhs(time, state, *args), "rhs")
        if derivative.shape != initial.shape:
            raise ValueError(f"rhs must return shape {initial.shape}; received {derivative.shape}.")
        return derivative

    if method == "odeint":
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", ODEintWarning)
                values, diagnostics = odeint(
                    evaluate, initial, grid, tfirst=True, full_output=True, rtol=rtol, atol=atol
                )
        except ODEintWarning as error:
            raise RuntimeError(f"odeint evolution failed: {error}") from error
        if diagnostics["message"] != "Integration successful.":
            raise RuntimeError(f"odeint evolution failed: {diagnostics['message']}")
    else:
        result = solve_ivp(
            evaluate,
            (float(grid[0]), float(grid[-1])),
            initial,
            t_eval=grid,
            method=method,
            rtol=1e-8 if rtol is None else rtol,
            atol=1e-10 if atol is None else atol,
        )
        if not result.success:
            raise RuntimeError(result.message)
        values = result.y.T
    values = _real_finite(values, "evolution output")
    if values.shape != (grid.size, initial.size):
        raise RuntimeError("ODE solver returned an incomplete output grid.")
    return values
