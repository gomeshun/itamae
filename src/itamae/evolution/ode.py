"""Generic ODE integration for model-supplied evolution equations."""

import warnings
from collections.abc import Callable, Mapping

import numpy as np
from scipy.integrate import ODEintWarning, odeint, solve_ivp


def _real_finite(value, name):
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
    The controller checks callback shape and finiteness and treats integration failure as an
    error. Solver tolerances concern the supplied differential equation; they do not
    validate its physical prescription or its interpolation inputs.
    """
    array = np.asarray(value)
    if np.iscomplexobj(array):
        raise ValueError(f"{name} must be real and finite.")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be real and finite.")
    return array


def _validated_odeint_options(options, initial, args):
    """Validate LSODA controls before integrating a physical state.

    Parameters
    ----------
    options : Mapping, optional
        Explicit numerical controls for the named solver or table. Unsupported or reserved
        keys are rejected by the receiving routine.
    initial : Mapping or array_like
        Per-node physical state for the indicated evolution stage; its fields and units are
        owned by the supplied physical model.
    args : tuple
        Additional positional parameters forwarded unchanged to the supplied numerical or
        physical callback.

    Returns
    -------
    dict
        Checked options, including a time-first Jacobian wrapper when requested.

    Notes
    -----
    The controller checks callback shape and finiteness and treats integration failure as an
    error. Solver tolerances concern the supplied differential equation; they do not
    validate its physical prescription or its interpolation inputs.
    """
    if options is None:
        return {}
    if not isinstance(options, Mapping):
        raise ValueError("odeint_options must be a mapping.")
    values = dict(options)
    allowed = {
        "Dfun",
        "col_deriv",
        "ml",
        "mu",
        "tcrit",
        "h0",
        "hmax",
        "hmin",
        "ixpr",
        "mxstep",
        "mxhnil",
        "mxordn",
        "mxords",
        "printmessg",
    }
    unknown = values.keys() - allowed
    if unknown:
        raise ValueError(f"Unsupported or reserved odeint options: {sorted(map(str, unknown))}.")
    for name in ("h0", "hmax", "hmin"):
        if name in values:
            value = _real_finite(values[name], name)
            if value.ndim != 0 or (name != "h0" and value < 0):
                raise ValueError(f"{name} must be a finite scalar with a valid step size.")
    for name in ("mxstep", "mxhnil", "mxordn", "mxords", "ml", "mu"):
        if name in values:
            value = values[name]
            if value is None and name in {"ml", "mu"}:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer.")
            upper = {"mxordn": 12, "mxords": 5, "ml": initial.size - 1, "mu": initial.size - 1}.get(
                name
            )
            if upper is not None and value > upper:
                raise ValueError(f"{name} exceeds the supported maximum {upper}.")
    for name in ("col_deriv", "ixpr", "printmessg"):
        if name in values:
            value = values[name]
            if not isinstance(value, (bool, np.bool_, int, np.integer)) or value not in (0, 1):
                raise ValueError(f"{name} must be a boolean or 0/1.")
    if "tcrit" in values:
        critical = _real_finite(values["tcrit"], "tcrit")
        if critical.ndim != 1:
            raise ValueError("tcrit must be a one-dimensional finite array.")
        values["tcrit"] = critical
    jacobian = values.get("Dfun")
    if jacobian is not None:
        if not callable(jacobian):
            raise ValueError("Dfun must be callable.")
        banded = values.get("ml") is not None or values.get("mu") is not None
        rows = (values.get("ml") or 0) + (values.get("mu") or 0) + 1 if banded else initial.size
        shape = (rows, initial.size)
        if values.get("col_deriv", False):
            shape = shape[::-1]

        def evaluate_jacobian(time, state):
            """Check the supplied Jacobian against the selected dense or banded state layout.

            Parameters
            ----------
            time : float or array_like
                Evolution coordinate or cosmic time, in the units of the supplied equation;
                canonical state times are in Gyr.
            state : Mapping or array_like
                Per-node physical state for the indicated evolution stage; its fields and units are
                owned by the supplied physical model.

            Returns
            -------
            numpy.ndarray
                Finite Jacobian with the exact requested shape and derivative orientation.

            Notes
            -----
            The controller checks callback shape and finiteness and treats integration failure as an
            error. Solver tolerances concern the supplied differential equation; they do not
            validate its physical prescription or its interpolation inputs.
            """
            result = _real_finite(jacobian(time, state, *args), "Dfun")
            if result.shape != shape:
                raise ValueError(f"Dfun must return shape {shape}; received {result.shape}.")
            return result

        values["Dfun"] = evaluate_jacobian
    return values


def solve_evolution(
    rhs: Callable,
    y0,
    t_eval,
    *,
    args=(),
    method="RK45",
    rtol=None,
    atol=None,
    odeint_options=None,
    allow_repeated_times=False,
):
    """Integrate a real evolution equation on a strictly monotonic output grid.

    ``rhs(t, y, *args)`` returns one derivative for each initial state. Output
    has shape ``(len(t_eval), len(y0))``. Scalar initial states are accepted.

    ``method='RK45'`` preserves ITAMAE's existing solve_ivp defaults
    (rtol=1e-8, atol=1e-10). ``method='odeint'`` explicitly selects ODEPACK
    LSODA with SciPy's original odeint tolerance defaults when None is passed.
    Other solve_ivp method names are forwarded explicitly. Tolerances must
    be finite positive scalars or arrays matching y0. No failure or non-finite
    callback/output is silently returned as a successful evolution.

    ``odeint_options`` forwards explicit LSODA step/order limits, critical
    times and an optional time-first Jacobian ``Dfun(t, y, *args)``. The
    controller owns ``args``, ``tfirst``, ``full_output`` and tolerances;
    duplicate/unknown options fail before the physical callback. Options for
    this backend cannot be silently passed to another solver.

    Explicit ``allow_repeated_times=True`` preserves odeint's repeated output
    coordinates, including a constant grid with no evolution. Integration
    runs only on distinct coordinates and results retain requested order.

    Parameters
    ----------
    rhs : callable
        Time-first derivative rhs(t,y,*args) in state units per integration-coordinate unit.
        The generic controller does not choose a tidal law.
    y0 : float or array_like
        Finite initial state at the first integration coordinate; its units must match the
        derivative convention.
    t_eval : array_like
        One-dimensional increasing or decreasing output grid. The first coordinate is the
        initial condition; redshift evolution usually decreases.
    args : tuple
        Additional positional parameters forwarded unchanged to the supplied numerical or
        physical callback.
    method : str
        Explicit SciPy evolution method; odeint uses LSODA, while other names are passed to
        solve_ivp.
    rtol : float or array_like, optional
        Positive relative or absolute local ODE error tolerance. Absolute tolerance uses the
        integrated state units; it is not a physical-model uncertainty.
    atol : float or array_like, optional
        Positive relative or absolute local ODE error tolerance. Absolute tolerance uses the
        integrated state units; it is not a physical-model uncertainty.
    odeint_options : Mapping, optional
        Explicit numerical controls for the named solver or table. Unsupported or reserved
        keys are rejected by the receiving routine.
    allow_repeated_times : bool
        Permit repeated output coordinates only for the odeint backend, retaining the
        requested ordering.

    Returns
    -------
    numpy.ndarray
        Trajectory of shape (len(t_eval), len(y0)); scalar initial states give one state
        column.

    Notes
    -----
    The controller checks callback shape and finiteness and treats integration failure as an
    error. Solver tolerances concern the supplied differential equation; they do not
    validate its physical prescription or its interpolation inputs.
    """
    grid = _real_finite(t_eval, "t_eval")
    initial = np.atleast_1d(_real_finite(y0, "y0"))
    if initial.ndim != 1 or initial.size == 0:
        raise ValueError("y0 must be a nonempty one-dimensional state vector.")
    if not isinstance(allow_repeated_times, (bool, np.bool_)):
        raise ValueError("allow_repeated_times must be boolean.")
    if allow_repeated_times and method != "odeint":
        raise ValueError("allow_repeated_times requires method='odeint'.")
    minimum_points = 1 if allow_repeated_times else 2
    if grid.ndim != 1 or grid.size < minimum_points:
        raise ValueError(f"t_eval must be one-dimensional with at least {minimum_points} points.")
    delta = np.diff(grid)
    monotonic = (
        (np.all(delta >= 0) or np.all(delta <= 0))
        if allow_repeated_times
        else (np.all(delta > 0) or np.all(delta < 0))
    )
    if not monotonic:
        raise ValueError("t_eval must be strictly monotonic.")
    restore = None
    if allow_repeated_times:
        distinct = np.concatenate(([True], delta != 0))
        restore = np.cumsum(distinct) - 1
        grid = grid[distinct]
    for name, value in (("rtol", rtol), ("atol", atol)):
        if value is not None:
            tolerance = _real_finite(value, name)
            if tolerance.shape not in ((), initial.shape) or np.any(tolerance <= 0):
                raise ValueError(f"{name} must be positive and scalar or match y0.")

    def evaluate(time, state):
        """Evaluate the physical derivative without admitting extra state axes or non-finite
        values.

        Parameters
        ----------
        time : float or array_like
            Evolution coordinate or cosmic time, in the units of the supplied equation;
            canonical state times are in Gyr.
        state : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.

        Returns
        -------
        numpy.ndarray
            Finite derivative vector matching the initial state shape.

        Notes
        -----
        The controller checks callback shape and finiteness and treats integration failure as an
        error. Solver tolerances concern the supplied differential equation; they do not
        validate its physical prescription or its interpolation inputs.
        """
        derivative = _real_finite(rhs(time, state, *args), "rhs")
        if derivative.shape != initial.shape:
            raise ValueError(f"rhs must return shape {initial.shape}; received {derivative.shape}.")
        return derivative

    if method == "odeint":
        options = _validated_odeint_options(odeint_options, initial, args)
        if grid.size == 1:
            assert restore is not None
            return np.broadcast_to(initial, (restore.size, initial.size)).copy()
        with warnings.catch_warnings(record=True) as observed:
            warnings.simplefilter("always", ODEintWarning)
            values, diagnostics = odeint(
                evaluate,
                initial,
                grid,
                tfirst=True,
                full_output=True,
                rtol=rtol,
                atol=atol,
                **options,
            )
        if diagnostics["message"] != "Integration successful.":
            raise RuntimeError(f"odeint evolution failed: {diagnostics['message']}")
        for warning in observed:
            if (
                issubclass(warning.category, ODEintWarning)
                and str(warning.message) != "Integration successful."
            ):
                cause = warning.message if isinstance(warning.message, BaseException) else None
                raise RuntimeError(f"odeint evolution failed: {warning.message}") from cause
            warnings.warn(warning.message, warning.category, stacklevel=2)
    else:
        if odeint_options is not None:
            raise ValueError("odeint_options require method='odeint'.")
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
    return values if restore is None else values[restore]
