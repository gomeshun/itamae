"""Preserve SciPy integration controls while rejecting ambiguous overrides."""

import numpy as np
import pytest
from scipy.integrate import ODEintWarning, odeint
from itamae.evolution import solve_evolution


@pytest.mark.parametrize("reverse", [False, True])
def test_selected_controls_preserve_direct_odeint_values(reverse):
    """Apply the caller's time-step and work limits to the actual physical integration.

    Parameters
    ----------
    reverse : bool
        Whether to reverse the integration grid while retaining consistent signed initial
        step.

    Notes
    -----
    Pass criteria: Forward and reversed forced-decay solutions exactly match direct SciPy
    with the same hmax, h0, mxstep and tolerances, leaving the options unchanged.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    time = np.linspace(0.0, 2.0, 23)
    if reverse:
        time = time[::-1]

    def rhs(t, y, rate):
        """Provide a forced linear decay equation with an independently called SciPy reference.

        Parameters
        ----------
        t : float
            Dimensionless integration time of the analytic test equation.
        y : numpy.ndarray
            Current state of the analytic differential equation, in the fixture convention.
        rate : float
            Positive linear decay coefficient in inverse fixture-time units.

        Returns
        -------
        numpy.ndarray
            -rate*y+sin(t), the controlled state derivative.

        Notes
        -----
        Apply the caller's time-step and work limits to the actual physical integration. Pass
        criteria: Forward and reversed forced-decay solutions exactly match direct SciPy with
        the same hmax, h0, mxstep and tolerances, leaving the options unchanged. Failure occurs
        when the stated comparison or expected rejection is not met. This controlled calculation
        does not calibrate a subhalo population against simulations or observations.
        """
        return -rate * y + np.sin(t)

    options = {"hmax": 0.03, "h0": -0.001 if reverse else 0.001, "mxstep": 5000}
    expected = odeint(
        rhs, [1.0, 2.0], time, args=(0.7,), tfirst=True, rtol=1e-9, atol=1e-11, **options
    )
    actual = solve_evolution(
        rhs,
        [1.0, 2.0],
        time,
        args=(0.7,),
        method="odeint",
        rtol=1e-9,
        atol=1e-11,
        odeint_options=options,
    )
    np.testing.assert_array_equal(actual, expected)
    assert options == {"hmax": 0.03, "h0": -0.001 if reverse else 0.001, "mxstep": 5000}


def test_stiff_jacobian_uses_canonical_time_first_signature():
    """Use the supplied Jacobian for a system with widely separated decay time-scales.

    Notes
    -----
    Pass criteria: Decay rates 1 and 1000 with a time-first Jacobian reproduce direct odeint
    exactly.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    time = np.linspace(0.0, 1.0, 21)
    rates = np.array([1.0, 1000.0])

    def rhs(t, y, rates):
        """Provide independent linear decay components spanning stiff time scales.

        Parameters
        ----------
        t : float
            Dimensionless integration time of the analytic test equation.
        y : numpy.ndarray
            Current state of the analytic differential equation, in the fixture convention.
        rates : numpy.ndarray
            Positive analytic decay rates; the large rate separation exercises stiff
            integration.

        Returns
        -------
        numpy.ndarray
            -rates*y for the controlled vector state.

        Notes
        -----
        Use the supplied Jacobian for a system with widely separated decay time-scales. Pass
        criteria: Decay rates 1 and 1000 with a time-first Jacobian reproduce direct odeint
        exactly. Failure occurs when the stated comparison or expected rejection is not met.
        This controlled calculation does not calibrate a subhalo population against simulations
        or observations.
        """
        return -rates * y

    def jac(t, y, rates):
        """Return the exact decay Jacobian in the dense or banded layout requested by the test.

        Parameters
        ----------
        t : float
            Dimensionless integration time of the analytic test equation.
        y : numpy.ndarray
            Current state of the analytic differential equation, in the fixture convention.
        rates : numpy.ndarray
            Positive analytic decay rates; the large rate separation exercises stiff
            integration.

        Returns
        -------
        numpy.ndarray
            Negative diagonal decay rates, with the tested storage orientation.

        Notes
        -----
        Use the supplied Jacobian for a system with widely separated decay time-scales. Pass
        criteria: Decay rates 1 and 1000 with a time-first Jacobian reproduce direct odeint
        exactly. Failure occurs when the stated comparison or expected rejection is not met.
        This controlled calculation does not calibrate a subhalo population against simulations
        or observations.
        """
        return -np.diag(rates)

    expected = odeint(
        rhs, [1.0, 1.0], time, args=(rates,), Dfun=jac, tfirst=True, rtol=1e-9, atol=1e-11
    )
    actual = solve_evolution(
        rhs,
        [1.0, 1.0],
        time,
        args=(rates,),
        method="odeint",
        rtol=1e-9,
        atol=1e-11,
        odeint_options={"Dfun": jac},
    )
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    "options",
    [
        {"hmax": -1.0},
        {"hmin": np.nan},
        {"h0": np.inf},
        {"mxstep": 1.5},
        {"mxstep": -1},
        {"mxordn": 13},
        {"mxords": 6},
        {"ixpr": 2},
        {"tcrit": [np.nan]},
        {"Dfun": 3},
        {"full_output": True},
        {"tfirst": False},
        {"args": ()},
        {"rtol": 1e-3},
        {"unknown": 1},
    ],
)
def test_invalid_controls_fail_before_rhs(options):
    """Reject ambiguous integration settings before invoking the physical model.

    Parameters
    ----------
    options : dict
        Numerical or provider settings for the controlled valid/invalid configuration case.

    Notes
    -----
    Pass criteria: Invalid step bounds, orders, work counts, Jacobian settings and
    unsupported overrides raise ValueError before the RHS runs.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    def rhs(t, y):
        """Fail immediately if invalid solver controls reach the physical derivative.

        Parameters
        ----------
        t : float
            Dimensionless integration time of the analytic test equation.
        y : numpy.ndarray
            Current state of the analytic differential equation, in the fixture convention.

        Notes
        -----
        Reject ambiguous integration settings before invoking the physical model. Pass criteria:
        Invalid step bounds, orders, work counts, Jacobian settings and unsupported overrides
        raise ValueError before the RHS runs. Failure occurs when the stated comparison or
        expected rejection is not met. This controlled calculation does not calibrate a subhalo
        population against simulations or observations.
        """
        pytest.fail("invalid integration controls reached the physical callback")

    with pytest.raises(ValueError):
        solve_evolution(rhs, [1.0], [0.0, 1.0], method="odeint", odeint_options=options)


def test_odeint_controls_are_not_silently_used_by_another_solver():
    """Prevent a user from believing LSODA-specific accuracy controls affected another solver.

    Notes
    -----
    Pass criteria: Passing odeint_options without selecting odeint raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="odeint"):
        solve_evolution(lambda t, y: -y, [1.0], [0.0, 1.0], odeint_options={"hmax": 0.1})


def test_failed_internal_step_limit_is_not_reported_as_success():
    """Do not return an uncompleted physical trajectory as a valid solution.

    Notes
    -----
    Pass criteria: For mxstep=1 on a nontrivial decay interval, the wrapper raises
    RuntimeError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(RuntimeError, match="odeint"):
        solve_evolution(
            lambda t, y: -y, [1.0], [0.0, 10.0], method="odeint", odeint_options={"mxstep": 1}
        )


def test_success_message_remains_a_diagnostic():
    """Distinguish an informational LSODA success warning from failed integration.

    Notes
    -----
    Pass criteria: printmessg emits the expected success warning while the final decay value
    agrees with exp(-1) to rtol=1e-7.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.warns(ODEintWarning, match="Integration successful"):
        result = solve_evolution(
            lambda t, y: -y, [1.0], [0.0, 1.0], method="odeint", odeint_options={"printmessg": True}
        )
    np.testing.assert_allclose(result[-1], np.exp(-1.0), rtol=1e-7)


@pytest.mark.parametrize("column_derivatives", [False, True])
def test_banded_jacobian_shape_and_orientation(column_derivatives):
    """Interpret row- and column-oriented banded Jacobians consistently in stiff evolution.

    Parameters
    ----------
    column_derivatives : bool
        Requested row/column storage orientation for a banded Jacobian.

    Notes
    -----
    Pass criteria: Both diagonal-band orientations reproduce direct odeint for decay rates 1
    and 1000.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    rates = np.array([1.0, 1000.0])
    time = np.linspace(0.0, 1.0, 21)

    def rhs(t, y):
        """Provide independent linear decay components spanning stiff time scales.

        Parameters
        ----------
        t : float
            Dimensionless integration time of the analytic test equation.
        y : numpy.ndarray
            Current state of the analytic differential equation, in the fixture convention.

        Returns
        -------
        numpy.ndarray
            -rates*y for the controlled vector state.

        Notes
        -----
        Interpret row- and column-oriented banded Jacobians consistently in stiff evolution.
        Pass criteria: Both diagonal-band orientations reproduce direct odeint for decay rates 1
        and 1000. Failure occurs when the stated comparison or expected rejection is not met.
        This controlled calculation does not calibrate a subhalo population against simulations
        or observations.
        """
        return -rates * y

    def jac(t, y):
        """Return the exact decay Jacobian in the dense or banded layout requested by the test.

        Parameters
        ----------
        t : float
            Dimensionless integration time of the analytic test equation.
        y : numpy.ndarray
            Current state of the analytic differential equation, in the fixture convention.

        Returns
        -------
        numpy.ndarray
            Negative diagonal decay rates, with the tested storage orientation.

        Notes
        -----
        Interpret row- and column-oriented banded Jacobians consistently in stiff evolution.
        Pass criteria: Both diagonal-band orientations reproduce direct odeint for decay rates 1
        and 1000. Failure occurs when the stated comparison or expected rejection is not met.
        This controlled calculation does not calibrate a subhalo population against simulations
        or observations.
        """
        result = -rates[None, :]
        return result.T if column_derivatives else result

    options = {"Dfun": jac, "ml": 0, "mu": 0, "col_deriv": column_derivatives}
    expected = odeint(rhs, [1.0, 1.0], time, tfirst=True, **options)
    actual = solve_evolution(rhs, [1.0, 1.0], time, method="odeint", odeint_options=options)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    "bad_jacobian", [np.ones((1, 2)), np.full((2, 2), np.nan), np.ones((2, 2), dtype=complex)]
)
def test_invalid_jacobian_is_not_returned_as_a_solution(bad_jacobian):
    """Reject an invalid linearization of the physical evolution equations.

    Parameters
    ----------
    bad_jacobian : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Wrong-shape, NaN or complex Jacobians raise ValueError rather than
    returning a trajectory.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="Dfun"):
        solve_evolution(
            lambda t, y: -np.array([1.0, 1000.0]) * y,
            [1.0, 1.0],
            [0.0, 1.0],
            method="odeint",
            odeint_options={"Dfun": lambda t, y: bad_jacobian},
        )


@pytest.mark.parametrize("time", [[1.0, 1.0, 0.5, 0.5, 0.0], [1.0, 1.0, 1.0], [1.0]])
def test_explicit_repeated_outputs_preserve_scipy_shape_and_values(time):
    """Allow repeated observation epochs without dropping or moving physical samples.

    Parameters
    ----------
    time : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Descending, constant and singleton time grids match SciPy shapes and
    values when repeats are explicitly enabled.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    expected = odeint(lambda y, t: -0.2 * y, [1.0, 3.0], time)
    actual = solve_evolution(
        lambda t, y: -0.2 * y, [1.0, 3.0], time, method="odeint", allow_repeated_times=True
    )
    np.testing.assert_array_equal(actual, expected)


def test_constant_grid_still_validates_solver_options():
    """Apply the same numerical-control validation to zero-duration evolution.

    Notes
    -----
    Pass criteria: A negative hmax is rejected even when every requested epoch is identical.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="hmax"):
        solve_evolution(
            lambda t, y: -y,
            [1.0],
            [0.0, 0.0],
            method="odeint",
            allow_repeated_times=True,
            odeint_options={"hmax": -1.0},
        )
