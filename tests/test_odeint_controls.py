"""Preserve SciPy integration controls while rejecting ambiguous overrides."""

import numpy as np
import pytest
from scipy.integrate import ODEintWarning, odeint
from itamae.evolution import solve_evolution


@pytest.mark.parametrize("reverse", [False, True])
def test_selected_controls_preserve_direct_odeint_values(reverse):
    time = np.linspace(0.0, 2.0, 23)
    if reverse:
        time = time[::-1]

    def rhs(t, y, rate):
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
    time = np.linspace(0.0, 1.0, 21)
    rates = np.array([1.0, 1000.0])

    def rhs(t, y, rates):
        return -rates * y

    def jac(t, y, rates):
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
    def rhs(t, y):
        pytest.fail("invalid integration controls reached the physical callback")

    with pytest.raises(ValueError):
        solve_evolution(rhs, [1.0], [0.0, 1.0], method="odeint", odeint_options=options)


def test_odeint_controls_are_not_silently_used_by_another_solver():
    with pytest.raises(ValueError, match="odeint"):
        solve_evolution(lambda t, y: -y, [1.0], [0.0, 1.0], odeint_options={"hmax": 0.1})


def test_failed_internal_step_limit_is_not_reported_as_success():
    with pytest.raises(RuntimeError, match="odeint"):
        solve_evolution(
            lambda t, y: -y, [1.0], [0.0, 10.0], method="odeint", odeint_options={"mxstep": 1}
        )


def test_success_message_remains_a_diagnostic():
    with pytest.warns(ODEintWarning, match="Integration successful"):
        result = solve_evolution(
            lambda t, y: -y, [1.0], [0.0, 1.0], method="odeint", odeint_options={"printmessg": True}
        )
    np.testing.assert_allclose(result[-1], np.exp(-1.0), rtol=1e-7)


@pytest.mark.parametrize("column_derivatives", [False, True])
def test_banded_jacobian_shape_and_orientation(column_derivatives):
    rates = np.array([1.0, 1000.0])
    time = np.linspace(0.0, 1.0, 21)

    def rhs(t, y):
        return -rates * y

    def jac(t, y):
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
    with pytest.raises(ValueError, match="Dfun"):
        solve_evolution(
            lambda t, y: -np.array([1.0, 1000.0]) * y,
            [1.0, 1.0],
            [0.0, 1.0],
            method="odeint",
            odeint_options={"Dfun": lambda t, y: bad_jacobian},
        )
