"""Independent solver and failure-boundary checks."""

import numpy as np
import pytest
from scipy.integrate import odeint

from itamae.evolution import solve_evolution


@pytest.mark.parametrize("direction", [1.0, -1.0])
def test_odeint_preserves_scipy_defaults_and_time_order(direction):
    grid = direction * np.linspace(0.0, 1.0, 19)
    initial = np.array([2.0, 7.0])
    rate = np.array([0.3, 1.2])
    expected = odeint(lambda state, time: -rate * state, initial, grid)
    actual = solve_evolution(
        lambda time, state, rates: -rates * state,
        initial,
        grid,
        args=(rate,),
        method="odeint",
    )
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_allclose(actual, initial * np.exp(-grid[:, None] * rate), rtol=1e-7)


@pytest.mark.parametrize("method", ["RK45", "odeint"])
@pytest.mark.parametrize("bad_rhs", [lambda t, y: [np.nan], lambda t, y: [1.0, 2.0]])
def test_invalid_rhs_is_rejected(method, bad_rhs):
    with pytest.raises(ValueError, match="rhs"):
        solve_evolution(bad_rhs, [1.0], [0.0, 1.0], method=method)


@pytest.mark.parametrize("method", ["RK45", "odeint"])
@pytest.mark.parametrize(
    "initial,grid,tolerance",
    [
        ([np.inf], [0, 1], None),
        ([], [0, 1], None),
        ([1], [0, np.inf], None),
        ([1], [0, 1], -1.0),
        ([1], [0, 1], np.nan),
        ([1j], [0, 1], None),
    ],
)
def test_invalid_inputs_are_rejected(method, initial, grid, tolerance):
    with pytest.raises(ValueError):
        solve_evolution(lambda t, y: -y, initial, grid, method=method, rtol=tolerance)
