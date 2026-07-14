import numpy as np
import pytest

from itamae.evolution import PerturbativeEvolutionSolver, shanks_transform, solve_evolution


def test_ode_exponential_decay():
    t = np.linspace(0.0, 2.0, 21)
    y = solve_evolution(lambda time, state, rate: -rate * state, [1.0], t, args=(2.0,))
    assert y[:, 0] == pytest.approx(np.exp(-2.0 * t), rel=1e-6)


def test_shanks_geometric_sequence():
    assert shanks_transform(1.0, 1.5, 1.75) == pytest.approx(2.0)


def test_perturbative_solver_constant_rhs_and_exponential_sequence():
    time = np.linspace(0.0, 1.0, 401)
    constant = PerturbativeEvolutionSolver(order=2).solve(
        lambda grid, state: -np.ones_like(state), 2.0, time
    )
    np.testing.assert_allclose(constant, 2.0 - time)

    exponential = PerturbativeEvolutionSolver(order=8).solve(lambda grid, state: -state, 1.0, time)
    np.testing.assert_allclose(exponential, np.exp(-time), rtol=1.0e-6, atol=1.0e-7)


def test_perturbative_solver_validates_grid_and_rhs():
    solver = PerturbativeEvolutionSolver()
    with pytest.raises(ValueError, match="strictly monotonic"):
        solver.solve(lambda grid, state: state, 1.0, [0.0, 1.0, 0.5])
    with pytest.raises(ValueError, match="non-finite"):
        solver.solve(lambda grid, state: np.full_like(state, np.nan), 1.0, [0.0, 1.0])
