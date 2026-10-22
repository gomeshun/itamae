import numpy as np
import pytest

from itamae.evolution import PerturbativeEvolutionSolver, shanks_transform, solve_evolution


def test_ode_exponential_decay():
    """Recover a known loss history with the generic evolution integrator.

    Notes
    -----
    Pass criteria: For dy/dt=-2y and y(0)=1, 21 outputs over t=0..2 agree with exp(-2t) to
    relative tolerance 1e-6.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    t = np.linspace(0.0, 2.0, 21)
    y = solve_evolution(lambda time, state, rate: -rate * state, [1.0], t, args=(2.0,))
    assert y[:, 0] == pytest.approx(np.exp(-2.0 * t), rel=1e-6)


def test_shanks_geometric_sequence():
    """Check the sequence-acceleration formula used by the historical stripping approximation.

    Notes
    -----
    Pass criteria: Applying Shanks to 1,1.5,1.75 recovers the geometric-series limit 2.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    assert shanks_transform(1.0, 1.5, 1.75) == pytest.approx(2.0)


def test_perturbative_solver_constant_rhs_and_exponential_sequence():
    """Recover analytic constant-loss and exponential histories through Picard iteration.

    Notes
    -----
    Pass criteria: The constant RHS yields 2-t; order-eight iteration follows exp(-t) with
    rtol=1e-6 and atol=1e-7 on 401 nodes.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    time = np.linspace(0.0, 1.0, 401)
    constant = PerturbativeEvolutionSolver(order=2).solve(
        lambda grid, state: -np.ones_like(state), 2.0, time
    )
    np.testing.assert_allclose(constant, 2.0 - time)

    exponential = PerturbativeEvolutionSolver(order=8).solve(lambda grid, state: -state, 1.0, time)
    np.testing.assert_allclose(exponential, np.exp(-time), rtol=1.0e-6, atol=1.0e-7)


def test_perturbative_solver_validates_grid_and_rhs():
    """Stop a numerical evolution whose time order or derivative is undefined.

    Notes
    -----
    Pass criteria: A nonmonotonic time grid or NaN right-hand side raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    solver = PerturbativeEvolutionSolver()
    with pytest.raises(ValueError, match="strictly monotonic"):
        solver.solve(lambda grid, state: state, 1.0, [0.0, 1.0, 0.5])
    with pytest.raises(ValueError, match="non-finite"):
        solver.solve(lambda grid, state: np.full_like(state, np.nan), 1.0, [0.0, 1.0])
