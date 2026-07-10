import numpy as np
import pytest

from itamae.evolution import shanks_transform, solve_evolution


def test_ode_exponential_decay():
    t = np.linspace(0.0, 2.0, 21)
    y = solve_evolution(lambda time, state, rate: -rate * state, [1.0], t, args=(2.0,))
    assert y[:, 0] == pytest.approx(np.exp(-2.0 * t), rel=1e-6)


def test_shanks_geometric_sequence():
    assert shanks_transform(1.0, 1.5, 1.75) == pytest.approx(2.0)
