"""Independent integral and derivative limits for a cutoff spectrum."""

from dataclasses import dataclass, replace

import numpy as np
from scipy.special import gammainc, gamma

from itamae.power import SharpKWindow
from itamae.variance import IntegratedVarianceModel


@dataclass(frozen=True)
class ExponentialPower:
    identifier: str = "test:exp(-k)"

    def __call__(self, k):
        return np.exp(-np.asarray(k))


def model(n_k=257, chunk_size=256):
    return IntegratedVarianceModel(
        ExponentialPower(),
        SharpKWindow(),
        rho_mean=3 / (4 * np.pi),
        k_min=1e-6,
        k_max=1e4,
        n_k=n_k,
        chunk_size=chunk_size,
    )


def test_cutoff_plateau_is_monotonic_without_projection():
    m = np.geomspace(1e-15, 1e12, 2000)
    calculation = model()
    s = calculation.variance(m)
    assert np.all(np.diff(s) <= 0.0)
    cutoff = np.minimum(m ** (-1 / 3), calculation.k_max)
    expected = (gammainc(3, cutoff) - gammainc(3, calculation.k_min)) * gamma(3) / (2 * np.pi**2)
    np.testing.assert_allclose(s, expected, rtol=3e-12, atol=3e-16)


def test_moving_boundary_derivative_agrees_with_independent_integral():
    m = np.geomspace(1e-3, 1e3, 30)
    calculation = model()
    step = 1e-4
    numerical = (
        calculation.variance(m * np.exp(step)) - calculation.variance(m * np.exp(-step))
    ) / (m * (np.exp(step) - np.exp(-step)))
    np.testing.assert_allclose(calculation.dvariance_dmass(m), numerical, rtol=4e-8, atol=1e-15)
    np.testing.assert_array_equal(
        calculation.variance(m), replace(calculation, chunk_size=3).variance(m)
    )
    np.testing.assert_allclose(calculation.variance(m), model(513).variance(m), rtol=2e-13)
