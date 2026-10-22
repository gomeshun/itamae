"""Independent integral and derivative limits for a cutoff spectrum."""

from dataclasses import dataclass, replace

import numpy as np
from scipy.special import gammainc, gamma

from itamae.power import SharpKWindow, TabulatedPowerSpectrum, TransferModifiedPowerSpectrum
from itamae.variance import IntegratedVarianceModel


@dataclass(frozen=True)
class ExponentialPower:
    """Supply P(k)=exp(-k) with an analytic sharp-k integral.

    Parameters
    ----------
    identifier : str
        Stable identifier of the physical and numerical prescription, used to distinguish
        calculations and cached data.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """

    identifier: str = "test:exp(-k)"

    def __call__(self, k):
        """Evaluate the positive analytic spectrum used to separate integration error from
        interpolation.

        Parameters
        ----------
        k : float or array_like
            Positive Fourier wavenumber, in the inverse-length convention of the supplied
            spectrum; factors of h are owned by that spectrum.

        Returns
        -------
        numpy.ndarray
            exp(-k) with the input wavenumber shape.

        Notes
        -----
        This fixture/helper supports the controlled checks in this test module; its synthetic
        values are not a calibrated halo model.
        """
        return np.exp(-np.asarray(k))


def model(n_k=257, chunk_size=256):
    """Choose a sharp-k variance fixture with R=M**(1/3) and an analytic cutoff.

    Parameters
    ----------
    n_k : int
        Number of logarithmic integration-grid nodes, at least three. Accuracy requires
        refinement checks for the chosen spectrum.
    chunk_size : int
        Maximum mass or quadrature-cell count per allocation; controls memory use without
        changing the intended integral.

    Returns
    -------
    IntegratedVarianceModel
        Analytic-spectrum variance provider with explicit log-k resolution and memory
        chunking.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
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
    """Resolve a suppressed spectrum's variance plateau without imposing monotonicity
    afterward.

    Notes
    -----
    Pass criteria: For P(k)=exp(-k), variance decreases with mass and matches the
    incomplete-gamma integral at rtol=3e-12, atol=3e-16.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    m = np.geomspace(1e-15, 1e12, 2000)
    calculation = model()
    s = calculation.variance(m)
    assert np.all(np.diff(s) <= 0.0)
    cutoff = np.minimum(m ** (-1 / 3), calculation.k_max)
    expected = (gammainc(3, cutoff) - gammainc(3, calculation.k_min)) * gamma(3) / (2 * np.pi**2)
    np.testing.assert_allclose(s, expected, rtol=3e-12, atol=3e-16)


def test_moving_boundary_derivative_agrees_with_independent_integral():
    """Make the sharp-k derivative describe the same continuous variance integral.

    Notes
    -----
    Pass criteria: A finite difference agrees with the moving-boundary derivative to
    rtol=4e-8; chunking preserves values exactly; doubling k-grid resolution agrees to
    2e-13.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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


def test_tabulated_spectrum_knots_have_exact_piecewise_polynomial_integrals():
    """Integrate across input-spectrum knots without missing changes of interpolation slope.

    Notes
    -----
    Pass criteria: Variance at six off-knot cutoffs matches direct piecewise polynomial
    integration to rtol=3e-13.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    k = np.array([0.01, 0.1234, 0.568, 3.141, 100.0])
    p = np.array([1.0, 8.0, 0.3, 2.0, 0.001])
    base = TabulatedPowerSpectrum(k, p, interpolation="linear")
    power = TransferModifiedPowerSpectrum(base, lambda x: np.ones_like(x), ratio_identifier="unity")
    calculation = IntegratedVarianceModel(
        power,
        SharpKWindow(),
        rho_mean=3 / (4 * np.pi),
        k_min=k[0],
        k_max=k[-1],
        n_k=101,
    )
    cutoffs = np.array([0.113, 0.234, 0.777, 2.555, 7.777, 99.0])
    expected = np.zeros(cutoffs.size)
    for i, cutoff in enumerate(cutoffs):
        for left, right, p0, p1 in zip(k[:-1], k[1:], p[:-1], p[1:]):
            if cutoff <= left:
                break
            upper = min(right, cutoff)
            slope = (p1 - p0) / (right - left)
            intercept = p0 - slope * left
            expected[i] += (
                slope * (upper**4 - left**4) / 4 + intercept * (upper**3 - left**3) / 3
            ) / (2 * np.pi**2)
    np.testing.assert_allclose(calculation.variance(cutoffs**-3), expected, rtol=3e-13)
