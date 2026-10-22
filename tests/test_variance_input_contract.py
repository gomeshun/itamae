"""Public variance boundaries must not conceal invalid values or broadcasting."""

from dataclasses import replace

import numpy as np
import pytest

from itamae.power import SharpKWindow, SphericalTopHatWindow
from itamae.variance import CallableVarianceModel, IntegratedVarianceModel


class ConstantPower:
    """Supply unit power to isolate variance input and shape checks.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """

    identifier = "constant:test"

    def __call__(self, k):
        """Evaluate unit power at every requested Fourier coordinate.

        Parameters
        ----------
        k : float or array_like
            Positive Fourier wavenumber, in the inverse-length convention of the supplied
            spectrum; factors of h are owned by that spectrum.

        Returns
        -------
        numpy.ndarray
            Ones with exactly the wavenumber shape.

        Notes
        -----
        This fixture/helper supports the controlled checks in this test module; its synthetic
        values are not a calibrated halo model.
        """
        return np.ones_like(k)


def integrated(**options):
    """Create a finite sharp-k integral with controlled provider and coordinate bounds.

    Parameters
    ----------
    **options : Mapping, optional
        Explicit numerical controls for the named solver or table. Unsupported or reserved
        keys are rejected by the receiving routine.

    Returns
    -------
    IntegratedVarianceModel
        Default unit-spectrum provider with test-specific configuration overrides.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    return IntegratedVarianceModel(
        **dict(
            power=ConstantPower(),
            window=SharpKWindow(),
            rho_mean=1.0,
            k_min=0.01,
            k_max=10.0,
            n_k=17,
        )
        | options
    )


def wrapped(**options):
    """Create controlled sigma and dS/dM callbacks for public-input validation.

    Parameters
    ----------
    **options : Mapping, optional
        Explicit numerical controls for the named solver or table. Unsupported or reserved
        keys are rejected by the receiving routine.

    Returns
    -------
    CallableVarianceModel
        Shape-consistent toy callbacks; these are contract fixtures, not a mutually
        calibrated physical variance law.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    return CallableVarianceModel(
        **dict(
            identifier="callable:test",
            sigma_function=lambda m, z: np.ones(np.broadcast_shapes(np.shape(m), np.shape(z))),
            derivative_function=lambda m, z: (
                -np.ones(np.broadcast_shapes(np.shape(m), np.shape(z)))
            ),
        )
        | options
    )


@pytest.mark.parametrize("factory", [wrapped, integrated])
@pytest.mark.parametrize("method", ["sigma", "variance", "dvariance_dmass"])
@pytest.mark.parametrize(
    "mass,z", [(1 + 1j, 0), (1, 1 + 1j), (0, 0), (-1, 0), (np.inf, 0), (1, np.nan)]
)
def test_invalid_coordinates_fail(factory, method, mass, z):
    """Reject nonphysical mass and undefined redshift at every variance entry point.

    Parameters
    ----------
    factory : callable
        Alternative public model/adapter constructor that must obey the same stated
        boundary.
    method : str
        Public evaluation or numerical solver path whose behavior is checked.
    mass : object
        Halo mass or deliberately invalid mass input for this controlled domain case.
    z : object
        Dimensionless redshift or deliberately invalid epoch.

    Notes
    -----
    Pass criteria: Callable and integrated sigma, variance and derivative reject complex
    coordinates, nonpositive mass, infinite mass and NaN redshift.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError):
        getattr(factory(), method)(mass, z)


@pytest.mark.parametrize(
    "options", [{"identifier": " "}, {"sigma_function": None}, {"derivative_function": 2}]
)
def test_callable_configuration_is_validated(options):
    """Require identifiable and callable prescriptions for model-supplied variance.

    Parameters
    ----------
    options : dict
        Numerical or provider settings for the controlled valid/invalid configuration case.

    Notes
    -----
    Pass criteria: Blank identity or noncallable sigma/derivative providers raises an error.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises((ValueError, TypeError)):
        wrapped(**options)


@pytest.mark.parametrize("result", [np.nan, np.inf, 1 + 2j, [1, 2, 3]])
@pytest.mark.parametrize(
    "method,keyword", [("sigma", "sigma_function"), ("dvariance_dmass", "derivative_function")]
)
def test_callback_results_cannot_corrupt_or_expand_output(result, method, keyword):
    """Keep returned variance values finite and aligned with the requested halo masses.

    Parameters
    ----------
    result : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.
    method : str
        Public evaluation or numerical solver path whose behavior is checked.
    keyword : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: NaN, infinity, complex values and an incompatible result length raise
    ValueError for sigma and derivative callbacks.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError):
        getattr(wrapped(**{keyword: lambda m, z: result}), method)([1, 2])


def test_negative_sigma_and_square_overflow_fail():
    """Keep the rms density contrast nonnegative and its variance finite.

    Notes
    -----
    Pass criteria: A negative sigma or sigma whose square overflows raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError):
        wrapped(sigma_function=lambda m, z: -1.0).sigma(1.0)
    with pytest.raises(ValueError):
        wrapped(sigma_function=lambda m, z: 1.0e200).variance(1.0)


def test_scalar_callback_can_broadcast_but_invalid_input_never_reaches_callback():
    """Allow a constant physical prescription while rejecting invalid masses early.

    Notes
    -----
    Pass criteria: A scalar sigma broadcasts over a 2-by-3 mass/redshift grid; negative mass
    fails before calling the supplied model.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    def forbidden(m, z):
        """Fail if a malformed mass or redshift reaches a variance callback.

        Parameters
        ----------
        m : array_like
            Deliberately invalid mass coordinate that should be rejected before callback
            evaluation.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Notes
        -----
        Allow a constant physical prescription while rejecting invalid masses early. Pass
        criteria: A scalar sigma broadcasts over a 2-by-3 mass/redshift grid; negative mass
        fails before calling the supplied model. Failure occurs when the stated comparison or
        expected rejection is not met. This controlled calculation does not calibrate a subhalo
        population against simulations or observations.
        """
        pytest.fail("invalid input reached callback")

    with pytest.raises(ValueError):
        wrapped(sigma_function=forbidden).sigma(-1.0)
    value = wrapped(sigma_function=lambda m, z: 2.0).sigma([[1.0], [2.0]], [0.0, 1.0, 2.0])
    np.testing.assert_array_equal(value, np.full((2, 3), 2.0))


@pytest.mark.parametrize("window", [SharpKWindow(), SphericalTopHatWindow()])
def test_complex_power_is_rejected_before_cast(window):
    """Prevent complex power from being silently used in real density-variance quadrature.

    Parameters
    ----------
    window : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Sharp-k and top-hat paths reject complex spectra for both variance and
    its derivative.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    class ComplexPower(ConstantPower):
        """Supply deliberately complex power to detect silent removal of an imaginary component.

        Notes
        -----
        Prevent complex power from being silently used in real density-variance quadrature. Pass
        criteria: Sharp-k and top-hat paths reject complex spectra for both variance and its
        derivative. Failure occurs when the stated comparison or expected rejection is not met.
        This controlled calculation does not calibrate a subhalo population against simulations
        or observations.
        """

        def __call__(self, k):
            """Return complex power that the real density-variance calculation must reject.

            Parameters
            ----------
            k : float or array_like
                Positive Fourier wavenumber, in the inverse-length convention of the supplied
                spectrum; factors of h are owned by that spectrum.

            Returns
            -------
            numpy.ndarray
                Complex 1+i values with the query shape.

            Notes
            -----
            Prevent complex power from being silently used in real density-variance quadrature. Pass
            criteria: Sharp-k and top-hat paths reject complex spectra for both variance and its
            derivative. Failure occurs when the stated comparison or expected rejection is not met.
            This controlled calculation does not calibrate a subhalo population against simulations
            or observations.
            """
            return np.ones_like(k) * (1 + 1j)

    with pytest.raises(ValueError):
        integrated(power=ComplexPower(), window=window).variance(1.0)
    with pytest.raises(ValueError):
        integrated(power=ComplexPower(), window=window).dvariance_dmass(1.0)


@pytest.mark.parametrize("growth", [lambda z: 1 + 1j, lambda z: np.ones((2, 2)), lambda z: 1.0e200])
def test_invalid_growth_and_overflow_fail(growth):
    """Reject growth factors that cannot define real finite variance at each epoch.

    Parameters
    ----------
    growth : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Complex, wrong-shape or overflowing growth factors raise ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    model = integrated(growth_function=growth, growth_identifier="bad:test")
    with pytest.raises(ValueError):
        model.variance([1.0, 2.0])


def test_invalid_window_shape_is_not_silently_broadcast():
    """Keep smoothing-window values on the intended mass-by-wavenumber grid.

    Notes
    -----
    Pass criteria: A window returning an incompatible shape raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    class BadWindow:
        """Supply a wrong-shaped smoothing window to detect accidental extra population axes.

        Notes
        -----
        Keep smoothing-window values on the intended mass-by-wavenumber grid. Pass criteria: A
        window returning an incompatible shape raises ValueError. Failure occurs when the stated
        comparison or expected rejection is not met. This controlled calculation does not
        calibrate a subhalo population against simulations or observations.
        """

        identifier = "bad-window"

        def __call__(self, x):
            """Return a fixed invalid shape instead of one window value per query coordinate.

            Parameters
            ----------
            x : float or array_like
                Nonnegative dimensionless Fourier-window argument x=kR/c; output retains its shape.

            Returns
            -------
            numpy.ndarray
                Deliberately mismatched array of shape (3,1).

            Notes
            -----
            Keep smoothing-window values on the intended mass-by-wavenumber grid. Pass criteria: A
            window returning an incompatible shape raises ValueError. Failure occurs when the stated
            comparison or expected rejection is not met. This controlled calculation does not
            calibrate a subhalo population against simulations or observations.
            """
            return np.ones((3, 1))

    with pytest.raises(ValueError):
        integrated(window=BadWindow()).variance([1.0, 2.0, 3.0])


def test_complex_integration_knots_and_configuration_fail():
    """Require real quadrature breakpoints and mean density for mass smoothing.

    Notes
    -----
    Pass criteria: Complex integration knots or complex rho_mean raise ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    power = ConstantPower()
    power.integration_breakpoints = np.array([0.1 + 1j, 1 + 1j])
    with pytest.raises(ValueError):
        _ = integrated(power=power).identifier
    with pytest.raises(ValueError):
        replace(integrated(), rho_mean=np.complex128(1 + 2j))


def test_valid_broadcast_empty_and_derivative_signs_are_preserved():
    """Preserve legitimate variance broadcasting and user-supplied derivative signs.

    Notes
    -----
    Pass criteria: Analytic sigma and variance broadcast correctly, empty queries stay
    empty, and a positive supplied derivative is retained.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    model = wrapped(
        sigma_function=lambda m, z: np.asarray(m) ** -0.25 / (1 + np.asarray(z)),
        derivative_function=lambda m, z: (
            -(np.asarray(m) ** -0.5) / (2 * np.asarray(m) * (1 + np.asarray(z)) ** 2)
        ),
    )
    mass = np.array([[1.0], [16.0]])
    z = np.array([0.0, 1.0, 2.0])
    np.testing.assert_array_equal(model.sigma(mass, z), mass**-0.25 / (1 + z))
    np.testing.assert_array_equal(model.variance(mass, z), (mass**-0.25 / (1 + z)) ** 2)
    assert model.sigma([]).shape == (0,)
    assert integrated().variance([]).shape == (0,)
    assert wrapped(derivative_function=lambda m, z: 2.0).dvariance_dmass(1.0) == 2.0
