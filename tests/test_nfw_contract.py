"""NFW boundaries checked against high-precision standard-library arithmetic."""

from decimal import Decimal, localcontext
import warnings

import numpy as np
import pytest

from itamae.halo import NFWProfile, invert_nfw_mass_function, nfw_mass_function


def decimal_mass(x):
    """Evaluate the NFW enclosed-mass function with high-precision Decimal arithmetic.

    Parameters
    ----------
    x : float
        Nonnegative NFW radius ratio r/r_s, including the cancellation-sensitive
        small-radius regime.

    Returns
    -------
    float
        Independent f(x) reference rounded only after 340-digit evaluation.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    with localcontext() as context:
        context.prec = 340
        value = Decimal.from_float(float(x))
        return float((1 + value).ln() - value / (1 + value))


def test_enclosed_mass_retains_small_radius_precision():
    """Resolve the NFW central mass without cancellation in log(1+x)-x/(1+x).

    Notes
    -----
    Pass criteria: The dimensionless mass agrees with 340-digit Decimal evaluation at radii
    1e-150..1e6 with rtol=5e-14 and zero absolute tolerance.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    radius = np.array([1e-150, 1e-12, 1e-6, 1e-3, 0.1, 1.0, 10.0, 1e6])
    expected = np.array([decimal_mass(x) for x in radius])
    np.testing.assert_allclose(nfw_mass_function(radius), expected, rtol=5e-14, atol=0.0)


def test_inverse_preserves_small_positive_radii():
    """Avoid turning a small positive bound mass into a zero truncation radius.

    Notes
    -----
    Pass criteria: Inverting independent Decimal masses over radii 1e-150..1e6 recovers
    radii to rtol=5e-12; zero and empty inputs keep their meanings.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    radius = np.geomspace(1e-150, 1e6, 20)
    independent_mass = [decimal_mass(x) for x in radius]
    recovered = invert_nfw_mass_function(independent_mass)
    np.testing.assert_allclose(recovered, radius, rtol=5e-12, atol=0.0)
    assert invert_nfw_mass_function(0.0) == 0.0
    assert invert_nfw_mass_function(np.array([])).shape == (0,)


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf, -1.0])
def test_nfw_mass_rejects_nonfinite_or_negative_inputs(invalid):
    """Reject a radius or enclosed-mass target that has no physical NFW interpretation.

    Parameters
    ----------
    invalid : object
        Parametrized input outside the supported physical or numerical domain.

    Notes
    -----
    Pass criteria: NaN, infinities and negative values raise ValueError in both mass and
    inverse functions.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="finite|nonnegative"):
        nfw_mass_function(invalid)
    with pytest.raises(ValueError, match="finite|nonnegative"):
        invert_nfw_mass_function(invalid)


@pytest.mark.parametrize("field", ["r_s", "rho_s"])
@pytest.mark.parametrize("invalid", [0.0, -1.0, np.nan, np.inf])
def test_profile_requires_finite_positive_parameters(field, invalid):
    """Require a physical NFW scale radius and scale density.

    Parameters
    ----------
    field : str
        Physical field or provenance field deliberately varied in this test.
    invalid : object
        Parametrized input outside the supported physical or numerical domain.

    Notes
    -----
    Pass criteria: Zero, negative or nonfinite r_s or rho_s raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="finite.*positive"):
        NFWProfile(**{**{"r_s": 0.02, "rho_s": 1.0e15}, field: invalid})


def test_central_potential_uses_its_limit_without_evaluating_zero_over_zero():
    """Evaluate the finite NFW central potential by its analytic limit.

    Notes
    -----
    Pass criteria: At r=0 the potential equals -4*pi*G*rho_s*r_s**2 to rtol=2e-15 without
    RuntimeWarning; negative radius fails.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    model = NFWProfile(r_s=0.02, rho_s=1.0e15)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        expected = -4 * np.pi * 4.30091e-9 * model.rho_s * model.r_s**2
        assert model.potential(0.0) == pytest.approx(expected, rel=2e-15)
    with pytest.raises(ValueError, match="finite|nonnegative"):
        model.potential(-1.0)
