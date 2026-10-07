"""NFW boundaries checked against high-precision standard-library arithmetic."""

from decimal import Decimal, localcontext
import warnings

import numpy as np
import pytest

from itamae.halo import NFWProfile, invert_nfw_mass_function, nfw_mass_function


def decimal_mass(x):
    with localcontext() as context:
        context.prec = 340
        value = Decimal.from_float(float(x))
        return float((1 + value).ln() - value / (1 + value))


def test_enclosed_mass_retains_small_radius_precision():
    radius = np.array([1e-150, 1e-12, 1e-6, 1e-3, 0.1, 1.0, 10.0, 1e6])
    expected = np.array([decimal_mass(x) for x in radius])
    np.testing.assert_allclose(nfw_mass_function(radius), expected, rtol=5e-14, atol=0.0)


def test_inverse_preserves_small_positive_radii():
    radius = np.geomspace(1e-150, 1e6, 20)
    independent_mass = [decimal_mass(x) for x in radius]
    recovered = invert_nfw_mass_function(independent_mass)
    np.testing.assert_allclose(recovered, radius, rtol=5e-12, atol=0.0)
    assert invert_nfw_mass_function(0.0) == 0.0
    assert invert_nfw_mass_function(np.array([])).shape == (0,)


@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf, -1.0])
def test_nfw_mass_rejects_nonfinite_or_negative_inputs(invalid):
    with pytest.raises(ValueError, match="finite|nonnegative"):
        nfw_mass_function(invalid)
    with pytest.raises(ValueError, match="finite|nonnegative"):
        invert_nfw_mass_function(invalid)


@pytest.mark.parametrize("field", ["r_s", "rho_s"])
@pytest.mark.parametrize("invalid", [0.0, -1.0, np.nan, np.inf])
def test_profile_requires_finite_positive_parameters(field, invalid):
    with pytest.raises(ValueError, match="finite.*positive"):
        NFWProfile(**{**{"r_s": 0.02, "rho_s": 1.0e15}, field: invalid})


def test_central_potential_uses_its_limit_without_evaluating_zero_over_zero():
    model = NFWProfile(r_s=0.02, rho_s=1.0e15)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        expected = -4 * np.pi * 4.30091e-9 * model.rho_s * model.r_s**2
        assert model.potential(0.0) == pytest.approx(expected, rel=2e-15)
    with pytest.raises(ValueError, match="finite|nonnegative"):
        model.potential(-1.0)
