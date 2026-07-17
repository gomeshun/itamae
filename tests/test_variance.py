import numpy as np
import pytest

from itamae.protocols import VarianceModel
from itamae.power import SharpKWindow, TabulatedPowerSpectrum
from itamae.variance import (
    CallableVarianceModel,
    IntegratedVarianceModel,
    load_variance_cache,
    save_variance_cache,
    variance_cache_key,
)


def test_callable_variance_model_satisfies_protocol():
    model = CallableVarianceModel(
        identifier="power-law-test",
        sigma_function=lambda mass, z: np.asarray(mass) ** -0.25 / (1.0 + np.asarray(z)),
        derivative_function=lambda mass, z: (
            -0.5 * np.asarray(mass) ** -1.5 / (1.0 + np.asarray(z)) ** 2
        ),
    )

    assert isinstance(model, VarianceModel)
    mass = np.array([1.0, 16.0, 81.0])
    redshift = 1.0
    sigma = model.sigma(mass, redshift)
    np.testing.assert_allclose(model.variance(mass, redshift), sigma**2)
    np.testing.assert_allclose(
        model.dvariance_dmass(mass, redshift),
        -0.5 * mass**-1.5 / 4.0,
    )
    assert model.identifier == "power-law-test"


def test_integrated_sharp_k_variance_matches_analytic_power_law():
    k = np.logspace(-3, 3, 6001)
    power = TabulatedPowerSpectrum(k, np.ones_like(k), identifier="constant-power")
    model = IntegratedVarianceModel(
        power=power,
        window=SharpKWindow(),
        rho_mean=3.0 / (4.0 * np.pi),
        k_min=1.0e-3,
        k_max=1.0e3,
        n_k=6001,
        filter_scale=1.0,
        growth_function=lambda z: 1.0 / (1.0 + np.asarray(z)),
        growth_identifier="test-growth",
    )

    expected = (1.0 - 1.0e-9) / (6.0 * np.pi**2)
    expected_derivative = -1.0 / (6.0 * np.pi**2)
    assert model.variance(1.0, 0.0) == pytest.approx(expected, rel=2.0e-11)
    assert model.variance(1.0, 1.0) == pytest.approx(expected / 4.0, rel=2.0e-11)
    assert model.dvariance_dmass(1.0, 0.0) == pytest.approx(
        expected_derivative,
        rel=2.0e-14,
    )
    assert model.dvariance_dmass(1.0, 1.0) == pytest.approx(
        expected_derivative / 4.0,
        rel=2.0e-14,
    )
    assert isinstance(model, VarianceModel)
    assert "integrated-variance:v2" in model.identifier
    assert "constant-power" in model.identifier


def test_sharp_k_variance_and_derivative_are_continuous_between_grid_nodes():
    """A moving cutoff must not create zero derivatives or grid-boundary spikes."""
    k = np.geomspace(1.0e-6, 1.0e6, 301)
    power = TabulatedPowerSpectrum(k, k**-2, identifier="inverse-square-power")
    model = IntegratedVarianceModel(
        power=power,
        window=SharpKWindow(),
        rho_mean=3.0 / (4.0 * np.pi),
        k_min=k[0],
        k_max=k[-1],
        n_k=1001,
    )
    mass = np.geomspace(1.0e-12, 1.0e12, 49)
    cutoff = mass ** (-1.0 / 3.0)
    expected_variance = (cutoff - k[0]) / (2.0 * np.pi**2)
    expected_derivative = -cutoff / (6.0 * np.pi**2 * mass)

    np.testing.assert_allclose(
        model.variance(mass),
        expected_variance,
        rtol=2.0e-9,
        atol=1.0e-15,
    )
    np.testing.assert_allclose(
        model.dvariance_dmass(mass),
        expected_derivative,
        rtol=2.0e-14,
        atol=0.0,
    )
    assert np.all(model.dvariance_dmass(mass) < 0.0)


def test_sharp_k_truncated_domain_has_zero_boundary_derivative():
    """The configured finite power domain should define constant end regimes."""
    k = np.geomspace(1.0e-2, 1.0e2, 101)
    power = TabulatedPowerSpectrum(k, np.ones_like(k))
    model = IntegratedVarianceModel(
        power=power,
        window=SharpKWindow(),
        rho_mean=3.0 / (4.0 * np.pi),
        k_min=k[0],
        k_max=k[-1],
        n_k=501,
    )

    assert model.variance(1.0e-9) > 0.0
    assert model.dvariance_dmass(1.0e-9) == 0.0
    assert model.variance(1.0e9) == 0.0
    assert model.dvariance_dmass(1.0e9) == 0.0


def test_integrator_preserves_strict_tabulated_endpoints():
    k = np.logspace(-2, 2, 101)
    power = TabulatedPowerSpectrum(k, np.ones_like(k))
    model = IntegratedVarianceModel(
        power=power,
        window=SharpKWindow(),
        rho_mean=1.0,
        k_min=k[0],
        k_max=k[-1],
        n_k=1001,
    )
    assert np.isfinite(model.variance(1.0))


def test_variance_cache_is_content_addressed_and_round_trips(tmp_path):
    mass = np.logspace(4, 8, 9)
    variance = mass**-0.3
    key = variance_cache_key(
        "wdm:test",
        mass,
        backend_identifier="native-test",
        settings={"window": "sharp-k", "c": 2.5},
    )
    changed = variance_cache_key(
        "wdm:test",
        mass,
        backend_identifier="native-test",
        settings={"window": "sharp-k", "c": 2.6},
    )
    assert key != changed

    path = tmp_path / "variance.npz"
    save_variance_cache(path, key=key, mass=mass, variance=variance)
    loaded_mass, loaded_variance = load_variance_cache(path, expected_key=key)
    np.testing.assert_array_equal(loaded_mass, mass)
    np.testing.assert_array_equal(loaded_variance, variance)
    with pytest.raises(ValueError, match="key"):
        load_variance_cache(path, expected_key=changed)
