"""Public variance boundaries must not conceal invalid values or broadcasting."""

from dataclasses import replace

import numpy as np
import pytest

from itamae.power import SharpKWindow, SphericalTopHatWindow
from itamae.variance import CallableVarianceModel, IntegratedVarianceModel


class ConstantPower:
    identifier = "constant:test"

    def __call__(self, k):
        return np.ones_like(k)


def integrated(**options):
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
    with pytest.raises(ValueError):
        getattr(factory(), method)(mass, z)


@pytest.mark.parametrize(
    "options", [{"identifier": " "}, {"sigma_function": None}, {"derivative_function": 2}]
)
def test_callable_configuration_is_validated(options):
    with pytest.raises((ValueError, TypeError)):
        wrapped(**options)


@pytest.mark.parametrize("result", [np.nan, np.inf, 1 + 2j, [1, 2, 3]])
@pytest.mark.parametrize(
    "method,keyword", [("sigma", "sigma_function"), ("dvariance_dmass", "derivative_function")]
)
def test_callback_results_cannot_corrupt_or_expand_output(result, method, keyword):
    with pytest.raises(ValueError):
        getattr(wrapped(**{keyword: lambda m, z: result}), method)([1, 2])


def test_negative_sigma_and_square_overflow_fail():
    with pytest.raises(ValueError):
        wrapped(sigma_function=lambda m, z: -1.0).sigma(1.0)
    with pytest.raises(ValueError):
        wrapped(sigma_function=lambda m, z: 1.0e200).variance(1.0)


def test_scalar_callback_can_broadcast_but_invalid_input_never_reaches_callback():
    def forbidden(m, z):
        pytest.fail("invalid input reached callback")

    with pytest.raises(ValueError):
        wrapped(sigma_function=forbidden).sigma(-1.0)
    value = wrapped(sigma_function=lambda m, z: 2.0).sigma([[1.0], [2.0]], [0.0, 1.0, 2.0])
    np.testing.assert_array_equal(value, np.full((2, 3), 2.0))


@pytest.mark.parametrize("window", [SharpKWindow(), SphericalTopHatWindow()])
def test_complex_power_is_rejected_before_cast(window):
    class ComplexPower(ConstantPower):
        def __call__(self, k):
            return np.ones_like(k) * (1 + 1j)

    with pytest.raises(ValueError):
        integrated(power=ComplexPower(), window=window).variance(1.0)
    with pytest.raises(ValueError):
        integrated(power=ComplexPower(), window=window).dvariance_dmass(1.0)


@pytest.mark.parametrize("growth", [lambda z: 1 + 1j, lambda z: np.ones((2, 2)), lambda z: 1.0e200])
def test_invalid_growth_and_overflow_fail(growth):
    model = integrated(growth_function=growth, growth_identifier="bad:test")
    with pytest.raises(ValueError):
        model.variance([1.0, 2.0])


def test_invalid_window_shape_is_not_silently_broadcast():
    class BadWindow:
        identifier = "bad-window"

        def __call__(self, x):
            return np.ones((3, 1))

    with pytest.raises(ValueError):
        integrated(window=BadWindow()).variance([1.0, 2.0, 3.0])


def test_complex_integration_knots_and_configuration_fail():
    power = ConstantPower()
    power.integration_breakpoints = np.array([0.1 + 1j, 1 + 1j])
    with pytest.raises(ValueError):
        _ = integrated(power=power).identifier
    with pytest.raises(ValueError):
        replace(integrated(), rho_mean=np.complex128(1 + 2j))


def test_valid_broadcast_empty_and_derivative_signs_are_preserved():
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
