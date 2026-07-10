import numpy as np

from itamae.protocols import VarianceModel
from itamae.variance import CallableVarianceModel


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
