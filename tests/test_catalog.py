import numpy as np
import pytest

from itamae.types import WeightedSubhaloCatalog


def test_catalog_factorized_weights_and_selection():
    catalog = WeightedSubhaloCatalog(
        columns={"mass": np.array([1.0, 2.0, 3.0])},
        weights={"base": np.array([1.0, 2.0, 3.0]), "survival": np.array([1.0, 0.5, 0.0])},
    )
    assert np.allclose(catalog.weight_final, [1.0, 1.0, 0.0])
    assert catalog.weighted_sum(catalog.columns["mass"]) == pytest.approx(3.0)
    assert catalog.select(np.array([True, False, True])).shape == (2,)


def test_catalog_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        WeightedSubhaloCatalog(columns={"x": np.ones(2)}, weights={"w": np.ones(3)})
