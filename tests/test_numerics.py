import numpy as np
import pytest
from hypothesis import given, strategies as st

from itamae.numerics import gauss_hermite_lognormal


@given(st.floats(min_value=0.1, max_value=100.0), st.floats(min_value=0.0, max_value=0.5))
def test_lognormal_quadrature_is_positive_and_normalized(median, sigma):
    nodes, weights = gauss_hermite_lognormal(median, sigma, order=7)
    assert np.all(nodes > 0.0)
    assert weights.sum() == pytest.approx(1.0)


def test_lognormal_invalid_inputs():
    with pytest.raises(ValueError):
        gauss_hermite_lognormal(0.0, 0.1)
