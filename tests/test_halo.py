import numpy as np
import pytest
from hypothesis import given, strategies as st

from itamae.halo import NFWProfile, invert_nfw_mass_function, nfw_mass_function


@given(st.floats(min_value=0.0, max_value=1.0e4, allow_nan=False, allow_infinity=False))
def test_nfw_inverse_roundtrip(x):
    recovered = invert_nfw_mass_function(nfw_mass_function(x))
    assert recovered == pytest.approx(x, rel=1e-9, abs=1e-9)


def test_nfw_profile_mass_is_monotonic():
    profile = NFWProfile(r_s=0.02, rho_s=1.0e15)
    r = np.logspace(-5, 0, 100)
    assert np.all(np.diff(profile.enclosed_mass(r)) > 0.0)
