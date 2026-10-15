import numpy as np
import pytest
from hypothesis import given, strategies as st

from itamae.halo import NFWProfile, invert_nfw_mass_function, nfw_mass_function


@given(st.floats(min_value=0.0, max_value=1.0e4, allow_nan=False, allow_infinity=False))
def test_nfw_inverse_roundtrip(x):
    """Recover the dimensionless truncation radius from its NFW enclosed mass.

    Parameters
    ----------
    x : float
        Dimensionless radius r/r_s generated for the NFW inverse check.

    Notes
    -----
    Pass criteria: For generated radii 0..1e4, inversion recovers the radius with relative
    and absolute tolerance 1e-9.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    recovered = invert_nfw_mass_function(nfw_mass_function(x))
    assert recovered == pytest.approx(x, rel=1e-9, abs=1e-9)


def test_nfw_profile_mass_is_monotonic():
    """Keep an NFW halo's enclosed mass increasing with physical radius.

    Notes
    -----
    Pass criteria: The enclosed mass increases strictly over 100 radii from 1e-5 to 1 Mpc
    for the fixed positive profile.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    profile = NFWProfile(r_s=0.02, rho_s=1.0e15)
    r = np.logspace(-5, 0, 100)
    assert np.all(np.diff(profile.enclosed_mass(r)) > 0.0)
