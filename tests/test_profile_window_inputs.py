"""Non-real coordinates must not silently become real physical inputs."""

import warnings

import numpy as np
import pytest

from itamae.halo import NFWProfile, invert_nfw_mass_function, nfw_mass_function
from itamae.power import SharpKWindow, SphericalTopHatWindow


@pytest.mark.parametrize("value", [np.complex128(1 + 2j), np.array([1 + 0j]), "1.0", True])
@pytest.mark.parametrize(
    "evaluate",
    [
        nfw_mass_function,
        invert_nfw_mass_function,
        NFWProfile(1.0, 1.0).density,
        NFWProfile(1.0, 1.0).potential,
        NFWProfile(1.0, 1.0).enclosed_mass,
        SharpKWindow(),
        SphericalTopHatWindow(),
    ],
)
def test_nonreal_profile_and_window_coordinates_fail(value, evaluate):
    """Prevent lossy conversion of physical radii and filter arguments.

    Parameters
    ----------
    value : object
        Parametrized physical value or malformed input, as specified by this case.
    evaluate : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Complex, text and Boolean coordinates raise ValueError without warnings
    in NFW mass, inverse, density, potential and both windows.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(ValueError, match="real numeric"):
            evaluate(value)


@pytest.mark.parametrize("name", ["r_s", "rho_s"])
@pytest.mark.parametrize("value", [np.complex128(1 + 2j), "1.0", True])
def test_profile_parameters_do_not_accept_nonreal_values(name, value):
    """Reject non-real NFW scale radius and density before model construction.

    Parameters
    ----------
    name : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.
    value : object
        Parametrized physical value or malformed input, as specified by this case.

    Notes
    -----
    Pass criteria: Complex, text and Boolean r_s or rho_s raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="real numeric"):
        NFWProfile(**{**dict(r_s=1.0, rho_s=1.0), name: value})
