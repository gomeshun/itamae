import numpy as np
import pytest

from itamae.types import HostState, OrbitalState, ProfileParameters, SubhaloState


def test_host_and_profile_state_validate_shapes():
    """Keep mass, density and profile fields attached to the same population nodes.

    Notes
    -----
    Pass criteria: Aligned host arrays and an NFW subhalo state are accepted; a mismatched
    host-array length raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    values = np.ones(2)
    host = HostState(values, values, values, values, values, values, values)
    assert host.m200.shape == (2,)
    with pytest.raises(ValueError, match="common shape"):
        HostState(values, np.ones(3), values, values, values, values, values)

    profile = ProfileParameters("nfw", {"r_s": values, "rho_s": values})
    state = SubhaloState(values, profile, np.ones(2, dtype=bool), np.zeros(2, dtype=int))
    assert state.profile.name == "nfw"


def test_orbital_state_requires_representation():
    """Require an explicit meaning for a stored orbital radius.

    Notes
    -----
    Pass criteria: A radius-only state without a representation is rejected; an
    instantaneous-radius state is accepted.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="explicit orbital representation"):
        OrbitalState(radius=np.ones(2))
    state = OrbitalState(radius=np.ones(2), representation="instantaneous-radius")
    assert state.radius.shape == (2,)
