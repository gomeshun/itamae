import numpy as np
import pytest

from itamae.types import HostState, OrbitalState, ProfileParameters, SubhaloState


def test_host_and_profile_state_validate_shapes():
    values = np.ones(2)
    host = HostState(values, values, values, values, values, values, values)
    assert host.m200.shape == (2,)
    with pytest.raises(ValueError, match="common shape"):
        HostState(values, np.ones(3), values, values, values, values, values)

    profile = ProfileParameters("nfw", {"r_s": values, "rho_s": values})
    state = SubhaloState(values, profile, np.ones(2, dtype=bool), np.zeros(2, dtype=int))
    assert state.profile.name == "nfw"


def test_orbital_state_requires_representation():
    with pytest.raises(ValueError, match="explicit orbital representation"):
        OrbitalState(radius=np.ones(2))
    state = OrbitalState(radius=np.ones(2), representation="instantaneous-radius")
    assert state.radius.shape == (2,)
