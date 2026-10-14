import pytest

from itamae.backends import BackendConfig
from itamae.cosmology import NativeFlatLCDM
from itamae.units import NativeUnits


def test_backend_config_validates_and_exports_metadata():
    """Record the cosmology and unit implementation used for a physical calculation.

    Notes
    -----
    Pass criteria: Supported NumPy/native backends produce consistent metadata; selecting
    the unsupported JAX array backend raises ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    config = BackendConfig(cosmology=NativeFlatLCDM(), units=NativeUnits())
    assert config.identifier.startswith("array=numpy;cosmology=native-flatlcdm")
    assert config.metadata()["backend_identifier"] == config.identifier
    with pytest.raises(ValueError, match="NumPy"):
        BackendConfig(cosmology=NativeFlatLCDM(), units=NativeUnits(), array="jax")


def test_backend_config_rejects_incomplete_objects():
    """Fail before a calculation uses an object lacking cosmology or unit operations.

    Notes
    -----
    Pass criteria: Objects missing the CosmologyBackend or UnitBackend interfaces raise
    TypeError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(TypeError, match="CosmologyBackend"):
        BackendConfig(cosmology=object(), units=NativeUnits())
    with pytest.raises(TypeError, match="UnitBackend"):
        BackendConfig(cosmology=NativeFlatLCDM(), units=object())
