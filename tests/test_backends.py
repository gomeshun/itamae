import pytest

from itamae.backends import BackendConfig
from itamae.cosmology import NativeFlatLCDM
from itamae.units import NativeUnits


def test_backend_config_validates_and_exports_metadata():
    config = BackendConfig(cosmology=NativeFlatLCDM(), units=NativeUnits())
    assert config.identifier.startswith("array=numpy;cosmology=native-flatlcdm")
    assert config.metadata()["backend_identifier"] == config.identifier
    with pytest.raises(ValueError, match="NumPy"):
        BackendConfig(cosmology=NativeFlatLCDM(), units=NativeUnits(), array="jax")


def test_backend_config_rejects_incomplete_objects():
    with pytest.raises(TypeError, match="CosmologyBackend"):
        BackendConfig(cosmology=object(), units=NativeUnits())
    with pytest.raises(TypeError, match="UnitBackend"):
        BackendConfig(cosmology=NativeFlatLCDM(), units=object())
