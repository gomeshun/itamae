import numpy as np
import pytest

from itamae.units import NativeUnits


def test_native_units_reject_nonfinite():
    with pytest.raises(ValueError):
        NativeUnits().validate([1.0, np.nan], "mass")


def test_astropy_units_roundtrip():
    u = pytest.importorskip("astropy.units")
    from itamae.units.astropy import AstropyUnits

    backend = AstropyUnits()
    value = backend.to_internal(2.0e3 * u.kpc, "length")
    assert value == pytest.approx(2.0)
    with pytest.raises(u.UnitConversionError):
        backend.to_internal(1.0 * u.s, "mass")
