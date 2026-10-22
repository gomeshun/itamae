import numpy as np
import pytest

from itamae.units import NativeUnits


def test_native_units_reject_nonfinite():
    """Keep canonical units explicit when accepting or returning physical values.

    Notes
    -----
    Pass criteria: NaN mass, unknown physical dimensions, and an unsupported output-unit
    string are rejected.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError):
        NativeUnits().validate([1.0, np.nan], "mass")
    with pytest.raises(KeyError):
        NativeUnits().to_internal(1.0, "unknown")
    with pytest.raises(ValueError):
        NativeUnits().from_internal(1.0, "kpc")


def test_astropy_units_roundtrip():
    """Convert physical lengths without losing the factor between kpc and Mpc.

    Notes
    -----
    Pass criteria: 2000 kpc becomes 2 Mpc and converts back; dimensionless arrays are
    unchanged; seconds cannot be converted to mass.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    u = pytest.importorskip("astropy.units")
    from itamae.units.astropy import AstropyUnits

    backend = AstropyUnits()
    value = backend.to_internal(2.0e3 * u.kpc, "length")
    assert value == pytest.approx(2.0)
    output = backend.from_internal(2.0, u.kpc)
    assert output.to_value(u.kpc) == pytest.approx(2.0e3)
    dimensionless = backend.to_internal([0.5, 1.0] * u.dimensionless_unscaled, "dimensionless")
    np.testing.assert_array_equal(dimensionless, np.array([0.5, 1.0]))
    with pytest.raises(u.UnitConversionError):
        backend.to_internal(1.0 * u.s, "mass")
