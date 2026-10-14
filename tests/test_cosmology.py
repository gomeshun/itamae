import numpy as np
import pytest

from itamae.cosmology import NativeFlatLCDM


def test_native_cosmology_basic_relations():
    """Check the present-day normalization and sign of the native cosmological clock.

    Notes
    -----
    Pass criteria: H(0)=67.4 km/s/Mpc, D(0)=1, cosmic times at z=0 and 1 are positive, and
    the z=1 lookback time is positive.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    cosmo = NativeFlatLCDM()
    assert cosmo.H(0.0) == pytest.approx(67.4)
    assert cosmo.growth_factor(0.0) == pytest.approx(1.0)
    assert np.all(cosmo.cosmic_time(np.array([0.0, 1.0])) > 0.0)
    assert cosmo.lookback_time(1.0) > 0.0


def test_colossus_backend_matches_order_of_magnitude():
    """Detect gross background mismatches against the optional Colossus implementation.

    Notes
    -----
    Pass criteria: At z=0,1,3 the two Planck-like backgrounds agree within 3 percent in H
    and 4 percent in growth.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Only H and growth are compared, not critical/matter density. The current Colossus
    wrapper multiplies its Msun*h^2/kpc^3 density by 1e9 without applying h^2; this test
    cannot detect that canonical-unit mismatch.
    """
    pytest.importorskip("colossus")
    from itamae.cosmology.colossus import ColossusCosmology

    native = NativeFlatLCDM()
    colossus = ColossusCosmology("planck18")
    z = np.array([0.0, 1.0, 3.0])
    assert np.allclose(native.H(z), colossus.H(z), rtol=0.03)
    assert np.allclose(native.growth_factor(z), colossus.growth_factor(z), rtol=0.04)
