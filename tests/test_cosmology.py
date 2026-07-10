import numpy as np
import pytest

from itamae.cosmology import NativeFlatLCDM


def test_native_cosmology_basic_relations():
    cosmo = NativeFlatLCDM()
    assert cosmo.H(0.0) == pytest.approx(67.4)
    assert cosmo.growth_factor(0.0) == pytest.approx(1.0)
    assert np.all(cosmo.cosmic_time(np.array([0.0, 1.0])) > 0.0)
    assert cosmo.lookback_time(1.0) > 0.0


def test_colossus_backend_matches_order_of_magnitude():
    pytest.importorskip("colossus")
    from itamae.cosmology.colossus import ColossusCosmology

    native = NativeFlatLCDM()
    colossus = ColossusCosmology("planck18")
    z = np.array([0.0, 1.0, 3.0])
    assert np.allclose(native.H(z), colossus.H(z), rtol=0.03)
    assert np.allclose(native.growth_factor(z), colossus.growth_factor(z), rtol=0.04)
