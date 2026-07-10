import importlib.util, numpy as np, pytest
from itamae import NativeFlatLCDM, NativeUnits, AstropyUnits, BackendConfig

def test_native_cosmology():
    c=NativeFlatLCDM(); assert c.H(0)==pytest.approx(67.4); assert c.growth_factor(0)==pytest.approx(1); assert c.rho_crit(0)>0; assert c.cosmic_time(0)>13

def test_backend_metadata():
    cfg=BackendConfig(NativeFlatLCDM(),NativeUnits()); assert cfg.metadata()["cosmology_backend"]=="NativeFlatLCDM"

@pytest.mark.skipif(importlib.util.find_spec("astropy") is None,reason="optional astropy")
def test_astropy_units():
    import astropy.units as u
    b=AstropyUnits(); np.testing.assert_allclose(b.to_value(2*u.kpc,"pc"),2000); assert b.attach(2,"kpc").unit==u.kpc

@pytest.mark.skipif(importlib.util.find_spec("colossus") is None,reason="optional colossus")
def test_colossus_backend_close_to_native():
    from itamae import ColossusCosmology
    native=NativeFlatLCDM(); col=ColossusCosmology("planck18"); np.testing.assert_allclose(col.H([0,1]),native.H([0,1]),rtol=0.03)
