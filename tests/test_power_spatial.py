import numpy as np, pytest
from scipy.integrate import simpson
from itamae import TabulatedPowerSpectrum, WDMTransferFunction, VarianceIntegrator
from itamae import normalize_pdf_q, turning_points, radial_period, radial_shell_pdf

def test_power_and_variance():
    k=np.geomspace(1e-4,1e4,1000)
    p=TabulatedPowerSpectrum(k,k*np.exp(-k/10))
    tf=WDMTransferFunction(0.1)
    assert tf(100)<tf(1)
    v=VarianceIntegrator(lambda x:p(x)*tf(x)**2,rho_mean=1).variance([1.,10.])
    assert np.all(np.isfinite(v)) and np.all(v>0)

def test_radial_normalization():
    q=np.linspace(0.01,1,1000)
    p=normalize_pdf_q(q,q*(1-q))
    assert simpson(p,x=q)==pytest.approx(1)

def test_kepler_orbit_kernel():
    pot=lambda r:-1/np.asarray(r); E=-0.5; L=np.sqrt(0.75)
    rp,ra=turning_points(pot,E,L,1e-3,10)
    assert rp==pytest.approx(0.5,rel=1e-4)
    assert ra==pytest.approx(1.5,rel=1e-4)
    tr=radial_period(pot,E,L,rp,ra)
    r=np.linspace(rp+1e-6,ra-1e-6,20000)
    pdf=radial_shell_pdf(r,pot,E,L,rp,ra,tr)
    assert simpson(pdf,x=r)==pytest.approx(1,rel=2e-2)
