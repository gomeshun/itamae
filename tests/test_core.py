import numpy as np, pytest
from itamae import *

def test_grids_and_quadrature():
    assert len(log_grid(1,100,3))==3
    nodes,w=gauss_hermite_lognormal(10,0.2,5)
    assert np.sum(w)==pytest.approx(1)
    assert nodes.shape==(5,)

def test_nfw_inverse():
    x=np.geomspace(1e-4,100,30)
    np.testing.assert_allclose(inverse_nfw_f(nfw_f(x)),x,rtol=1e-10)

def test_nfw_mass():
    p=NFWProfile(1.,2.)
    assert p.enclosed_mass(2)>p.enclosed_mass(1)>0

def test_catalog():
    c=WeightedSubhaloCatalog({"m":np.array([1.,2.])},{"weight_base":np.array([2.,3.]),"weight_survival":np.array([1.,0.5])})
    np.testing.assert_allclose(c.final_weight,[2,1.5])
    assert c.weighted_sum(c.columns["m"])==pytest.approx(5)
    assert len(c.select([True,False]))==1
