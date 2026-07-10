import numpy as np
from itamae import ODEEvolutionSolver, PerturbativeEvolutionSolver, build_accretion_batch

def test_ode_exponential():
    t=np.linspace(0,1,50)
    y=ODEEvolutionSolver().solve(lambda t,y:-y,1.,t,rtol=1e-9,atol=1e-11)
    np.testing.assert_allclose(y,np.exp(-t),rtol=1e-6)

def test_perturbative_constant_rhs():
    t=np.linspace(0,1,20)
    y=PerturbativeEvolutionSolver().solve(lambda t,y:-np.ones_like(t),2.,t,order=2)
    np.testing.assert_allclose(y,2-t)

def test_batch_broadcast():
    b=build_accretion_batch(np.array([1,2]),0.5,10.,np.array([0.2,0.3]),1.)
    assert b.m200_acc.shape==(2,)
    np.testing.assert_allclose(b.weight_concentration,1)
