from dataclasses import dataclass, field
from typing import Any, Mapping
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

Array=np.ndarray

@dataclass(frozen=True)
class HostState:
    redshift:Array; time:Array; m200:Array; mvir:Array; r200:Array; rvir:Array; concentration:Array

@dataclass(frozen=True)
class AccretionBatch:
    m200_acc:Array; mvir_acc:Array; z_acc:Array; concentration_acc:Array; weight_base:Array; weight_concentration:Array
    metadata:Mapping[str,Any]=field(default_factory=dict)
    def __post_init__(self):
        shape=np.shape(self.m200_acc)
        for name in ("mvir_acc","z_acc","concentration_acc","weight_base","weight_concentration"):
            if np.shape(getattr(self,name))!=shape: raise ValueError(f"{name} has incompatible shape")

@dataclass
class SubhaloState:
    m_bound:Array; profile:Mapping[str,Array]; alive:Array; flags:Array; extra:Mapping[str,Array]=field(default_factory=dict)

@dataclass
class OrbitalState:
    energy:Array|None=None; angular_momentum:Array|None=None; radius:Array|None=None; radial_velocity:Array|None=None
    tangential_velocity:Array|None=None; pericenter:Array|None=None; apocenter:Array|None=None; phase:Array|None=None

@dataclass
class WeightedSubhaloCatalog:
    columns: Mapping[str,np.ndarray]
    weights: Mapping[str,np.ndarray]
    metadata: dict[str,Any]=field(default_factory=dict)
    def __post_init__(self):
        lengths={len(np.asarray(v)) for v in [*self.columns.values(),*self.weights.values()]}
        if len(lengths)>1: raise ValueError("all catalog columns and weights must have the same length")
    def __len__(self):
        if not self.columns and not self.weights: return 0
        source=self.columns if self.columns else self.weights
        return len(np.asarray(next(iter(source.values()))))
    @property
    def final_weight(self):
        if "weight_final" in self.weights: return np.asarray(self.weights["weight_final"],float)
        result=np.ones(len(self))
        for value in self.weights.values(): result*=np.asarray(value,float)
        return result
    def select(self,mask):
        m=np.asarray(mask,bool)
        return type(self)({k:np.asarray(v)[m] for k,v in self.columns.items()},{k:np.asarray(v)[m] for k,v in self.weights.items()},dict(self.metadata))
    def weighted_sum(self,values): return float(np.sum(np.asarray(values)*self.final_weight))
    def weighted_histogram(self,column,bins): return np.histogram(self.columns[column],bins=bins,weights=self.final_weight)
    def poisson_realization(self,rng=None):
        rng=np.random.default_rng() if rng is None else rng
        idx=np.repeat(np.arange(len(self)),rng.poisson(self.final_weight))
        return {k:np.asarray(v)[idx] for k,v in self.columns.items()}

def log_grid(vmin,vmax,n):
    if vmin<=0 or vmax<=vmin or n<2: raise ValueError("invalid logarithmic grid")
    return np.geomspace(vmin,vmax,n)

def redshift_grid(z_target,z_max,dz):
    if dz<=0 or z_max<=z_target: raise ValueError("invalid redshift grid")
    return np.arange(z_target+dz,z_max+0.5*dz,dz)

def gauss_hermite_lognormal(median,sigma_log10,n):
    if n<1: raise ValueError("n must be positive")
    x,w=hermgauss(n)
    return np.asarray(median)[...,None]*10**(np.sqrt(2)*sigma_log10*x),w/np.sqrt(np.pi)

def nfw_f(x):
    x=np.asarray(x,float); return np.log1p(x)-x/(1+x)

def inverse_nfw_f(y):
    arr=np.atleast_1d(np.asarray(y,float)); out=[]
    for val in arr:
        if val<0: raise ValueError("NFW enclosed-mass factor must be nonnegative")
        if val==0: out.append(0.0); continue
        hi=1.0
        while nfw_f(hi)<val: hi*=2
        out.append(brentq(lambda x:nfw_f(x)-val,0,hi,xtol=5e-15,rtol=1e-14))
    out=np.asarray(out); return float(out[0]) if np.ndim(y)==0 else out

@dataclass(frozen=True)
class NFWProfile:
    rs: float|np.ndarray; rhos: float|np.ndarray; G: float=4.30091e-6
    def density(self,r):
        x=np.asarray(r)/self.rs; return self.rhos/(x*(1+x)**2)
    def enclosed_mass(self,r): return 4*np.pi*np.asarray(self.rhos)*np.asarray(self.rs)**3*nfw_f(np.asarray(r)/self.rs)
    def potential(self,r):
        x=np.asarray(r)/self.rs; return -4*np.pi*self.G*self.rhos*self.rs**2*np.log1p(x)/x
    def vcirc(self,r): return np.sqrt(self.G*self.enclosed_mass(r)/np.asarray(r))

def radius_from_mass(mass,overdensity,rho_reference):
    return (3*np.asarray(mass)/(4*np.pi*overdensity*np.asarray(rho_reference)))**(1/3)

class ODEEvolutionSolver:
    def solve(self,rhs,y0,t_eval,**kwargs):
        t=np.asarray(t_eval,float)
        sol=solve_ivp(lambda x,y:np.atleast_1d(rhs(x,y)),(t[0],t[-1]),np.atleast_1d(y0).astype(float),t_eval=t,**kwargs)
        if not sol.success: raise RuntimeError(sol.message)
        return sol.y.squeeze()

class PerturbativeEvolutionSolver:
    """Successive Picard iteration with optional elementwise Shanks acceleration."""
    def solve(self,rhs,y0,t_grid,order=2,shanks=False):
        t=np.asarray(t_grid,float); values=[np.broadcast_to(np.asarray(y0,float),t.shape).copy()]
        for _ in range(order+1):
            prev=values[-1]; f=np.asarray(rhs(t,prev),float); cur=np.empty_like(prev); cur[0]=np.asarray(y0,float)
            cur[1:]=cur[0]+np.cumsum(0.5*(f[1:]+f[:-1])*np.diff(t)); values.append(cur)
        if shanks and len(values)>=3:
            a,b,c=values[-3:]; den=c-2*b+a
            return np.where(np.abs(den)>1e-14,a-(b-a)**2/den,c)
        return values[-1]

def build_accretion_batch(m200,z,concentration,base_weight,concentration_weight,mvir=None,metadata=None):
    arrays=np.broadcast_arrays(m200,z,concentration,base_weight,concentration_weight)
    m200,z,c,wb,wc=[np.asarray(a,float).reshape(-1) for a in arrays]
    mv=m200.copy() if mvir is None else np.broadcast_to(mvir,arrays[0].shape).reshape(-1)
    return AccretionBatch(m200,mv,z,c,wb,wc,metadata or {})
