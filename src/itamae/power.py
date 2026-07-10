from dataclasses import dataclass
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simpson

class TabulatedPowerSpectrum:
    def __init__(self,k,pk):
        k=np.asarray(k,float); pk=np.asarray(pk,float)
        if np.any(k<=0) or np.any(pk<0) or np.any(np.diff(k)<=0): raise ValueError("invalid power spectrum table")
        self.k=k; self.pk=pk
        self._interp=interp1d(np.log(k),np.log(np.maximum(pk,np.finfo(float).tiny)),bounds_error=True)
    def __call__(self,k): return np.exp(self._interp(np.log(k)))

@dataclass(frozen=True)
class WDMTransferFunction:
    alpha: float
    nu: float=1.12
    def __call__(self,k): return (1+(self.alpha*np.asarray(k))**(2*self.nu))**(-5/self.nu)

@dataclass(frozen=True)
class FDMTransferFunction:
    k_half: float
    def __call__(self,k):
        x=1.61*np.asarray(k)/self.k_half
        return np.cos(x**3)/(1+x**8)

class VarianceIntegrator:
    def __init__(self,power,window="tophat",rho_mean=1.0,sharp_k_c=2.5):
        if window not in {"tophat","sharp-k"}: raise ValueError("unknown window")
        self.power=power; self.window=window; self.rho_mean=rho_mean; self.sharp_k_c=sharp_k_c
    def _radius(self,m): return (3*np.asarray(m)/(4*np.pi*self.rho_mean))**(1/3)
    def _w(self,x):
        if self.window=="sharp-k": return (x<=1).astype(float)
        out=np.ones_like(x,float); mask=np.abs(x)>1e-5; xm=x[mask]
        out[mask]=3*(np.sin(xm)-xm*np.cos(xm))/xm**3
        return out
    def variance(self,mass,kmin=1e-4,kmax=1e4,n=2048):
        scalar=np.ndim(mass)==0
        mass=np.atleast_1d(np.asarray(mass,float)); k=np.geomspace(kmin,kmax,n); r=self._radius(mass)
        if self.window=="sharp-k": r=r/self.sharp_k_c
        x=r[:,None]*k[None,:]
        integrand=k[None,:]**3*self.power(k)[None,:]*self._w(x)**2/(2*np.pi**2)
        out=simpson(integrand,x=np.log(k),axis=1)
        return float(out[0]) if scalar else out
