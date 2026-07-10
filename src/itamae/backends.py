from dataclasses import dataclass
from typing import Any, Protocol
import numpy as np
from scipy.integrate import quad

@dataclass(frozen=True)
class BackendConfig:
    cosmology: Any
    units: Any
    array_namespace: str = "numpy"
    def metadata(self) -> dict[str, str]:
        return {"cosmology_backend": type(self.cosmology).__name__, "unit_backend": type(self.units).__name__, "array_namespace": self.array_namespace}

class UnitBackend(Protocol):
    def to_value(self, value: Any, unit: str) -> np.ndarray: ...
    def attach(self, value: Any, unit: str) -> Any: ...

class NativeUnits:
    """Float-array backend. Values are assumed to already use requested units."""
    def to_value(self, value: Any, unit: str) -> np.ndarray:
        return np.asarray(value, dtype=float)
    def attach(self, value: Any, unit: str) -> np.ndarray:
        return np.asarray(value, dtype=float)

class AstropyUnits:
    def __init__(self):
        try:
            import astropy.units as u
        except ImportError as exc:
            raise ImportError("AstropyUnits requires the 'astropy' extra") from exc
        self._u = u
    def to_value(self, value: Any, unit: str) -> np.ndarray:
        target = self._u.Unit(unit)
        if hasattr(value, "to_value"):
            return np.asarray(value.to_value(target), dtype=float)
        return np.asarray(value, dtype=float)
    def attach(self, value: Any, unit: str) -> Any:
        return np.asarray(value) * self._u.Unit(unit)

@dataclass(frozen=True)
class NativeFlatLCDM:
    omega_m: float = 0.315
    h: float = 0.674
    G: float = 4.30091e-9
    @property
    def omega_l(self): return 1.0 - self.omega_m
    @property
    def H0(self): return 100.0 * self.h
    def E(self, z):
        z=np.asarray(z,float); return np.sqrt(self.omega_m*(1+z)**3+self.omega_l)
    def H(self,z): return self.H0*self.E(z)
    def rho_crit(self,z): return 3*self.H(z)**2/(8*np.pi*self.G)
    def growth_factor(self,z):
        z=np.asarray(z,float); omz=self.omega_m*(1+z)**3/self.E(z)**2; olz=self.omega_l/self.E(z)**2
        g=2.5*omz/(omz**(4/7)-olz+(1+omz/2)*(1+olz/70))
        g0=2.5*self.omega_m/(self.omega_m**(4/7)-self.omega_l+(1+self.omega_m/2)*(1+self.omega_l/70))
        return g/(g0*(1+z))
    def cosmic_time(self,z):
        conv=977.7922216807892
        arr=np.atleast_1d(np.asarray(z,float)); vals=[quad(lambda zp: 1/((1+zp)*self.H(zp)), zi, np.inf)[0]*conv for zi in arr]
        out=np.asarray(vals); return float(out[0]) if np.ndim(z)==0 else out

class ColossusCosmology:
    def __init__(self, name="planck18", params=None):
        try:
            from colossus.cosmology import cosmology
        except ImportError as exc:
            raise ImportError("ColossusCosmology requires the 'colossus' extra") from exc
        self._module=cosmology
        self._cosmo = cosmology.setCosmology(name, params) if params is not None else cosmology.setCosmology(name)
    def _activate(self):
        if self._module.getCurrent() is not self._cosmo: self._module.setCurrent(self._cosmo)
    def H(self,z):
        self._activate(); arr=np.asarray(z,float); out=np.asarray(self._cosmo.Hz(arr),float); return float(out) if arr.ndim==0 else out
    def rho_crit(self,z):
        self._activate(); arr=np.asarray(z,float); out=np.asarray(self._cosmo.rho_c(arr),float)*1e9; return float(out) if arr.ndim==0 else out
    def growth_factor(self,z):
        self._activate(); arr=np.asarray(z,float); out=np.asarray(self._cosmo.growthFactor(arr),float); return float(out) if arr.ndim==0 else out
    def cosmic_time(self,z):
        self._activate(); arr=np.asarray(z,float); out=np.asarray(self._cosmo.age(arr),float); return float(out) if arr.ndim==0 else out
