from dataclasses import dataclass
import numpy as np
from scipy.integrate import simpson, quad
from scipy.optimize import brentq

def normalize_pdf_q(q,pdf):
    q=np.asarray(q,float); p=np.asarray(pdf,float); norm=simpson(p,x=q)
    if norm<=0 or not np.isfinite(norm): raise ValueError("radial PDF has invalid normalization")
    return p/norm

@dataclass(frozen=True)
class RadialMeasure:
    q_nodes: np.ndarray
    weights: np.ndarray
    def __post_init__(self):
        if np.shape(self.q_nodes)!=np.shape(self.weights): raise ValueError("radial nodes and weights must match")
        if np.any(np.asarray(self.weights)<0): raise ValueError("radial weights must be nonnegative")

def turning_points(potential,E,L,rmin,rmax,nscan=2048):
    grid=np.geomspace(rmin,rmax,nscan); f=2*(E-potential(grid))-L**2/grid**2; roots=[]
    for a,b,fa,fb in zip(grid[:-1],grid[1:],f[:-1],f[1:]):
        if fa*fb<0: roots.append(brentq(lambda r:2*(E-potential(r))-L**2/r**2,a,b))
    if len(roots)<2: raise ValueError("bound orbit requires two turning points")
    return roots[0],roots[-1]

def radial_period(potential,E,L,rp,ra):
    mid=(rp+ra)/2; half=(ra-rp)/2
    value=quad(lambda th: half*np.sin(th)/np.sqrt(max(2*(E-potential(mid-half*np.cos(th)))-L**2/(mid-half*np.cos(th))**2,1e-300)),0,np.pi,limit=200)[0]
    return 2*value

def radial_shell_pdf(r,potential,E,L,rp,ra,period=None):
    r=np.asarray(r,float); tr=radial_period(potential,E,L,rp,ra) if period is None else period
    vr2=2*(E-potential(r))-L**2/r**2
    return np.where((r>=rp)&(r<=ra)&(vr2>0),2/(tr*np.sqrt(vr2)),0.0)
