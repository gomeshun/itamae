"""Structural interfaces implemented by SASHIMI physical model packages.

Protocols in this package describe how scientific components communicate with
ITAMAE numerical machinery. They intentionally avoid selecting a particular
CDM, WDM, SIDM, or FDM prescription.
"""

from .backends import CosmologyBackend, UnitBackend
from .power import PowerSpectrum, WindowFunction
from .variance import VarianceModel

__all__ = [
    "CosmologyBackend",
    "PowerSpectrum",
    "UnitBackend",
    "VarianceModel",
    "WindowFunction",
]
