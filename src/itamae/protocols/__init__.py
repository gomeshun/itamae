"""Structural interfaces implemented by SASHIMI physical model packages.

Protocols in this package describe how scientific components communicate with
ITAMAE numerical machinery. They intentionally avoid selecting a particular
CDM, WDM, SIDM, or FDM prescription.
"""

from .backends import CosmologyBackend, UnitBackend
from .execution import (
    AccretionRateModel,
    CatalogColumnBuilder,
    ConcentrationModel,
    HostHistoryModel,
    InitialStructureModel,
    MassLossLaw,
    PopulationEvolver,
    PopulationInitializer,
    PopulationSurvivalSelector,
    ProfileEvolutionModel,
    SurvivalModel,
)
from .power import PowerSpectrum, WindowFunction
from .variance import VarianceModel

__all__ = [
    "AccretionRateModel",
    "CatalogColumnBuilder",
    "ConcentrationModel",
    "CosmologyBackend",
    "HostHistoryModel",
    "InitialStructureModel",
    "MassLossLaw",
    "PopulationEvolver",
    "PopulationInitializer",
    "PopulationSurvivalSelector",
    "ProfileEvolutionModel",
    "SurvivalModel",
    "PowerSpectrum",
    "UnitBackend",
    "VarianceModel",
    "WindowFunction",
]
