"""Protocols for physical components used by the population executor."""

from collections.abc import Mapping
from typing import Any, Protocol

import numpy as np

from itamae.types import AccretionBatch

from .variance import VarianceModel

PopulationState = Mapping[str, np.ndarray]
SurvivalSelection = np.ndarray | Mapping[str, np.ndarray]


class HostHistoryModel(Protocol):
    """Supply a host mass history and its redshift derivative."""

    def m200(self, host_reference: Any, z: Any, cosmology: Any) -> Any:
        """Return the host ``M200`` at redshift ``z``."""
        ...

    def dmvir_dz(self, host_reference: Any, z: Any, cosmology: Any) -> Any:
        """Return the derivative of host virial mass with redshift."""
        ...


class AccretionRateModel(Protocol):
    """Supply a differential accretion abundance."""

    def differential_number(
        self, m_acc: Any, z_acc: Any, host: Any, variance: VarianceModel
    ) -> Any:
        """Return the differential number of accreted objects."""
        ...


class ConcentrationModel(Protocol):
    """Supply the median concentration relation."""

    def median(self, m200: Any, z: Any, cosmology: Any) -> Any:
        """Return median concentration values."""
        ...


class InitialStructureModel(Protocol):
    """Assign structure to accretion nodes."""

    def assign(self, m200: Any, z: Any, concentration_nodes: Any, context: Any) -> Any:
        """Return named initial-structure arrays."""
        ...


class MassLossLaw(Protocol):
    """Supply a mass-loss right-hand side for a generic solver."""

    def rhs(self, state: Any, host_state: Any, orbital_state: Any = None) -> Any:
        """Return the state derivative."""
        ...


class ProfileEvolutionModel(Protocol):
    """Evolve an initial profile after accretion."""

    def evolve(self, initial_profile: Any, mass_history: Any, context: Any) -> Any:
        """Return the evolved profile representation."""
        ...


class SurvivalModel(Protocol):
    """Evaluate survival or disruption conditions."""

    def evaluate(self, state: Any, context: Any) -> Any:
        """Return a boolean survival mask or named masks."""
        ...


class PopulationInitializer(Protocol):
    """Initialize aligned per-node state for one accretion batch.

    This is an execution-stage protocol. Scientific choices used to construct
    the initial state remain owned by the variant implementation.
    """

    def initialize(self, batch: AccretionBatch, context: Any) -> Mapping[str, Any]:
        """Return named array-like values aligned with ``batch``."""
        ...


class PopulationEvolver(Protocol):
    """Evolve validated initial population state for one accretion batch."""

    def evolve(
        self,
        batch: AccretionBatch,
        initial: PopulationState,
        context: Any,
    ) -> Mapping[str, Any]:
        """Return named array-like evolved values aligned with ``batch``."""
        ...


class PopulationSurvivalSelector(Protocol):
    """Select one or more named survival views from evolved population state."""

    def select(
        self,
        batch: AccretionBatch,
        initial: PopulationState,
        evolved: PopulationState,
        context: Any,
    ) -> SurvivalSelection:
        """Return one mask or a mapping of named masks aligned with ``batch``."""
        ...


class CatalogColumnBuilder(Protocol):
    """Build final aligned catalog columns from validated pipeline stage state."""

    def build(
        self,
        batch: AccretionBatch,
        initial: PopulationState,
        evolved: PopulationState,
        survival: PopulationState,
        context: Any,
    ) -> Mapping[str, Any]:
        """Return named array-like final columns aligned with ``batch``."""
        ...


__all__ = [
    "AccretionRateModel",
    "CatalogColumnBuilder",
    "ConcentrationModel",
    "HostHistoryModel",
    "InitialStructureModel",
    "MassLossLaw",
    "PopulationEvolver",
    "PopulationInitializer",
    "PopulationState",
    "PopulationSurvivalSelector",
    "ProfileEvolutionModel",
    "SurvivalModel",
    "SurvivalSelection",
]
