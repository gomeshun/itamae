"""Protocols for physical components used by the population executor."""

from typing import Any, Protocol

from .variance import VarianceModel


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


__all__ = [
    "AccretionRateModel",
    "ConcentrationModel",
    "HostHistoryModel",
    "InitialStructureModel",
    "MassLossLaw",
    "ProfileEvolutionModel",
    "SurvivalModel",
]
