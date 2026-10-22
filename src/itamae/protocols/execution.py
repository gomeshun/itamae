"""Protocols for physical components used by the population executor."""

from collections.abc import Mapping
from typing import Any, Protocol

import numpy as np

from itamae.types import AccretionBatch

from .variance import VarianceModel

PopulationState = Mapping[str, np.ndarray]
SurvivalSelection = np.ndarray | Mapping[str, np.ndarray]


class HostHistoryModel(Protocol):
    """Supply a host mass history and its redshift derivative.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def m200(self, host_reference: Any, z: Any, cosmology: Any) -> Any:
        """Return the host ``M200`` at redshift ``z``.

        Parameters
        ----------
        host_reference : object
            Reference host mass/state identifying the history being evaluated, with mass and
            epoch conventions supplied by the model.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.
        cosmology : CosmologyBackend
            Background provider for expansion, density, growth and time; its identity is part of
            the calculation provenance.

        Returns
        -------
        array_like
            Host M200 at each requested redshift, in the model mass unit.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...

    def dmvir_dz(self, host_reference: Any, z: Any, cosmology: Any) -> Any:
        """Return the derivative of host virial mass with redshift.

        Parameters
        ----------
        host_reference : object
            Reference host mass/state identifying the history being evaluated, with mass and
            epoch conventions supplied by the model.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.
        cosmology : CosmologyBackend
            Background provider for expansion, density, growth and time; its identity is part of
            the calculation provenance.

        Returns
        -------
        array_like
            Derivative of virial host mass per unit redshift, in the same mass unit.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class AccretionRateModel(Protocol):
    """Supply a differential accretion abundance.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def differential_number(
        self, m_acc: Any, z_acc: Any, host: Any, variance: VarianceModel
    ) -> Any:
        """Return the differential number of accreted objects.

        Parameters
        ----------
        m_acc : float or array_like
            Positive halo mass. Use solar masses with the canonical backend; a supplied variance
            prescription must use the same mass convention and factors of h.
        z_acc : float or array_like
            Dimensionless accretion redshift. History integration requires one scalar accretion
            epoch.
        host : object
            Host properties at the relevant epoch, in the same physical units as the subhalo
            model.
        variance : array_like
            Dimensionless linear variance S=sigma**2 aligned with the mass grid.

        Returns
        -------
        array_like
            Differential expected accretion count; the concrete model must declare the
            mass/redshift measure.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class ConcentrationModel(Protocol):
    """Supply the median concentration relation.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def median(self, m200: Any, z: Any, cosmology: Any) -> Any:
        """Return median concentration values.

        Parameters
        ----------
        m200 : array_like
            Halo mass within 200 times critical density, in solar masses under the canonical
            unit schema; aligned with the population coordinates.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.
        cosmology : CosmologyBackend
            Background provider for expansion, density, growth and time; its identity is part of
            the calculation provenance.

        Returns
        -------
        array_like
            Dimensionless median halo concentration, aligned with mass and redshift.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class InitialStructureModel(Protocol):
    """Assign structure to accretion nodes.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def assign(self, m200: Any, z: Any, concentration_nodes: Any, context: Any) -> Any:
        """Return named initial-structure arrays.

        Parameters
        ----------
        m200 : array_like
            Halo mass within 200 times critical density, in solar masses under the canonical
            unit schema; aligned with the population coordinates.
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.
        concentration_nodes : array_like
            Positive dimensionless halo concentration aligned with the accretion nodes; its
            radius definition (r200/r_s or rvir/r_s) must match the supplied mass and structural
            prescription.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        Mapping or object
            Initial structural parameters for each concentration/accretion node.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class MassLossLaw(Protocol):
    """Supply a mass-loss right-hand side for a generic solver.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def rhs(self, state: Any, host_state: Any, orbital_state: Any = None) -> Any:
        """Return the state derivative.

        Parameters
        ----------
        state : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.
        host_state : object
            Host properties at the relevant epoch, in the same physical units as the subhalo
            model.
        orbital_state : object, optional
            Optional orbit state used by the mass-loss prescription; omission does not imply a
            particular orbit distribution.

        Returns
        -------
        array_like or Mapping
            State derivative in the equation integration-coordinate convention.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class ProfileEvolutionModel(Protocol):
    """Evolve an initial profile after accretion.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def evolve(self, initial_profile: Any, mass_history: Any, context: Any) -> Any:
        """Return the evolved profile representation.

        Parameters
        ----------
        initial_profile : object
            Internal density profile at accretion, in the variant-defined parameterization.
        mass_history : array_like
            Bound-mass trajectory after accretion, evaluated on the time/redshift grid declared
            by context.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        object
            Evolved internal profile in the variant representation.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class SurvivalModel(Protocol):
    """Evaluate survival or disruption conditions.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def evaluate(self, state: Any, context: Any) -> Any:
        """Return a boolean survival mask or named masks.

        Parameters
        ----------
        state : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        array_like or Mapping
            Boolean survival decisions, possibly in several named physical views.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class PopulationInitializer(Protocol):
    """Initialize aligned per-node state for one accretion batch.

    This is an execution-stage protocol. Scientific choices used to construct
    the initial state remain owned by the variant implementation.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def initialize(self, batch: AccretionBatch, context: Any) -> Mapping[str, Any]:
        """Return named array-like values aligned with ``batch``.

        Parameters
        ----------
        batch : AccretionBatch
            Aligned accretion nodes with separate population, concentration and optional
            host/orbit factors.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        Mapping[str, array_like]
            Initial physical fields aligned with the accretion batch.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class PopulationEvolver(Protocol):
    """Evolve validated initial population state for one accretion batch.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def evolve(
        self,
        batch: AccretionBatch,
        initial: PopulationState,
        context: Any,
    ) -> Mapping[str, Any]:
        """Return named array-like evolved values aligned with ``batch``.

        Parameters
        ----------
        batch : AccretionBatch
            Aligned accretion nodes with separate population, concentration and optional
            host/orbit factors.
        initial : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        Mapping[str, array_like]
            Evolved physical fields aligned with the accretion batch.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class PopulationSurvivalSelector(Protocol):
    """Select one or more named survival views from evolved population state.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def select(
        self,
        batch: AccretionBatch,
        initial: PopulationState,
        evolved: PopulationState,
        context: Any,
    ) -> SurvivalSelection:
        """Return one mask or a mapping of named masks aligned with ``batch``.

        Parameters
        ----------
        batch : AccretionBatch
            Aligned accretion nodes with separate population, concentration and optional
            host/orbit factors.
        initial : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.
        evolved : bool
            Use present bound masses/profiles when true or accretion masses/profiles when false,
            under the observable population selection.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        array_like or Mapping
            One boolean/binary mask or several named survival masks.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        ...


class CatalogColumnBuilder(Protocol):
    """Build final aligned catalog columns from validated pipeline stage state.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    def build(
        self,
        batch: AccretionBatch,
        initial: PopulationState,
        evolved: PopulationState,
        survival: PopulationState,
        context: Any,
    ) -> Mapping[str, Any]:
        """Return named array-like final columns aligned with ``batch``.

        Parameters
        ----------
        batch : AccretionBatch
            Aligned accretion nodes with separate population, concentration and optional
            host/orbit factors.
        initial : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.
        evolved : bool
            Use present bound masses/profiles when true or accretion masses/profiles when false,
            under the observable population selection.
        survival : callable or Mapping
            Model survival selector, or the named boolean masks it produced; survival enters the
            final weight once.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        Mapping[str, array_like]
            Final physical columns aligned with the same node axis; weights remain separate.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
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
