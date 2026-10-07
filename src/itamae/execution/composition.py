"""Explicit composition of variant-owned population execution stages."""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

from itamae.protocols.execution import (
    CatalogColumnBuilder,
    PopulationEvolver,
    PopulationInitializer,
    PopulationSurvivalSelector,
)
from itamae.types import AccretionBatch

from .pipeline import PopulationExecution, PopulationPipeline


@dataclass(frozen=True, slots=True)
class PopulationComponents:
    """Compose explicit model components through :class:`PopulationPipeline`.

    ITAMAE owns stage ordering, validation, concatenation, and weight transport.
    The supplied objects own the scientific prescriptions used at each stage.
    This adapter deliberately preserves the existing callback-driven pipeline:
    it converts component methods into those callbacks rather than introducing
    a second execution implementation.

    Parameters
    ----------
    initializer
        Variant-owned component that creates aligned initial state.
    evolver
        Variant-owned component that evolves the initial state.
    survival
        Variant-owned component that selects one or more survival views.
    columns
        Variant-owned component that constructs final catalog columns.
    """

    initializer: PopulationInitializer
    evolver: PopulationEvolver
    survival: PopulationSurvivalSelector
    columns: CatalogColumnBuilder

    def to_pipeline(self) -> PopulationPipeline:
        """Return the canonical callback executor wired to these components."""
        return PopulationPipeline(
            initialize=self.initializer.initialize,
            evolve=self.evolver.evolve,
            survival=self.survival.select,
            columns=self.columns.build,
        )

    def execute(
        self,
        batches: Iterable[AccretionBatch],
        *,
        contexts: Iterable[Any] | None = None,
        diagnostics: Mapping[str, Any] | None = None,
    ) -> PopulationExecution:
        """Execute the composed components through the canonical pipeline."""
        return self.to_pipeline().execute(
            batches,
            contexts=contexts,
            diagnostics=diagnostics,
        )


__all__ = ["PopulationComponents"]
