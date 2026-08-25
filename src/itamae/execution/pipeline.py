"""Shared population transport and weighted-catalog assembly."""

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import numpy as np

from itamae.types import AccretionBatch, CatalogMetadata, WeightedSubhaloCatalog

StageInitializer = Callable[[AccretionBatch, Any], Mapping[str, Any]]
StageEvolution = Callable[[AccretionBatch, Mapping[str, np.ndarray], Any], Mapping[str, Any]]
StageSurvival = Callable[
    [
        AccretionBatch,
        Mapping[str, np.ndarray],
        Mapping[str, np.ndarray],
        Any,
    ],
    Any,
]
StageColumns = Callable[
    [
        AccretionBatch,
        Mapping[str, np.ndarray],
        Mapping[str, np.ndarray],
        Mapping[str, np.ndarray],
        Any,
    ],
    Mapping[str, Any],
]


def _concatenate(values: Iterable[np.ndarray], name: str) -> np.ndarray:
    """Concatenate one aligned stage field along the population axis."""
    arrays = tuple(np.asarray(value) for value in values)
    if not arrays:
        raise ValueError(f"Cannot concatenate empty field {name!r}.")
    return np.concatenate(arrays, axis=0)


def _validate_stage(
    values: Mapping[str, Any], shape: tuple[int, ...], stage: str
) -> Mapping[str, np.ndarray]:
    """Validate one callback's aligned finite arrays."""
    if not isinstance(values, Mapping):
        raise TypeError(f"{stage} stage must return a mapping of named arrays.")
    arrays = {}
    for name, value in values.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"{stage} field names must be non-empty strings.")
        array = np.asarray(value)
        if array.shape != shape:
            raise ValueError(f"{stage} field {name!r} has shape {array.shape}; expected {shape}.")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{stage} field {name!r} contains non-finite values.")
        arrays[name] = array
    return MappingProxyType(arrays)


def _validate_survival(value: Any, shape: tuple[int, ...]) -> Mapping[str, np.ndarray]:
    """Normalize one or more survival masks without changing physical columns."""
    values = {"default": value} if not isinstance(value, Mapping) else dict(value)
    if not values:
        raise ValueError("The survival stage must return at least one named mask.")
    masks = {}
    for name, mask in values.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Survival view names must be non-empty strings.")
        raw = np.asarray(mask)
        if raw.shape != shape:
            raise ValueError(f"Survival mask {name!r} has shape {raw.shape}; expected {shape}.")
        if raw.dtype.kind != "b" and not np.all(np.isfinite(raw)):
            raise ValueError(f"Survival mask {name!r} contains non-finite values.")
        masks[name] = raw.astype(bool, copy=False)
    return MappingProxyType(masks)


def concatenate_accretion_batches(
    batches: Iterable[AccretionBatch],
) -> AccretionBatch:
    """Concatenate aligned accretion batches while preserving optional factors."""
    values = tuple(batches)
    if not values:
        raise ValueError("At least one accretion batch is required.")

    def optional_factor(name: str) -> np.ndarray | None:
        factors = tuple(getattr(batch, name) for batch in values)
        if all(factor is None for factor in factors):
            return None
        if any(factor is None for factor in factors):
            raise ValueError(f"Optional factor {name!r} must exist in every batch.")
        return _concatenate(factors, name)  # type: ignore[arg-type]

    return AccretionBatch(
        m200_acc=_concatenate((batch.m200_acc for batch in values), "m200_acc"),
        mvir_acc=_concatenate((batch.mvir_acc for batch in values), "mvir_acc"),
        z_acc=_concatenate((batch.z_acc for batch in values), "z_acc"),
        concentration_acc=_concatenate(
            (batch.concentration_acc for batch in values), "concentration_acc"
        ),
        weight_base=_concatenate((batch.weight_base for batch in values), "weight_base"),
        weight_concentration=_concatenate(
            (batch.weight_concentration for batch in values), "weight_concentration"
        ),
        weight_host_history=optional_factor("weight_host_history"),
        weight_orbit=optional_factor("weight_orbit"),
        metadata=dict(values[0].metadata),
    )


@dataclass(frozen=True, slots=True)
class PopulationExecution:
    """Hold a transported population before selecting a catalog state view."""

    batch: AccretionBatch
    columns: Mapping[str, np.ndarray]
    survival: Mapping[str, np.ndarray]
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate and freeze the transported result."""
        columns = _validate_stage(self.columns, self.batch.m200_acc.shape, "columns")
        survival = _validate_survival(self.survival, self.batch.m200_acc.shape)
        object.__setattr__(self, "columns", columns)
        object.__setattr__(self, "survival", survival)
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))

    @property
    def weight_factors(self) -> Mapping[str, np.ndarray]:
        """Return generation-stage factors supplied by the accretion batch."""
        values: dict[str, np.ndarray] = {
            "weight_base": self.batch.weight_base,
            "weight_concentration": self.batch.weight_concentration,
        }
        if self.batch.weight_host_history is not None:
            values["weight_host_history"] = self.batch.weight_host_history
        if self.batch.weight_orbit is not None:
            values["weight_orbit"] = self.batch.weight_orbit
        return MappingProxyType(values)

    def to_catalog(
        self,
        metadata: CatalogMetadata | Mapping[str, Any],
        *,
        view: str = "default",
    ) -> WeightedSubhaloCatalog:
        """Assemble one named survival view as a weighted catalog."""
        if view not in self.survival:
            raise ValueError(
                f"Unknown survival view {view!r}; available views are {sorted(self.survival)}."
            )
        weights = dict(self.weight_factors)
        weights["weight_survival"] = self.survival[view].astype(float)
        return WeightedSubhaloCatalog(
            columns=self.columns,
            weights=weights,
            metadata=metadata,
        )


@dataclass(frozen=True, slots=True)
class PopulationPipeline:
    """Execute model-supplied stages over aligned accretion batches.

    The executor owns ordering, shape validation, concatenation, and weight
    transport. Callbacks own all physical prescriptions and may return several
    named survival views for models such as SIDM.
    """

    initialize: StageInitializer
    evolve: StageEvolution
    survival: StageSurvival
    columns: StageColumns

    def execute(
        self,
        batches: Iterable[AccretionBatch],
        *,
        contexts: Iterable[Any] | None = None,
        diagnostics: Mapping[str, Any] | None = None,
    ) -> PopulationExecution:
        """Run all stages in batch order and return a transport result."""
        batch_values = tuple(batches)
        if not batch_values:
            raise ValueError("At least one accretion batch is required.")
        if contexts is None:
            context_values = (None,) * len(batch_values)
        else:
            context_values = tuple(contexts)
            if len(context_values) != len(batch_values):
                raise ValueError("contexts must have one entry per accretion batch.")

        column_values: list[Mapping[str, np.ndarray]] = []
        survival_values: list[Mapping[str, np.ndarray]] = []
        for batch, context in zip(batch_values, context_values, strict=True):
            shape = batch.m200_acc.shape
            initial = _validate_stage(self.initialize(batch, context), shape, "initial")
            evolved = _validate_stage(self.evolve(batch, initial, context), shape, "evolved")
            survival = _validate_survival(self.survival(batch, initial, evolved, context), shape)
            columns = _validate_stage(
                self.columns(batch, initial, evolved, survival, context),
                shape,
                "columns",
            )
            column_values.append(columns)
            survival_values.append(survival)

        column_names = tuple(column_values[0])
        if any(tuple(values) != column_names for values in column_values[1:]):
            raise ValueError("The columns stage must return the same fields for every batch.")
        survival_names = tuple(survival_values[0])
        if any(tuple(values) != survival_names for values in survival_values[1:]):
            raise ValueError("The survival stage must return the same views for every batch.")

        return PopulationExecution(
            batch=concatenate_accretion_batches(batch_values),
            columns={
                name: _concatenate((values[name] for values in column_values), name)
                for name in column_names
            },
            survival={
                name: _concatenate((values[name] for values in survival_values), name)
                for name in survival_names
            },
            diagnostics={} if diagnostics is None else diagnostics,
        )


def execute_population(
    batches: Iterable[AccretionBatch],
    *,
    initialize: StageInitializer,
    evolve: StageEvolution,
    survival: StageSurvival,
    columns: StageColumns,
    contexts: Iterable[Any] | None = None,
    diagnostics: Mapping[str, Any] | None = None,
) -> PopulationExecution:
    """Convenience wrapper around :class:`PopulationPipeline`."""
    return PopulationPipeline(
        initialize=initialize,
        evolve=evolve,
        survival=survival,
        columns=columns,
    ).execute(batches, contexts=contexts, diagnostics=diagnostics)


__all__ = [
    "PopulationExecution",
    "PopulationPipeline",
    "concatenate_accretion_batches",
    "execute_population",
]
