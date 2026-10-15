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
    """Concatenate one aligned stage field along the population axis.

    Parameters
    ----------
    values : Mapping or array_like
        Named physical fields or aligned numerical values supplied by the caller; each
        population field must refer to the same nodes.
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.

    Returns
    -------
    numpy.ndarray
        Stage values joined along the population axis, retaining their order.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """
    arrays = tuple(np.asarray(value) for value in values)
    if not arrays:
        raise ValueError(f"Cannot concatenate empty field {name!r}.")
    return np.concatenate(arrays, axis=0)


def _metadata_difference(left: Any, right: Any, path: str = "metadata") -> str | None:
    """Find the first exact mismatch, including nested mappings and array values.

    Parameters
    ----------
    left : object
        Two metadata values to compare exactly before mixing results from different batches.
    right : object
        Two metadata values to compare exactly before mixing results from different batches.
    path : str or pathlib.Path
        Location of the scientific data archive; no pickle evaluation is used by the NPZ
        readers.

    Returns
    -------
    str or None
        Path to the first exact metadata difference, or None when compatible.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        for key in sorted(left.keys() | right.keys(), key=str):
            field = f"{path}.{key}"
            if key not in left or key not in right:
                return field
            mismatch = _metadata_difference(left[key], right[key], field)
            if mismatch is not None:
                return mismatch
        return None
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        if len(left) != len(right):
            return path
        for index, (a, b) in enumerate(zip(left, right, strict=True)):
            mismatch = _metadata_difference(a, b, f"{path}[{index}]")
            if mismatch is not None:
                return mismatch
        return None
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        return None if np.array_equal(left, right) else path
    return None if type(left) is type(right) and left == right else path


def _validate_batches(values: tuple[AccretionBatch, ...]) -> None:
    """Reject incompatible contracts before invoking any physical component.

    Parameters
    ----------
    values : Mapping or array_like
        Named physical fields or aligned numerical values supplied by the caller; each
        population field must refer to the same nodes.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """
    if not values:
        raise ValueError("At least one accretion batch is required.")
    for index, batch in enumerate(values):
        if not isinstance(batch, AccretionBatch):
            raise TypeError(f"batch {index} must be an AccretionBatch.")
        if batch.m200_acc.ndim != 1:
            raise ValueError(f"batch {index} requires one flat population axis.")
        mismatch = _metadata_difference(values[0].metadata, batch.metadata)
        if mismatch is not None:
            raise ValueError(f"batch {index} has incompatible {mismatch}.")
        for name in ("weight_host_history", "weight_orbit"):
            if (getattr(batch, name) is None) != (getattr(values[0], name) is None):
                raise ValueError(
                    f"batch {index}: optional factor {name!r} must exist in every batch."
                )


def _validate_stage(
    values: Mapping[str, Any], shape: tuple[int, ...], stage: str
) -> Mapping[str, np.ndarray]:
    """Validate one callback's aligned finite arrays.

    Parameters
    ----------
    values : Mapping or array_like
        Named physical fields or aligned numerical values supplied by the caller; each
        population field must refer to the same nodes.
    shape : tuple of int
        Required population or coordinate shape. Extra axes cannot be silently treated as
        additional objects.
    stage : str
        Population stage name used to report the location of a failed physical callback.

    Returns
    -------
    Mapping[str, numpy.ndarray]
        Named finite arrays with the required population shape.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """
    if not isinstance(values, Mapping):
        raise TypeError(f"{stage} stage must return a mapping of named arrays.")
    if stage == "columns" and not values:
        raise ValueError("The columns stage must return at least one physical column.")
    arrays = {}
    for name, value in values.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"{stage} field names must be non-empty strings.")
        array = np.asarray(value)
        if array.shape != shape:
            raise ValueError(f"{stage} field {name!r} has shape {array.shape}; expected {shape}.")
        if array.dtype.kind not in "buif":
            raise TypeError(f"{stage} field {name!r} must contain real numeric or boolean values.")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{stage} field {name!r} contains non-finite values.")
        arrays[name] = array
    return MappingProxyType(arrays)


def _validate_survival(value: Any, shape: tuple[int, ...]) -> Mapping[str, np.ndarray]:
    """Normalize one or more survival masks without changing physical columns.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.
    shape : tuple of int
        Required population or coordinate shape. Extra axes cannot be silently treated as
        additional objects.

    Returns
    -------
    Mapping[str, numpy.ndarray]
        Named boolean survival masks with the required population shape.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """
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
        if raw.dtype.kind not in "buif" or not np.all((raw == 0) | (raw == 1)):
            raise ValueError(
                f"survival mask {name!r} must contain only boolean or binary 0/1 values."
            )
        masks[name] = raw.astype(bool, copy=False)
    return MappingProxyType(masks)


def concatenate_accretion_batches(
    batches: Iterable[AccretionBatch],
) -> AccretionBatch:
    """Concatenate aligned accretion batches while preserving optional factors.

    Parameters
    ----------
    batches : Iterable[AccretionBatch]
        Accretion batches with the same physical metadata and weight-factor schema; order is
        preserved during concatenation.

    Returns
    -------
    AccretionBatch
        Combined accretion nodes retaining every independent weight factor.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """
    values = tuple(batches)
    _validate_batches(values)

    def optional_factor(name: str) -> np.ndarray | None:
        """Transport one optional population-weight factor without inventing missing probabilities.

        Parameters
        ----------
        name : str
            Name of the quantity or component included in an error or provenance record so an
            invalid scientific input can be located.

        Returns
        -------
        numpy.ndarray or None
            Joined optional factor, or None only when it is absent in every batch.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        factors = tuple(getattr(batch, name) for batch in values)
        if all(factor is None for factor in factors):
            return None
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
    """Hold a transported population before selecting a catalog state view.

    Parameters
    ----------
    batch : AccretionBatch
        Aligned accretion nodes with separate population, concentration and optional
        host/orbit factors.
    columns : callable or Mapping
        Model column builder, or aligned physical arrays it produced, kept separately from
        statistical weights.
    survival : callable or Mapping
        Model survival selector, or the named boolean masks it produced; survival enters the
        final weight once.
    diagnostics : Mapping, optional
        Recorded numerical settings or counters for the execution, separate from physical
        catalogue columns.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """

    batch: AccretionBatch
    columns: Mapping[str, np.ndarray]
    survival: Mapping[str, np.ndarray]
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate and freeze the transported result.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        columns = _validate_stage(self.columns, self.batch.m200_acc.shape, "columns")
        survival = _validate_survival(self.survival, self.batch.m200_acc.shape)
        object.__setattr__(self, "columns", columns)
        object.__setattr__(self, "survival", survival)
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))

    @property
    def weight_factors(self) -> Mapping[str, np.ndarray]:
        """Return generation-stage factors supplied by the accretion batch.

        Returns
        -------
        Mapping[str, numpy.ndarray]
            Generation-stage factors; the survival factor is attached only when a view is
            selected.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
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
        """Assemble one named survival view as a weighted catalog.

        Parameters
        ----------
        metadata : CatalogMetadata or Mapping
            Physical model, backend, unit and source identities needed to interpret or combine
            the nodes.
        view : str
            Name of the survival mask used in this catalogue, for example a CDM reference or
            SIDM selection.

        Returns
        -------
        WeightedSubhaloCatalog
            Selected survival view with explicit generation and survival factors.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
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
    named survival views for models such as SIDM. All batch metadata must match
    exactly; batch-local coordinates or labels belong in arrays or ``contexts``.
    Nodes and output fields use a flat population axis. An empty iterable is an
    error; a zero-length batch runs the stages and retains their output schema.

    Parameters
    ----------
    initialize : callable or PopulationInitializer
        Model stage assigning the initial physical structure to each accretion node.
    evolve : callable or PopulationEvolver
        Model stage evolving the initialized structure under its tidal and dark-matter laws.
    survival : callable or Mapping
        Model survival selector, or the named boolean masks it produced; survival enters the
        final weight once.
    columns : callable or Mapping
        Model column builder, or aligned physical arrays it produced, kept separately from
        statistical weights.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
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
        """Run all stages in batch order and return a transport result.

        Parameters
        ----------
        batches : Iterable[AccretionBatch]
            Accretion batches with the same physical metadata and weight-factor schema; order is
            preserved during concatenation.
        contexts : Iterable, optional
            One variant-owned context per accretion batch, in the same order; omitted entries
            are passed as None.
        diagnostics : Mapping, optional
            Recorded numerical settings or counters for the execution, separate from physical
            catalogue columns.

        Returns
        -------
        PopulationExecution
            Transported columns and named survival masks in input batch order.

        Notes
        -----
        Physical initialization, evolution, disruption and column definitions belong to the
        SASHIMI variant. The executor retains aligned nodes and multiplies each independent
        abundance/probability/survival factor once. It does not introduce a new accretion or
        halo-evolution prescription.
        """
        batch_values = tuple(batches)
        _validate_batches(batch_values)
        if contexts is None:
            context_values = (None,) * len(batch_values)
        else:
            context_values = tuple(contexts)
            if len(context_values) != len(batch_values):
                raise ValueError("contexts must have one entry per accretion batch.")

        column_values: list[Mapping[str, np.ndarray]] = []
        survival_values: list[Mapping[str, np.ndarray]] = []
        for index, (batch, context) in enumerate(zip(batch_values, context_values, strict=True)):
            shape = batch.m200_acc.shape

            def run(
                stage: str, callback: Callable[..., Any], *args: Any
            ) -> Mapping[str, np.ndarray]:
                """Run one physical stage and attach its batch identity to any failure.

                Parameters
                ----------
                stage : str
                    Population stage name used to report the location of a failed physical callback.
                callback : callable
                    Model callback for the selected execution stage; must return arrays aligned with the
                    batch.
                *args : tuple
                    Additional positional parameters forwarded unchanged to the supplied numerical or
                    physical callback.

                Returns
                -------
                Mapping[str, numpy.ndarray]
                    Validated stage arrays; callback errors retain the batch/stage location and original
                    cause.

                Notes
                -----
                Physical initialization, evolution, disruption and column definitions belong to the
                SASHIMI variant. The executor retains aligned nodes and multiplies each independent
                abundance/probability/survival factor once. It does not introduce a new accretion or
                halo-evolution prescription.
                """
                try:
                    output = callback(*args)
                    if stage == "survival":
                        return _validate_survival(output, shape)
                    return _validate_stage(output, shape, stage)
                except Exception as exc:
                    name = getattr(callback, "__qualname__", type(callback).__name__)
                    raise ValueError(f"batch {index}, {stage} stage ({name}): {exc}") from exc

            initial = run("initial", self.initialize, batch, context)
            evolved = run("evolved", self.evolve, batch, initial, context)
            survival = run("survival", self.survival, batch, initial, evolved, context)
            columns = run("columns", self.columns, batch, initial, evolved, survival, context)
            if column_values and set(columns) != set(column_values[0]):
                raise ValueError(f"batch {index}: columns stage fields differ from batch 0.")
            if survival_values and set(survival) != set(survival_values[0]):
                raise ValueError(f"batch {index}: survival stage views differ from batch 0.")
            column_values.append(columns)
            survival_values.append(survival)

        column_names = tuple(column_values[0])
        survival_names = tuple(survival_values[0])

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
    """Transport a population through explicit physical stages and retain factorized expected
    counts.

    Parameters
    ----------
    batches : Iterable[AccretionBatch]
        Accretion batches with the same physical metadata and weight-factor schema; order is
        preserved during concatenation.
    initialize : callable or PopulationInitializer
        Model stage assigning the initial physical structure to each accretion node.
    evolve : callable or PopulationEvolver
        Model stage evolving the initialized structure under its tidal and dark-matter laws.
    survival : callable or Mapping
        Model survival selector, or the named boolean masks it produced; survival enters the
        final weight once.
    columns : callable or Mapping
        Model column builder, or aligned physical arrays it produced, kept separately from
        statistical weights.
    contexts : Iterable, optional
        One variant-owned context per accretion batch, in the same order; omitted entries
        are passed as None.
    diagnostics : Mapping, optional
        Recorded numerical settings or counters for the execution, separate from physical
        catalogue columns.

    Returns
    -------
    PopulationExecution
        Transported columns and named survival masks in input batch order.

    Notes
    -----
    Physical initialization, evolution, disruption and column definitions belong to the
    SASHIMI variant. The executor retains aligned nodes and multiplies each independent
    abundance/probability/survival factor once. It does not introduce a new accretion or
    halo-evolution prescription.
    """
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
