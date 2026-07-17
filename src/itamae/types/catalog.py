"""Versioned weighted-subhalo catalog representation.

Catalogs store deterministic population nodes and independent, explicitly
named weight factors. Physical model packages own the column definitions while
ITAMAE owns shape validation, weight semantics, and reproducibility metadata.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import os
from pathlib import Path
import tempfile
from types import MappingProxyType
from typing import Any, ClassVar, Mapping, Sequence

import numpy as np

CATALOG_SCHEMA_VERSION = "1.0"
CANONICAL_WEIGHT_FACTORS = (
    "weight_base",
    "weight_host_history",
    "weight_concentration",
    "weight_orbit",
    "weight_survival",
)
_REQUIRED_METADATA = ("schema_version", "model_identifier", "backend_identifier")
_RESERVED_METADATA = {*_REQUIRED_METADATA, "source_identifier"}


@dataclass(frozen=True, slots=True)
class CatalogMetadata:
    """Describe the model and computational backend that produced a catalog.

    Parameters
    ----------
    model_identifier
        Stable identifier for the SASHIMI physical-model composition.
    backend_identifier
        Stable identifier from :class:`itamae.backends.BackendConfig`.
    source_identifier
        Optional provenance identifier, such as a git commit or fixture name.
    schema_version
        ITAMAE catalog-schema version.
    extra
        Additional JSON-compatible provenance fields.
    """

    model_identifier: str
    backend_identifier: str
    source_identifier: str | None = None
    schema_version: str = CATALOG_SCHEMA_VERSION
    extra: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate stable identifiers and freeze additional metadata."""
        for name in ("model_identifier", "backend_identifier", "schema_version"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string.")
        if self.source_identifier is not None and (
            not isinstance(self.source_identifier, str) or not self.source_identifier.strip()
        ):
            raise ValueError("source_identifier must be None or a non-empty string.")
        collisions = sorted(_RESERVED_METADATA.intersection(self.extra))
        if collisions:
            raise ValueError(
                f"CatalogMetadata.extra may not override reserved fields: {collisions}."
            )
        object.__setattr__(self, "extra", MappingProxyType(dict(self.extra)))

    def as_mapping(self) -> Mapping[str, Any]:
        """Return a flat immutable mapping suitable for serialization."""
        values = {
            "schema_version": self.schema_version,
            "model_identifier": self.model_identifier,
            "backend_identifier": self.backend_identifier,
            **dict(self.extra),
        }
        if self.source_identifier is not None:
            values["source_identifier"] = self.source_identifier
        return MappingProxyType(values)


@dataclass(frozen=True, slots=True)
class WeightedSubhaloCatalog:
    """Store aligned columns and factorized statistical weights.

    Parameters
    ----------
    columns
        Mapping of column names to arrays sharing one shape.
    weights
        Mapping of weight-factor names to nonnegative arrays with the same
        shape. ``weight_base`` is required. Other names must begin with
        ``weight_``; the derived name ``weight_final`` is reserved.
    metadata
        Reproducibility metadata. A :class:`CatalogMetadata` object is
        preferred; mappings must provide ``schema_version``,
        ``model_identifier``, and ``backend_identifier``.
    """

    canonical_weight_factors: ClassVar[tuple[str, ...]] = CANONICAL_WEIGHT_FACTORS

    columns: Mapping[str, np.ndarray]
    weights: Mapping[str, np.ndarray]
    metadata: CatalogMetadata | Mapping[str, Any]

    def __post_init__(self) -> None:
        """Validate shapes, weight semantics, and reproducibility metadata."""
        if not self.columns:
            raise ValueError("A catalog must contain at least one physical column.")
        invalid_columns = [
            name for name in self.columns if not isinstance(name, str) or not name.strip()
        ]
        if invalid_columns:
            raise ValueError(f"Catalog column names must be non-empty strings: {invalid_columns}.")
        if "weight_base" not in self.weights:
            raise ValueError("Catalog weights must include the canonical 'weight_base' factor.")
        invalid_names = [
            name
            for name in self.weights
            if not isinstance(name, str) or not name.startswith("weight_") or name == "weight_final"
        ]
        if invalid_names:
            raise ValueError(
                "Weight names must start with 'weight_' and may not be 'weight_final'; "
                f"got {invalid_names}."
            )

        columns = {name: np.asarray(value) for name, value in self.columns.items()}
        weights = {name: np.asarray(value, dtype=float) for name, value in self.weights.items()}
        shapes = {value.shape for value in [*columns.values(), *weights.values()]}
        if len(shapes) != 1:
            raise ValueError(f"All catalog arrays must have the same shape; got {shapes}.")
        for name, value in weights.items():
            if not np.all(np.isfinite(value)):
                raise ValueError(f"Weight factor {name!r} contains non-finite values.")
            if np.any(value < 0.0):
                raise ValueError(f"Weight factor {name!r} contains negative values.")

        metadata = (
            self.metadata.as_mapping()
            if isinstance(self.metadata, CatalogMetadata)
            else MappingProxyType(dict(self.metadata))
        )
        missing = [name for name in _REQUIRED_METADATA if not metadata.get(name)]
        if missing:
            raise ValueError(f"Catalog metadata is missing required fields: {missing}.")
        invalid_metadata = [
            name
            for name in _REQUIRED_METADATA
            if not isinstance(metadata[name], str) or not metadata[name].strip()
        ]
        if invalid_metadata:
            raise ValueError(
                f"Catalog metadata fields must be non-empty strings: {invalid_metadata}."
            )
        if metadata["schema_version"] != CATALOG_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported catalog schema {metadata['schema_version']!r}; "
                f"expected {CATALOG_SCHEMA_VERSION!r}."
            )

        object.__setattr__(self, "columns", MappingProxyType(columns))
        object.__setattr__(self, "weights", MappingProxyType(weights))
        object.__setattr__(self, "metadata", metadata)

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the common catalog shape."""
        return next(iter(self.columns.values())).shape

    def __len__(self) -> int:
        """Return the size of the leading catalog axis."""
        return self.shape[0] if self.shape else 1

    @property
    def weight_final(self) -> np.ndarray:
        """Multiply every independent weight factor."""
        result = np.ones(self.shape, dtype=float)
        for value in self.weights.values():
            result *= value
        return result

    def select(self, mask: Any) -> WeightedSubhaloCatalog:
        """Return a catalog subset while preserving metadata."""
        return type(self)(
            columns={name: value[mask] for name, value in self.columns.items()},
            weights={name: value[mask] for name, value in self.weights.items()},
            metadata=self.metadata,
        )

    def weighted_sum(self, values: Any) -> float:
        """Return the weighted sum of an array aligned with the catalog."""
        array = np.asarray(values, dtype=float)
        if array.shape != self.shape:
            raise ValueError("Values must have the catalog shape.")
        return float(np.sum(array * self.weight_final))

    def weighted_histogram(
        self, column: str, bins: Any, **kwargs: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return a histogram weighted by ``weight_final``.

        Parameters
        ----------
        column
            Name of the catalog column to histogram.
        bins
            Bin specification forwarded to :func:`numpy.histogram`.
        **kwargs
            Additional arguments forwarded to :func:`numpy.histogram`.
        """
        if "weights" in kwargs:
            raise TypeError("weighted_histogram determines weights from the catalog.")
        return np.histogram(self.columns[column], bins=bins, weights=self.weight_final, **kwargs)

    def poisson_realization(self, rng: np.random.Generator) -> Mapping[str, np.ndarray]:
        """Draw a Poisson realization of the weighted catalog.

        Parameters
        ----------
        rng
            Explicit random-number generator. Requiring a generator avoids
            hidden global random state.

        Returns
        -------
        Mapping[str, numpy.ndarray]
            Realized physical columns. Deterministic population weights are not
            copied into the stochastic catalog.
        """
        if not isinstance(rng, np.random.Generator):
            raise TypeError("rng must be an explicit numpy.random.Generator.")
        if len(self.shape) != 1:
            raise ValueError("Poisson realization currently requires a one-dimensional catalog.")
        multiplicity = rng.poisson(self.weight_final)
        indices = np.repeat(np.arange(len(self)), multiplicity)
        return MappingProxyType({name: value[indices] for name, value in self.columns.items()})

    def to_npz(self, path: str | os.PathLike[str]) -> None:
        """Atomically serialize the catalog without Python pickle objects.

        Parameters
        ----------
        path
            Destination ``.npz`` path. Metadata must be finite
            JSON-compatible data.

        Notes
        -----
        Physical arrays retain their NumPy dtypes. Weight arrays are already
        canonical floating-point values. A manifest maps names to numbered
        archive fields so model-specific names never become executable object
        payloads.
        """
        metadata = (
            self.metadata.as_mapping()
            if isinstance(self.metadata, CatalogMetadata)
            else self.metadata
        )
        manifest = {
            "archive_schema": "itamae-weighted-catalog-npz:1.0",
            "catalog_metadata": dict(metadata),
            "columns": list(self.columns),
            "weights": list(self.weights),
        }
        try:
            manifest_json = json.dumps(
                manifest,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
        except (TypeError, ValueError) as exc:
            raise ValueError("Catalog metadata must be finite JSON-compatible data.") from exc

        payload = {"manifest_json": np.asarray(manifest_json)}
        payload.update(
            {f"column_{index}": value for index, value in enumerate(self.columns.values())}
        )
        payload.update(
            {f"weight_{index}": value for index, value in enumerate(self.weights.values())}
        )
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary_name: str | None = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=destination.parent,
                prefix=f".{destination.name}.",
                suffix=".npz",
                delete=False,
            ) as temporary:
                temporary_name = temporary.name
                np.savez_compressed(temporary, **payload)  # type: ignore[arg-type]
            os.replace(temporary_name, destination)
        finally:
            if temporary_name is not None and os.path.exists(temporary_name):
                os.unlink(temporary_name)

    @classmethod
    def from_npz(cls, path: str | os.PathLike[str]) -> WeightedSubhaloCatalog:
        """Load and validate a catalog produced by :meth:`to_npz`."""
        with np.load(Path(path), allow_pickle=False) as archive:
            try:
                manifest = json.loads(str(archive["manifest_json"]))
            except (KeyError, json.JSONDecodeError, TypeError) as exc:
                raise ValueError("Catalog archive has an invalid manifest.") from exc
            if manifest.get("archive_schema") != "itamae-weighted-catalog-npz:1.0":
                raise ValueError("Catalog archive schema is unsupported.")
            column_names = manifest.get("columns")
            weight_names = manifest.get("weights")
            metadata = manifest.get("catalog_metadata")
            if (
                not isinstance(column_names, list)
                or not all(isinstance(name, str) for name in column_names)
                or not isinstance(weight_names, list)
                or not all(isinstance(name, str) for name in weight_names)
                or not isinstance(metadata, dict)
            ):
                raise ValueError("Catalog archive manifest fields are invalid.")
            expected_fields = {
                "manifest_json",
                *(f"column_{index}" for index in range(len(column_names))),
                *(f"weight_{index}" for index in range(len(weight_names))),
            }
            if set(archive.files) != expected_fields:
                raise ValueError("Catalog archive field set does not match its manifest.")
            columns = {
                name: np.asarray(archive[f"column_{index}"])
                for index, name in enumerate(column_names)
            }
            weights = {
                name: np.asarray(archive[f"weight_{index}"], dtype=float)
                for index, name in enumerate(weight_names)
            }
        return cls(columns=columns, weights=weights, metadata=metadata)

    @classmethod
    def concatenate(cls, catalogs: Sequence[WeightedSubhaloCatalog]) -> WeightedSubhaloCatalog:
        """Concatenate compatible catalogs along their leading axis."""
        if not catalogs:
            raise ValueError("At least one catalog is required for concatenation.")
        first = catalogs[0]
        column_names = set(first.columns)
        weight_names = set(first.weights)
        first_metadata = (
            first.metadata.as_mapping()
            if isinstance(first.metadata, CatalogMetadata)
            else first.metadata
        )
        metadata = dict(first_metadata)
        for catalog in catalogs[1:]:
            if set(catalog.columns) != column_names or set(catalog.weights) != weight_names:
                raise ValueError("Catalog columns and weight factors must match for concatenation.")
            catalog_metadata = (
                catalog.metadata.as_mapping()
                if isinstance(catalog.metadata, CatalogMetadata)
                else catalog.metadata
            )
            if dict(catalog_metadata) != metadata:
                raise ValueError("Catalog metadata must match for concatenation.")
        return cls(
            columns={
                name: np.concatenate([catalog.columns[name] for catalog in catalogs], axis=0)
                for name in first.columns
            },
            weights={
                name: np.concatenate([catalog.weights[name] for catalog in catalogs], axis=0)
                for name in first.weights
            },
            metadata=metadata,
        )


__all__ = [
    "CANONICAL_WEIGHT_FACTORS",
    "CATALOG_SCHEMA_VERSION",
    "CatalogMetadata",
    "WeightedSubhaloCatalog",
]
