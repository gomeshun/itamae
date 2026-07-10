"""Weighted subhalo catalog representation."""

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True, slots=True)
class WeightedSubhaloCatalog:
    """Store columns and factorized statistical weights.

    Parameters
    ----------
    columns
        Mapping of column names to arrays sharing one leading shape.
    weights
        Mapping of weight-factor names to arrays with the same shape.
    metadata
        Reproducibility metadata, including backend and model identifiers.
    """

    columns: Mapping[str, np.ndarray]
    weights: Mapping[str, np.ndarray]
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        shapes = {np.asarray(v).shape for v in [*self.columns.values(), *self.weights.values()]}
        if len(shapes) > 1:
            raise ValueError(f"All catalog arrays must have the same shape; got {shapes}.")

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the common catalog shape."""
        source = next(iter(self.columns.values()), next(iter(self.weights.values()), np.array([])))
        return np.asarray(source).shape

    @property
    def weight_final(self) -> np.ndarray:
        """Multiply all independent weight factors."""
        result = np.ones(self.shape, dtype=float)
        for value in self.weights.values():
            result *= np.asarray(value, dtype=float)
        return result

    def select(self, mask) -> "WeightedSubhaloCatalog":
        """Return a catalog subset while preserving metadata."""
        return WeightedSubhaloCatalog(
            columns={k: np.asarray(v)[mask] for k, v in self.columns.items()},
            weights={k: np.asarray(v)[mask] for k, v in self.weights.items()},
            metadata=self.metadata,
        )

    def weighted_sum(self, values) -> float:
        """Return the weighted sum of an array aligned with the catalog."""
        values = np.asarray(values, dtype=float)
        if values.shape != self.shape:
            raise ValueError("Values must have the catalog shape.")
        return float(np.sum(values * self.weight_final))
