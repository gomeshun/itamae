"""Normalized radial distributions for deterministic spatial quadrature."""

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import trapezoid


@dataclass(frozen=True, slots=True)
class RadialMeasure:
    """Store aligned radius nodes and normalized nonnegative probabilities."""

    radius: np.ndarray
    weight: np.ndarray
    representation: str = "shell-probability"

    def __post_init__(self) -> None:
        """Validate shape, domain, normalization, and representation."""
        radius = np.asarray(self.radius, dtype=float)
        weight = np.asarray(self.weight, dtype=float)
        if (
            radius.ndim != 1
            or radius.size == 0
            or weight.shape != radius.shape
            or not np.all(np.isfinite(radius))
            or np.any(radius < 0.0)
            or np.any(np.diff(radius) <= 0.0)
        ):
            raise ValueError("Radius nodes must be aligned, finite, and strictly increasing.")
        if (
            not np.all(np.isfinite(weight))
            or np.any(weight < 0.0)
            or not np.isclose(np.sum(weight), 1.0, rtol=0.0, atol=1.0e-12)
        ):
            raise ValueError("Radial weights must be finite, nonnegative, and normalized.")
        if not isinstance(self.representation, str) or not self.representation.strip():
            raise ValueError("Radial representation must be a non-empty string.")
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "weight", weight)


def normalize_radial_pdf(q, pdf_q):
    """Normalize a probability density defined per unit ``q``."""
    q = np.asarray(q, dtype=float)
    pdf_q = np.asarray(pdf_q, dtype=float)
    if q.ndim != 1 or pdf_q.shape != q.shape:
        raise ValueError("q and pdf_q must be aligned one-dimensional arrays.")
    if np.any(np.diff(q) <= 0.0) or np.any(pdf_q < 0.0):
        raise ValueError("q must increase and pdf_q must be nonnegative.")
    norm = trapezoid(pdf_q, q)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("Radial density has a nonpositive or nonfinite integral.")
    return pdf_q / norm


def shell_probabilities(q_edges, pdf_q, q_centers):
    """Approximate shell probabilities from a density sampled at shell centers."""
    q_edges = np.asarray(q_edges, dtype=float)
    q_centers = np.asarray(q_centers, dtype=float)
    pdf_q = np.asarray(pdf_q, dtype=float)
    if q_edges.size != q_centers.size + 1 or pdf_q.shape != q_centers.shape:
        raise ValueError("Edges must have one more element than centers and density values.")
    probabilities = pdf_q * np.diff(q_edges)
    total = probabilities.sum()
    if total <= 0.0:
        raise ValueError("Shell probabilities have nonpositive total weight.")
    return probabilities / total


def radial_measure(radius: Any, weight: Any) -> RadialMeasure:
    """Construct a normalized shell-probability measure.

    Input weights may have any positive finite normalization; the returned
    object always sums to one.
    """
    values = np.asarray(weight, dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("Radial weights must be finite and nonnegative.")
    total = float(np.sum(values))
    if total <= 0.0:
        raise ValueError("Radial weights must have positive total measure.")
    return RadialMeasure(radius=np.asarray(radius, dtype=float), weight=values / total)


__all__ = [
    "RadialMeasure",
    "normalize_radial_pdf",
    "radial_measure",
    "shell_probabilities",
]
