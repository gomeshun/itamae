"""Normalized radial distributions for deterministic spatial quadrature."""

import numpy as np
from scipy.integrate import trapezoid


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
