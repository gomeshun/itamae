import numpy as np
import pytest
from scipy.integrate import trapezoid

from itamae.spatial import normalize_radial_pdf, shell_probabilities


def test_radial_pdf_normalization():
    q = np.linspace(0.01, 1.0, 1000)
    pdf = normalize_radial_pdf(q, q / (q + 0.1) ** 2)
    assert trapezoid(pdf, q) == pytest.approx(1.0)


def test_shell_probabilities_normalize():
    edges = np.linspace(0.0, 1.0, 11)
    centers = 0.5 * (edges[:-1] + edges[1:])
    probabilities = shell_probabilities(edges, np.ones_like(centers), centers)
    assert probabilities.sum() == pytest.approx(1.0)
