"""Quadrature rules used by SASHIMI-family population integrals."""

import numpy as np
from numpy.polynomial.hermite import hermgauss


def gauss_hermite_lognormal(median, sigma_log10: float, order: int = 5):
    """Return nodes and normalized weights for log-normal concentration scatter.

    Parameters
    ----------
    median
        Median of the positive log-normal variable.
    sigma_log10
        Standard deviation in base-10 logarithm.
    order
        Gauss-Hermite quadrature order.

    Returns
    -------
    nodes, weights
        Arrays with the quadrature axis prepended to the input shape.
    """
    if order < 1:
        raise ValueError("Quadrature order must be positive.")
    median = np.asarray(median, dtype=float)
    if np.any(median <= 0.0):
        raise ValueError("Median must be positive.")
    x, w = hermgauss(order)
    shape = (order,) + (1,) * median.ndim
    nodes = 10.0 ** (np.log10(median)[None, ...] + np.sqrt(2.0) * sigma_log10 * x.reshape(shape))
    weights = np.broadcast_to((w / np.sqrt(np.pi)).reshape(shape), nodes.shape)
    return nodes, weights
