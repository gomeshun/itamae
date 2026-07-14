"""Builders for aligned accretion nodes and factorized weights."""

from typing import Any, Mapping

import numpy as np

from itamae.types import AccretionBatch


def build_accretion_batch(
    m200_acc: Any,
    z_acc: Any,
    concentration_acc: Any,
    weight_base: Any,
    weight_concentration: Any,
    *,
    mvir_acc: Any,
    weight_host_history: Any | None = None,
    weight_orbit: Any | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> AccretionBatch:
    """Broadcast population nodes into a validated accretion batch.

    Parameters
    ----------
    m200_acc, z_acc, concentration_acc
        Accretion mass, redshift, and concentration nodes.
    weight_base, weight_concentration
        Independent population and concentration-quadrature weights.
    mvir_acc
        Virial accretion mass. It is required explicitly because ITAMAE does
        not assume equality between mass definitions.
    weight_host_history, weight_orbit
        Optional independent weight factors.
    metadata
        Model and backend provenance attached to the batch.

    Returns
    -------
    AccretionBatch
        Flattened, aligned arrays with validated nonnegative weights.
    """
    values = [m200_acc, z_acc, concentration_acc, weight_base, weight_concentration]
    values.append(mvir_acc)
    if weight_host_history is not None:
        values.append(weight_host_history)
    if weight_orbit is not None:
        values.append(weight_orbit)
    broadcast = np.broadcast_arrays(*values)
    iterator = iter(broadcast)
    m200 = next(iterator).reshape(-1)
    redshift = next(iterator).reshape(-1)
    concentration = next(iterator).reshape(-1)
    base = next(iterator).reshape(-1)
    concentration_weight = next(iterator).reshape(-1)
    mvir = next(iterator).reshape(-1)
    host_weight = None if weight_host_history is None else next(iterator).reshape(-1)
    orbit_weight = None if weight_orbit is None else next(iterator).reshape(-1)
    return AccretionBatch(
        m200_acc=m200,
        mvir_acc=mvir,
        z_acc=redshift,
        concentration_acc=concentration,
        weight_base=base,
        weight_concentration=concentration_weight,
        weight_host_history=host_weight,
        weight_orbit=orbit_weight,
        metadata={} if metadata is None else metadata,
    )
