"""State containers shared by SASHIMI-family calculations."""

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

Array = np.ndarray


@dataclass(frozen=True, slots=True)
class HostState:
    """Describe a host halo at one or more redshifts."""

    redshift: Array
    time: Array
    m200: Array
    mvir: Array
    r200: Array
    rvir: Array
    concentration: Array
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class AccretionBatch:
    """Store accretion nodes and their independent quadrature weights."""

    m200_acc: Array
    mvir_acc: Array
    z_acc: Array
    concentration_acc: Array
    weight_base: Array
    weight_concentration: Array
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class ProfileParameters:
    """Store named profile parameters for arbitrary halo profiles."""

    name: str
    values: Mapping[str, Array]


@dataclass(slots=True)
class SubhaloState:
    """Describe evolved internal subhalo properties."""

    m_bound: Array
    profile: ProfileParameters
    alive: Array
    flags: Array
    extra: Mapping[str, Array] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class OrbitalState:
    """Represent optional instantaneous or orbit-averaged spatial information."""

    energy: Array | None = None
    angular_momentum: Array | None = None
    radius: Array | None = None
    radial_velocity: Array | None = None
    tangential_velocity: Array | None = None
    pericenter: Array | None = None
    apocenter: Array | None = None
    phase: Array | None = None
    representation: str = "unspecified"
