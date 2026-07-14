"""Validated state containers shared by SASHIMI-family calculations."""

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

Array = np.ndarray


def _aligned_arrays(owner: str, values: Mapping[str, Any]) -> dict[str, Array]:
    """Convert state values to arrays and require one common shape."""
    arrays = {name: np.asarray(value) for name, value in values.items()}
    shapes = {value.shape for value in arrays.values()}
    if len(shapes) > 1:
        raise ValueError(f"{owner} arrays must have one common shape; got {shapes}.")
    for name, value in arrays.items():
        if not np.all(np.isfinite(value)):
            raise ValueError(f"{owner}.{name} contains non-finite values.")
    return arrays


def _validate_weight(name: str, value: Any, shape: tuple[int, ...]) -> Array:
    """Return one finite, nonnegative weight array with the expected shape."""
    array = np.asarray(value, dtype=float)
    if array.shape != shape:
        raise ValueError(f"{name} has shape {array.shape}; expected {shape}.")
    if not np.all(np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{name} must contain finite, nonnegative values.")
    return array


@dataclass(frozen=True, slots=True)
class HostState:
    """Describe a host halo at one or more aligned redshifts."""

    redshift: Array
    time: Array
    m200: Array
    mvir: Array
    r200: Array
    rvir: Array
    concentration: Array
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate shape alignment and freeze metadata."""
        arrays = _aligned_arrays(
            type(self).__name__,
            {
                "redshift": self.redshift,
                "time": self.time,
                "m200": self.m200,
                "mvir": self.mvir,
                "r200": self.r200,
                "rvir": self.rvir,
                "concentration": self.concentration,
            },
        )
        for name, value in arrays.items():
            object.__setattr__(self, name, value)
        if np.any(arrays["redshift"] < 0.0) or np.any(arrays["time"] < 0.0):
            raise ValueError("Host redshift and cosmic time must be nonnegative.")
        for name in ("m200", "mvir", "r200", "rvir", "concentration"):
            if np.any(arrays[name] <= 0.0):
                raise ValueError(f"HostState.{name} must be positive.")
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True, slots=True)
class AccretionBatch:
    """Store accretion nodes and independent quadrature weights."""

    m200_acc: Array
    mvir_acc: Array
    z_acc: Array
    concentration_acc: Array
    weight_base: Array
    weight_concentration: Array
    metadata: Mapping[str, Any] = field(default_factory=dict)
    weight_host_history: Array | None = None
    weight_orbit: Array | None = None

    def __post_init__(self) -> None:
        """Validate node alignment, weights, and metadata."""
        arrays = _aligned_arrays(
            type(self).__name__,
            {
                "m200_acc": self.m200_acc,
                "mvir_acc": self.mvir_acc,
                "z_acc": self.z_acc,
                "concentration_acc": self.concentration_acc,
            },
        )
        shape = arrays["m200_acc"].shape
        for name, value in arrays.items():
            object.__setattr__(self, name, value)
        if np.any(arrays["z_acc"] < 0.0):
            raise ValueError("Accretion redshift must be nonnegative.")
        for name in ("m200_acc", "mvir_acc", "concentration_acc"):
            if np.any(arrays[name] <= 0.0):
                raise ValueError(f"AccretionBatch.{name} must be positive.")
        for name in ("weight_base", "weight_concentration"):
            object.__setattr__(self, name, _validate_weight(name, getattr(self, name), shape))
        for name in ("weight_host_history", "weight_orbit"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _validate_weight(name, value, shape))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True, slots=True)
class ProfileParameters:
    """Store aligned named parameters for an arbitrary halo profile."""

    name: str
    values: Mapping[str, Array]

    def __post_init__(self) -> None:
        """Require a named profile with at least one aligned parameter."""
        if not self.name.strip():
            raise ValueError("Profile name must be a non-empty string.")
        if not self.values:
            raise ValueError("ProfileParameters must contain at least one value array.")
        values = _aligned_arrays(type(self).__name__, self.values)
        object.__setattr__(self, "values", MappingProxyType(values))


@dataclass(slots=True)
class SubhaloState:
    """Describe evolved internal subhalo properties."""

    m_bound: Array
    profile: ProfileParameters
    alive: Array
    flags: Array
    extra: Mapping[str, Array] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate state arrays against the profile batch shape."""
        profile_shape = next(iter(self.profile.values.values())).shape
        reserved = {"m_bound", "alive", "flags"}.intersection(self.extra)
        if reserved:
            raise ValueError(f"SubhaloState.extra uses reserved fields: {sorted(reserved)}.")
        arrays = {
            "m_bound": np.asarray(self.m_bound),
            "alive": np.asarray(self.alive, dtype=bool),
            "flags": np.asarray(self.flags),
            **{name: np.asarray(value) for name, value in self.extra.items()},
        }
        for name, value in arrays.items():
            if value.shape != profile_shape:
                raise ValueError(f"{name} has shape {value.shape}; expected {profile_shape}.")
            if not np.all(np.isfinite(value)):
                raise ValueError(f"{name} contains non-finite values.")
        if np.any(arrays["m_bound"] < 0.0):
            raise ValueError("m_bound must be nonnegative.")
        self.m_bound = arrays.pop("m_bound")
        self.alive = arrays.pop("alive")
        self.flags = arrays.pop("flags")
        self.extra = MappingProxyType(arrays)


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

    def __post_init__(self) -> None:
        """Validate alignment and require explicit spatial semantics."""
        populated = {
            name: value
            for name in (
                "energy",
                "angular_momentum",
                "radius",
                "radial_velocity",
                "tangential_velocity",
                "pericenter",
                "apocenter",
                "phase",
            )
            if (value := getattr(self, name)) is not None
        }
        if populated and self.representation == "unspecified":
            raise ValueError("Spatial fields require an explicit orbital representation.")
        if populated:
            arrays = _aligned_arrays(type(self).__name__, populated)
            for name, value in arrays.items():
                object.__setattr__(self, name, value)


__all__ = [
    "AccretionBatch",
    "HostState",
    "OrbitalState",
    "ProfileParameters",
    "SubhaloState",
]
