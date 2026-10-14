"""Validated state containers shared by SASHIMI-family calculations."""

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

Array = np.ndarray


def _aligned_arrays(owner: str, values: Mapping[str, Any]) -> dict[str, Array]:
    """Convert state values to arrays and require one common shape.

    Parameters
    ----------
    owner : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.
    values : Mapping or array_like
        Named physical fields or aligned numerical values supplied by the caller; each
        population field must refer to the same nodes.

    Returns
    -------
    dict[str, numpy.ndarray]
        Finite physical fields with one common shape.

    Notes
    -----
    All populated fields describe the same aligned host or subhalo nodes. Canonical
    dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
    meanings must be declared in their names or metadata. Containers validate
    representation, not the evolution law.
    """
    arrays = {name: np.asarray(value) for name, value in values.items()}
    shapes = {value.shape for value in arrays.values()}
    if len(shapes) > 1:
        raise ValueError(f"{owner} arrays must have one common shape; got {shapes}.")
    for name, value in arrays.items():
        if not np.all(np.isfinite(value)):
            raise ValueError(f"{owner}.{name} contains non-finite values.")
    return arrays


def _validate_weight(name: str, value: Any, shape: tuple[int, ...]) -> Array:
    """Return one finite, nonnegative weight array with the expected shape.

    Parameters
    ----------
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.
    shape : tuple of int
        Required population or coordinate shape. Extra axes cannot be silently treated as
        additional objects.

    Returns
    -------
    numpy.ndarray
        Finite nonnegative weight values with exactly the required node shape.

    Notes
    -----
    All populated fields describe the same aligned host or subhalo nodes. Canonical
    dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
    meanings must be declared in their names or metadata. Containers validate
    representation, not the evolution law.
    """
    array = np.asarray(value, dtype=float)
    if array.shape != shape:
        raise ValueError(f"{name} has shape {array.shape}; expected {shape}.")
    if not np.all(np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{name} must contain finite, nonnegative values.")
    return array


@dataclass(frozen=True, slots=True)
class HostState:
    """Describe a host halo at one or more aligned redshifts.

    Parameters
    ----------
    redshift : float or array_like
        Dimensionless redshift. The cosmological background requires z > -1; population
        nodes use nonnegative redshifts.
    time : float or array_like
        Evolution coordinate or cosmic time, in the units of the supplied equation;
        canonical state times are in Gyr.
    m200 : array_like
        Halo mass within 200 times critical density, in solar masses under the canonical
        unit schema; aligned with the population coordinates.
    mvir : array_like
        Virial halo mass in solar masses under the canonical unit schema. It is distinct
        from M200 and must be supplied with the model conversion.
    r200 : float or array_like
        Physical radius in Mpc under the canonical unit schema; use the same length
        convention as the profile.
    rvir : float or array_like
        Physical radius in Mpc under the canonical unit schema; use the same length
        convention as the profile.
    concentration : array_like
        Positive dimensionless halo concentration aligned with the accretion nodes; its
        radius definition (r200/r_s or rvir/r_s) must match the supplied mass and structural
        prescription.
    metadata : CatalogMetadata or Mapping
        Physical model, backend, unit and source identities needed to interpret or combine
        the nodes.

    Notes
    -----
    All populated fields describe the same aligned host or subhalo nodes. Canonical
    dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
    meanings must be declared in their names or metadata. Containers validate
    representation, not the evolution law.
    """

    redshift: Array
    time: Array
    m200: Array
    mvir: Array
    r200: Array
    rvir: Array
    concentration: Array
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate shape alignment and freeze metadata.

        Notes
        -----
        All populated fields describe the same aligned host or subhalo nodes. Canonical
        dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
        meanings must be declared in their names or metadata. Containers validate
        representation, not the evolution law.
        """
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
    """Store accretion nodes and independent quadrature weights.

    Parameters
    ----------
    m200_acc : array_like
        Halo mass within 200 times critical density, in solar masses under the canonical
        unit schema; aligned with the population coordinates.
    mvir_acc : array_like
        Virial halo mass in solar masses under the canonical unit schema. It is distinct
        from M200 and must be supplied with the model conversion.
    z_acc : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    concentration_acc : array_like
        Positive dimensionless halo concentration aligned with the accretion nodes; its
        radius definition (r200/r_s or rvir/r_s) must match the supplied mass and structural
        prescription.
    weight_base : array_like
        Expected population abundance from the accretion quadrature, before independent
        concentration, orbit and survival factors.
    weight_concentration : array_like
        Normalized concentration-quadrature probability, kept distinct from the accretion
        abundance.
    metadata : CatalogMetadata or Mapping
        Physical model, backend, unit and source identities needed to interpret or combine
        the nodes.
    weight_host_history : array_like, optional
        Independent nonnegative host-history or orbit quadrature factor, included once in
        the final population weight.
    weight_orbit : array_like, optional
        Independent nonnegative host-history or orbit quadrature factor, included once in
        the final population weight.

    Notes
    -----
    All populated fields describe the same aligned host or subhalo nodes. Canonical
    dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
    meanings must be declared in their names or metadata. Containers validate
    representation, not the evolution law.
    """

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
        """Validate node alignment, weights, and metadata.

        Notes
        -----
        All populated fields describe the same aligned host or subhalo nodes. Canonical
        dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
        meanings must be declared in their names or metadata. Containers validate
        representation, not the evolution law.
        """
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
    """Store aligned named parameters for an arbitrary halo profile.

    Parameters
    ----------
    name : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.
    values : Mapping or array_like
        Named physical fields or aligned numerical values supplied by the caller; each
        population field must refer to the same nodes.

    Notes
    -----
    All populated fields describe the same aligned host or subhalo nodes. Canonical
    dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
    meanings must be declared in their names or metadata. Containers validate
    representation, not the evolution law.
    """

    name: str
    values: Mapping[str, Array]

    def __post_init__(self) -> None:
        """Require a named profile with at least one aligned parameter.

        Notes
        -----
        All populated fields describe the same aligned host or subhalo nodes. Canonical
        dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
        meanings must be declared in their names or metadata. Containers validate
        representation, not the evolution law.
        """
        if not self.name.strip():
            raise ValueError("Profile name must be a non-empty string.")
        if not self.values:
            raise ValueError("ProfileParameters must contain at least one value array.")
        values = _aligned_arrays(type(self).__name__, self.values)
        object.__setattr__(self, "values", MappingProxyType(values))


@dataclass(slots=True)
class SubhaloState:
    """Describe evolved internal subhalo properties.

    Parameters
    ----------
    m_bound : array_like
        Nonnegative bound subhalo mass, in solar masses under the canonical state
        convention.
    profile : ProfileParameters
        Aligned named internal halo-profile parameters and their declared physical meaning.
    alive : array_like of bool
        Survival state aligned with each subhalo profile; a false value denotes disruption
        in the model convention.
    flags : array_like
        Finite model-defined diagnostic flags aligned with the evolved subhalo population.
    extra : Mapping, optional
        Additional named model/provenance information; reserved standard fields cannot be
        overridden.

    Notes
    -----
    All populated fields describe the same aligned host or subhalo nodes. Canonical
    dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
    meanings must be declared in their names or metadata. Containers validate
    representation, not the evolution law.
    """

    m_bound: Array
    profile: ProfileParameters
    alive: Array
    flags: Array
    extra: Mapping[str, Array] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate state arrays against the profile batch shape.

        Notes
        -----
        All populated fields describe the same aligned host or subhalo nodes. Canonical
        dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
        meanings must be declared in their names or metadata. Containers validate
        representation, not the evolution law.
        """
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
    """Represent optional instantaneous or orbit-averaged spatial information.

    Parameters
    ----------
    energy : float or array_like
        Specific orbital energy, in the same velocity-squared units as the supplied
        potential.
    angular_momentum : float or array_like
        Nonnegative specific angular momentum, in length times velocity units consistent
        with the potential.
    radius : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    radial_velocity : array_like, optional
        Orbital velocity component in the common length/time convention of the state.
    tangential_velocity : array_like, optional
        Orbital velocity component in the common length/time convention of the state.
    pericenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    apocenter : float or array_like
        Positive radius or radial bound in one consistent length unit. Pericenter and
        apocenter bound the allowed radial motion.
    phase : array_like, optional
        Orbital phase in the convention named by representation; no orbit model is assumed
        by this container.
    representation : str
        Explicit meaning of the spatial fields or radial weights, such as shell-probability;
        not an implicit volume-density convention.

    Notes
    -----
    All populated fields describe the same aligned host or subhalo nodes. Canonical
    dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
    meanings must be declared in their names or metadata. Containers validate
    representation, not the evolution law.
    """

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
        """Validate alignment and require explicit spatial semantics.

        Notes
        -----
        All populated fields describe the same aligned host or subhalo nodes. Canonical
        dimensional fields use Msun, Mpc, km/s and Gyr; model-specific profile and orbit
        meanings must be declared in their names or metadata. Containers validate
        representation, not the evolution law.
        """
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
