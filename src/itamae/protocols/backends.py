"""Structural interfaces for cosmology and unit backends."""

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class CosmologyBackend(Protocol):
    """Provide background-cosmology quantities in canonical ITAMAE units.

    Notes
    -----
    The same background and unit convention must be used for halo masses, densities, growth
    and time. These interfaces identify and transport those choices; satisfying an interface
    does not establish physical equivalence between two providers.
    """

    @property
    def identifier(self) -> str:
        """Return a stable identifier suitable for metadata and cache keys.

        Returns
        -------
        str
            Stable prescription identity for cache separation and catalogue provenance.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def H(self, z: Any) -> Any:
        """Return the Hubble rate in km s^-1 Mpc^-1.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Hubble rate H(z), in km/s/Mpc, with the input redshift shape.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def rho_crit(self, z: Any) -> Any:
        """Return critical density in Msun Mpc^-3.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Physical critical density in Msun/Mpc^3, with the input redshift shape.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def rho_m(self, z: Any) -> Any:
        """Return physical matter density in Msun Mpc^-3.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Physical matter density in Msun/Mpc^3, with the input redshift shape.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def omega_m(self, z: Any) -> Any:
        """Return the matter density fraction.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Dimensionless matter fraction Omega_m(z).

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def growth_factor(self, z: Any) -> Any:
        """Return the linear growth factor normalized at redshift zero.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Dimensionless linear growth D(z), normalized at z=0.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def collapse_threshold(self, z: Any) -> Any:
        """Return the spherical-collapse threshold at redshift ``z``.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Dimensionless linearly extrapolated spherical-collapse barrier 1.686/D(z).

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def cosmic_time(self, z: Any) -> Any:
        """Return cosmic age in Gyr.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Cosmic age in Gyr, with the input redshift shape.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def lookback_time(self, z: Any) -> Any:
        """Return lookback time in Gyr.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Present age minus age at z, in Gyr.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...


@runtime_checkable
class UnitBackend(Protocol):
    """Convert public values to and from the canonical internal unit system.

    Notes
    -----
    The same background and unit convention must be used for halo masses, densities, growth
    and time. These interfaces identify and transport those choices; satisfying an interface
    does not establish physical equivalence between two providers.
    """

    @property
    def identifier(self) -> str:
        """Return a stable identifier including the unit-schema version.

        Returns
        -------
        str
            Stable prescription identity for cache separation and catalogue provenance.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def to_internal(self, value: Any, physical_type: str) -> Any:
        """Convert a value to a plain array in canonical units.

        Parameters
        ----------
        value : array_like
            Numerical input to be checked before use; shape and physical interpretation are set
            by the caller and retained by the validation helper.
        physical_type : str
            Quantity name in the canonical unit schema, such as mass, length, time or density.

        Returns
        -------
        numpy.ndarray
            Finite plain numerical array in the canonical unit for physical_type.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def from_internal(self, value: Any, unit: Any = None) -> Any:
        """Convert a canonical value to a requested public unit.

        Parameters
        ----------
        value : array_like
            Numerical input to be checked before use; shape and physical interpretation are set
            by the caller and retained by the validation helper.
        unit : object, optional
            Requested display unit. NativeUnits only accepts None or an empty unit; AstropyUnits
            requires a dimensionally equivalent unit.

        Returns
        -------
        numpy.ndarray or astropy.units.Quantity
            Values in the requested unit when supported; native output remains a plain canonical
            array.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...

    def validate(self, value: Any, physical_type: str) -> None:
        """Validate dimensional compatibility and finite values.

        Parameters
        ----------
        value : array_like
            Numerical input to be checked before use; shape and physical interpretation are set
            by the caller and retained by the validation helper.
        physical_type : str
            Quantity name in the canonical unit schema, such as mass, length, time or density.

        Notes
        -----
        The same background and unit convention must be used for halo masses, densities, growth
        and time. These interfaces identify and transport those choices; satisfying an interface
        does not establish physical equivalence between two providers.
        """
        ...


__all__ = ["CosmologyBackend", "UnitBackend"]
