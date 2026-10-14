"""Colossus cosmology adapter."""

from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
from colossus.cosmology import cosmology as col_cosmology


@dataclass(frozen=True, slots=True)
class ColossusCosmology:
    r"""Wrap a named Colossus cosmology while restoring global state after calls.

    Parameters
    ----------
    name : str
        Named Colossus cosmology, such as planck18, used for every background query.

    Notes
    -----
    The default planck18 preset is supplied by Colossus and is based on Planck 2018. The
    preset is external library data; its full parameter set is not inferred from the rounded
    native defaults.

    Each query activates the named Colossus cosmology and restores the previous global
    selection. The wrapper preserves the implemented provider conversion; comparisons must
    also check factors of h, not only matching parameter labels.

    References
    ----------

    .. [cite-colossus] Diemer (2018), COLOSSUS: A Python Toolkit for Cosmology, Large-scale Structure, and Dark Matter Halos.
        https://arxiv.org/abs/1712.04512

    .. [cite-planck] Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Table 2.
        https://arxiv.org/abs/1807.06209
    """

    name: str = "planck18"
    _params: dict = field(init=False, repr=False)

    def __post_init__(self):
        r"""Snapshot the named cosmology parameters for this immutable provider.

        Notes
        -----
        The default planck18 preset is supplied by Colossus and is based on Planck 2018. The
        preset is external library data; its full parameter set is not inferred from the rounded
        native defaults.

        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.

        References
        ----------

        .. [cite-colossus] Diemer (2018), COLOSSUS: A Python Toolkit for Cosmology, Large-scale Structure, and Dark Matter Halos.
            https://arxiv.org/abs/1712.04512

        .. [cite-planck] Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Table 2.
            https://arxiv.org/abs/1807.06209
        """
        object.__setattr__(self, "_params", dict(col_cosmology.cosmologies[self.name]))

    @property
    def identifier(self) -> str:
        """Return the backend identifier.

        Returns
        -------
        str
            Stable prescription identity for cache separation and catalogue provenance.

        Notes
        -----
        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.
        """
        return f"colossus:{self.name}"

    @contextmanager
    def _active(self):
        """Temporarily select the requested cosmology for a background query.

        Yields
        ------
        colossus.cosmology.cosmology.Cosmology
            Active provider; the previously selected cosmology is restored on exit.

        Notes
        -----
        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.
        """
        try:
            previous = col_cosmology.getCurrent()
        except Exception:
            previous = None
        previous_name = previous.name if previous is not None else None
        current = col_cosmology.setCosmology(self.name)
        try:
            yield current
        finally:
            if previous_name is not None:
                col_cosmology.setCosmology(previous_name)
            else:
                col_cosmology.current_cosmo = None

    def H(self, z):
        """Return Hubble rate in km s^-1 Mpc^-1.

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
        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.
        """
        with self._active() as cosmo:
            return np.asarray(cosmo.Hz(z))

    def rho_crit(self, z):
        """Return the Colossus critical density with kpc converted to Mpc.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        numpy.ndarray
            Density in Msun*h^2/Mpc^3 under the implemented conversion, retaining redshift
            shape.

        Notes
        -----
        Colossus supplies density in Msun*h^2/kpc^3. This implementation multiplies by 1e9 only
        and therefore retains the h^2 convention. It does not yet satisfy the canonical backend
        density contract in Msun/Mpc^3; multiplying these numerical values by the selected
        cosmology h**2 is needed for that conversion. The H/growth comparison test does not
        cover this density conversion.
        """
        with self._active() as cosmo:
            return np.asarray(cosmo.rho_c(z)) * 1.0e9

    def rho_m(self, z):
        """Return the Colossus matter density with kpc converted to Mpc.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        numpy.ndarray
            Density in Msun*h^2/Mpc^3 under the implemented conversion, retaining redshift
            shape.

        Notes
        -----
        Colossus supplies density in Msun*h^2/kpc^3. This implementation multiplies by 1e9 only
        and therefore retains the h^2 convention. It does not yet satisfy the canonical backend
        density contract in Msun/Mpc^3; multiplying these numerical values by the selected
        cosmology h**2 is needed for that conversion. The H/growth comparison test does not
        cover this density conversion.
        """
        with self._active() as cosmo:
            return np.asarray(cosmo.rho_m(z)) * 1.0e9

    def omega_m(self, z):
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
        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.
        """
        with self._active() as cosmo:
            return np.asarray(cosmo.Om(z))

    def growth_factor(self, z):
        """Return normalized linear growth factor.

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
        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.
        """
        with self._active() as cosmo:
            return np.asarray(cosmo.growthFactor(z))

    def collapse_threshold(self, z):
        r"""Return 1.686 divided by the growth factor.

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
        The adopted spherical-collapse barrier is delta_c(z)=1.686/D(z), as used in Ludlow et
        al., Appendix C. This approximation is separate from any mass-dependent WDM collapse
        correction.

        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.

        References
        ----------

        .. [cite-ludlow] Ludlow et al. (2016), The mass-concentration-redshift relation of cold and warm dark matter haloes, Appendix C.
            https://arxiv.org/abs/1601.02624
        """
        return 1.686 / self.growth_factor(z)

    def cosmic_time(self, z):
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
        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.
        """
        with self._active() as cosmo:
            return np.asarray(cosmo.age(z))

    def lookback_time(self, z):
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
        Each query activates the named Colossus cosmology and restores the previous global
        selection. The wrapper preserves the implemented provider conversion; comparisons must
        also check factors of h, not only matching parameter labels.
        """
        return self.cosmic_time(0.0) - self.cosmic_time(z)
