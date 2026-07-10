"""Colossus cosmology adapter."""

from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
from colossus.cosmology import cosmology as col_cosmology


@dataclass(frozen=True, slots=True)
class ColossusCosmology:
    """Wrap a named Colossus cosmology while restoring global state after calls."""

    name: str = "planck18"
    _params: dict = field(init=False, repr=False)

    def __post_init__(self):
        object.__setattr__(self, "_params", dict(col_cosmology.cosmologies[self.name]))

    @property
    def identifier(self) -> str:
        """Return the backend identifier."""
        return f"colossus:{self.name}"

    @contextmanager
    def _active(self):
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
        """Return Hubble rate in km s^-1 Mpc^-1."""
        with self._active() as cosmo:
            return np.asarray(cosmo.Hz(z))

    def rho_crit(self, z):
        """Return critical density in Msun Mpc^-3."""
        with self._active() as cosmo:
            return np.asarray(cosmo.rho_c(z)) * 1.0e9

    def rho_m(self, z):
        """Return matter density in Msun Mpc^-3."""
        with self._active() as cosmo:
            return np.asarray(cosmo.rho_m(z)) * 1.0e9

    def omega_m(self, z):
        """Return the matter density fraction."""
        with self._active() as cosmo:
            return np.asarray(cosmo.Om(z))

    def growth_factor(self, z):
        """Return normalized linear growth factor."""
        with self._active() as cosmo:
            return np.asarray(cosmo.growthFactor(z))

    def collapse_threshold(self, z):
        """Return 1.686 divided by the growth factor."""
        return 1.686 / self.growth_factor(z)

    def cosmic_time(self, z):
        """Return cosmic age in Gyr."""
        with self._active() as cosmo:
            return np.asarray(cosmo.age(z))

    def lookback_time(self, z):
        """Return lookback time in Gyr."""
        return self.cosmic_time(0.0) - self.cosmic_time(z)
