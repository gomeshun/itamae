"""Native flat-LCDM cosmology implementation."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

_G_MPC_KMS2_MSUN = 4.30091e-9
_KMS_MPC_TO_GYR_INV = 1.0227121650537077e-3


def _redshift(value):
    """Require a real expansion epoch with positive scale factor.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.

    Returns
    -------
    numpy.ndarray
        Finite real redshifts greater than -1, retaining input shape.

    Notes
    -----
    This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
    Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
    calibration of a SASHIMI halo prescription.
    """
    array = np.asarray(value)
    if array.dtype.kind not in "biuf":
        raise ValueError("Redshift must be real, finite and greater than -1.")
    array = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array <= -1.0):
        raise ValueError("Redshift must be real, finite and greater than -1.")
    return array


def _finite_background(evaluate, quantity):
    """Evaluate a cosmological quantity only within its finite numerical domain.

    Parameters
    ----------
    evaluate : callable
        Zero-argument numerical evaluation whose overflow or non-finite result must be
        reported explicitly.
    quantity : str
        Name of the quantity or component included in an error or provenance record so an
        invalid scientific input can be located.

    Returns
    -------
    float or numpy.ndarray
        Finite background quantity returned by evaluate in its documented units.

    Notes
    -----
    This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
    Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
    calibration of a SASHIMI halo prescription.
    """
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            result = evaluate()
    except FloatingPointError as error:
        raise ValueError(f"{quantity} exceeds the finite numerical background domain.") from error
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{quantity} must be finite.")
    return result


@dataclass(frozen=True, slots=True)
class NativeFlatLCDM:
    r"""Flat-LCDM background with NumPy/SciPy only.

    Parameters
    ----------
    omega_m0
        Present-day matter density fraction.
    h
        Reduced Hubble constant.

    Notes
    -----
    The defaults Omega_m0=0.315 and h=0.674 are rounded Planck-2018 values. User-supplied
    alternatives remain explicit background choices.

    This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
    Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
    calibration of a SASHIMI halo prescription.

    References
    ----------

    .. [cite-planck] Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Table 2.
        https://arxiv.org/abs/1807.06209
    """

    omega_m0: float = 0.315
    h: float = 0.674

    def __post_init__(self):
        r"""Validate the matter-plus-nonnegative-Lambda background domain.

        Notes
        -----
        The defaults Omega_m0=0.315 and h=0.674 are rounded Planck-2018 values. User-supplied
        alternatives remain explicit background choices.

        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.

        References
        ----------

        .. [cite-planck] Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Table 2.
            https://arxiv.org/abs/1807.06209
        """
        for name in ("omega_m0", "h"):
            value = np.asarray(getattr(self, name))
            if value.ndim != 0 or value.dtype.kind not in "iuf" or not np.isfinite(value):
                raise ValueError(f"{name} must be a finite real scalar.")
            object.__setattr__(self, name, float(value))
        if not 0.0 < self.omega_m0 <= 1.0:
            raise ValueError("omega_m0 must satisfy 0 < omega_m0 <= 1.")
        if self.h <= 0.0:
            raise ValueError("h must be positive.")

    @property
    def identifier(self) -> str:
        """Return a stable identifier including cosmological parameters.

        Returns
        -------
        str
            Stable prescription identity for cache separation and catalogue provenance.

        Notes
        -----
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        return f"native-flatlcdm:Om={self.omega_m0}:h={self.h}"

    @property
    def omega_lambda0(self) -> float:
        """Return the present-day dark-energy density fraction.

        Returns
        -------
        float
            Present dark-energy fraction 1-Omega_m0.

        Notes
        -----
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        return 1.0 - self.omega_m0

    def e2(self, z):
        """Return the dimensionless squared expansion rate ``E(z)^2``.

        Parameters
        ----------
        z : float or array_like
            Dimensionless redshift. The cosmological background requires z > -1; population
            nodes use nonnegative redshifts.

        Returns
        -------
        float or numpy.ndarray
            Dimensionless squared expansion rate E(z)**2=Omega_m0*(1+z)**3+Omega_Lambda0.

        Notes
        -----
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        z = _redshift(z)
        try:
            with np.errstate(over="raise", invalid="raise"):
                result = self.omega_m0 * (1.0 + z) ** 3 + self.omega_lambda0
        except FloatingPointError as error:
            raise ValueError(
                "Redshift exceeds the finite numerical expansion-rate domain."
            ) from error
        return result

    def H(self, z):
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
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        return _finite_background(lambda: 100.0 * self.h * np.sqrt(self.e2(z)), "Hubble rate")

    def rho_crit(self, z):
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
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        return _finite_background(
            lambda: 3.0 * self.H(z) ** 2 / (8.0 * np.pi * _G_MPC_KMS2_MSUN),
            "Critical density",
        )

    def rho_m(self, z):
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
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        return _finite_background(
            lambda: self.omega_m0 * self.rho_crit(0.0) * (1.0 + _redshift(z)) ** 3,
            "Matter density",
        )

    def omega_m(self, z):
        """Return the redshift-dependent matter density fraction.

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
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        z = _redshift(z)
        return self.omega_m0 * (1.0 + z) ** 3 / self.e2(z)

    def growth_factor(self, z):
        r"""Return the Carroll-Press-Turner growth approximation normalized at z=0.

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
        The factors 5/2, 4/7 and 1/70 belong to the Carroll--Press--Turner growth approximation.
        Redshift derivatives differentiate that approximation; they are not a separate growth
        calibration.

        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.

        References
        ----------

        .. [cite-growth] Carroll, Press & Turner (1992), The Cosmological Constant, Section 3.5.
            https://doi.org/10.1146/annurev.aa.30.090192.002435
        """
        z = _redshift(z)
        om = self.omega_m(z)
        ol = 1.0 - om
        g = 2.5 * om / (om ** (4.0 / 7.0) - ol + (1.0 + om / 2.0) * (1.0 + ol / 70.0))
        om0 = self.omega_m0
        ol0 = self.omega_lambda0
        g0 = 2.5 * om0 / (om0 ** (4.0 / 7.0) - ol0 + (1.0 + om0 / 2.0) * (1.0 + ol0 / 70.0))
        return g / (g0 * (1.0 + z))

    def collapse_threshold(self, z):
        r"""Return the spherical-collapse threshold scaled by the growth factor.

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

        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.

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
        The integral is evaluated independently for scalar inputs. Array inputs
        are vectorized to preserve a simple public API.

        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """

        def one(zi: float) -> float:
            """Integrate cosmic age from one redshift to the infinite-redshift limit.

            Parameters
            ----------
            zi : float or array_like
                Dimensionless redshift. The cosmological background requires z > -1; population
                nodes use nonnegative redshifts.

            Returns
            -------
            float
                Cosmic age in Gyr from quadrature of |dt/dz|.

            Notes
            -----
            This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
            Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
            calibration of a SASHIMI halo prescription.
            """

            def integrand(zp: float) -> float:
                """Evaluate the positive lookback-time element |dt/dz|.

                Parameters
                ----------
                zp : float or array_like
                    Dimensionless redshift. The cosmological background requires z > -1; population
                    nodes use nonnegative redshifts.

                Returns
                -------
                float
                    Gyr per unit redshift, 1/[(1+z) H(z)] after Hubble-unit conversion.

                Notes
                -----
                This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
                Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
                calibration of a SASHIMI halo prescription.
                """
                return 1.0 / ((1.0 + zp) * self.H(zp) * _KMS_MPC_TO_GYR_INV)

            return quad(integrand, zi, np.inf, epsabs=1e-9, epsrel=1e-9)[0]

        arr = _redshift(z)
        out = np.vectorize(one, otypes=[float])(arr)
        return float(out) if out.ndim == 0 else out

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
        This matter-plus-Lambda, spatially flat background neglects radiation. Growth uses the
        Carroll-Press-Turner approximation normalized to D(0)=1; the generic background is not a
        calibration of a SASHIMI halo prescription.
        """
        return self.cosmic_time(0.0) - self.cosmic_time(z)
