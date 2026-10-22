import numpy as np
import pytest
from scipy.integrate import trapezoid

from itamae.spatial import (
    normalize_radial_pdf,
    orbit_radial_measure,
    radial_period,
    radial_shell_pdf,
    shell_probabilities,
    turning_points,
)


def test_radial_pdf_normalization():
    """Normalize a spatial subhalo distribution with the intended radial measure.

    Notes
    -----
    Pass criteria: Trapezoidal integration of the normalized toy radial density equals one.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    q = np.linspace(0.01, 1.0, 1000)
    pdf = normalize_radial_pdf(q, q / (q + 0.1) ** 2)
    assert trapezoid(pdf, q) == pytest.approx(1.0)


def test_shell_probabilities_normalize():
    """Keep the total probability of all radial shells equal to unity.

    Notes
    -----
    Pass criteria: The ten-shell toy measure sums to one.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    edges = np.linspace(0.0, 1.0, 11)
    centers = 0.5 * (edges[:-1] + edges[1:])
    probabilities = shell_probabilities(edges, np.ones_like(centers), centers)
    assert probabilities.sum() == pytest.approx(1.0)


def test_kepler_orbit_turning_points_period_and_measure():
    # In units GM=1, an orbit with a=1 and e=0.5 has E=-1/2,
    # L=sqrt(1-e^2), (rp, ra)=(0.5, 1.5), and period 2*pi.
    """Recover an analytic bound orbit and its residence-time distribution.

    Notes
    -----
    Pass criteria: For GM=a=1 and eccentricity 0.5, turning points are 0.5 and 1.5, period
    is 2*pi, shell weights are nonnegative and normalized, and missing endpoints are
    rejected.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    def potential(radius):
        """Evaluate the Kepler potential with GM=1 for an analytic orbital reference.

        Parameters
        ----------
        radius : float or array_like
            Positive radius or radial bound in one consistent length unit. Pericenter and
            apocenter bound the allowed radial motion.

        Returns
        -------
        numpy.ndarray
            Specific potential -1/r in consistent fixture units.

        Notes
        -----
        Recover an analytic bound orbit and its residence-time distribution. Pass criteria: For
        GM=a=1 and eccentricity 0.5, turning points are 0.5 and 1.5, period is 2*pi, shell
        weights are nonnegative and normalized, and missing endpoints are rejected. Failure
        occurs when the stated comparison or expected rejection is not met. This controlled
        calculation does not calibrate a subhalo population against simulations or observations.
        """
        return -1.0 / np.asarray(radius)

    energy = -0.5
    angular_momentum = np.sqrt(0.75)
    pericenter, apocenter = turning_points(
        potential,
        energy,
        angular_momentum,
        0.1,
        3.0,
    )
    assert pericenter == pytest.approx(0.5, rel=1.0e-10)
    assert apocenter == pytest.approx(1.5, rel=1.0e-10)
    period = radial_period(
        potential,
        energy,
        angular_momentum,
        pericenter,
        apocenter,
    )
    assert period == pytest.approx(2.0 * np.pi, rel=1.0e-9)

    measure = orbit_radial_measure(
        np.linspace(pericenter, apocenter, 17),
        potential,
        energy,
        angular_momentum,
        pericenter,
        apocenter,
    )
    assert np.sum(measure.weight) == pytest.approx(1.0)
    assert np.all(measure.weight >= 0.0)
    density = radial_shell_pdf(
        measure.radius,
        potential,
        energy,
        angular_momentum,
        pericenter,
        apocenter,
        period=period,
    )
    assert np.all(np.isfinite(density))
    assert np.all(density > 0.0)
    with pytest.raises(ValueError, match="pericenter"):
        orbit_radial_measure(
            np.linspace(pericenter + 0.01, apocenter, 17),
            potential,
            energy,
            angular_momentum,
            pericenter,
            apocenter,
        )


def test_turning_points_reject_unbound_scan():
    """Avoid assigning a bound radial period to an orbit without two turning points.

    Notes
    -----
    Pass criteria: The positive-energy radial Kepler case raises ValueError on the supplied
    scan interval.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    def potential(radius):
        """Evaluate the Kepler potential with GM=1 for an analytic orbital reference.

        Parameters
        ----------
        radius : float or array_like
            Positive radius or radial bound in one consistent length unit. Pericenter and
            apocenter bound the allowed radial motion.

        Returns
        -------
        numpy.ndarray
            Specific potential -1/r in consistent fixture units.

        Notes
        -----
        Avoid assigning a bound radial period to an orbit without two turning points. Pass
        criteria: The positive-energy radial Kepler case raises ValueError on the supplied scan
        interval. Failure occurs when the stated comparison or expected rejection is not met.
        This controlled calculation does not calibrate a subhalo population against simulations
        or observations.
        """
        return -1.0 / np.asarray(radius)

    with pytest.raises(ValueError, match="two turning points"):
        turning_points(potential, 1.0, 0.0, 0.1, 10.0)
