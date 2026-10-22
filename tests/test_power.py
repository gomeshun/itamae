import numpy as np
import pytest

from itamae.power import (
    SharpKWindow,
    SphericalTopHatWindow,
    TabulatedPowerSpectrum,
    TransferModifiedPowerSpectrum,
)
from itamae.protocols import PowerSpectrum, WindowFunction


def test_tabulated_power_spectrum_log_interpolation_and_provenance():
    """Recover a power-law spectrum without silently extending its support.

    Notes
    -----
    Pass criteria: A 7*k**-2 table reproduces queried powers, exposes its source label and
    protocol, and rejects a point below k_min.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    k = np.logspace(-3, 3, 61)
    spectrum = TabulatedPowerSpectrum(k, 7.0 * k**-2, identifier="power-law")

    assert isinstance(spectrum, PowerSpectrum)
    assert "power-law" in spectrum.identifier
    np.testing.assert_allclose(
        spectrum(np.array([0.01, 1.0, 100.0])), 7.0 * np.array([0.01, 1.0, 100.0]) ** -2
    )
    with pytest.raises(ValueError, match="outside"):
        spectrum(1.0e-4)


def test_transfer_modified_spectrum_takes_explicit_power_ratio():
    """Apply a power suppression once, without confusing transfer amplitude with its square.

    Notes
    -----
    Pass criteria: Multiplying a unit base spectrum by exp(-k**2) yields that power ratio
    and retains its identifier.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    k = np.logspace(-2, 2, 41)
    base = TabulatedPowerSpectrum(k, np.ones_like(k), identifier="base")
    modified = TransferModifiedPowerSpectrum(
        base,
        lambda value: np.exp(-(np.asarray(value) ** 2)),
        ratio_identifier="test-cutoff",
    )

    np.testing.assert_allclose(modified(k), np.exp(-(k**2)))
    assert "test-cutoff" in modified.identifier


def test_linear_power_interpolation_preserves_physical_zero():
    """Represent a genuine spectral node without imposing a positive floor.

    Notes
    -----
    Pass criteria: Linear interpolation preserves a zero tabulated power; log-log
    interpolation rejects that zero.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    spectrum = TabulatedPowerSpectrum(
        [1.0, 2.0, 3.0],
        [1.0, 0.0, 1.0],
        interpolation="linear",
        identifier="oscillatory-zero",
    )
    assert spectrum(2.0) == pytest.approx(0.0)
    with pytest.raises(ValueError, match="strictly positive"):
        TabulatedPowerSpectrum([1.0, 2.0], [1.0, 0.0])


def test_window_limits_and_protocols():
    """Use correctly normalized smoothing windows at long wavelength.

    Notes
    -----
    Pass criteria: Top-hat tends to one at x=0 and 1e-6; sharp-k gives [1,1,1,0] at
    x=[0,1e-6,1,2].

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    top_hat = SphericalTopHatWindow()
    sharp_k = SharpKWindow()
    assert isinstance(top_hat, WindowFunction)
    assert isinstance(sharp_k, WindowFunction)

    x = np.array([0.0, 1.0e-6, 1.0, 2.0])
    assert top_hat(x[:2]) == pytest.approx([1.0, 1.0])
    np.testing.assert_array_equal(sharp_k(x), [1.0, 1.0, 1.0, 0.0])


@pytest.mark.parametrize(
    ("wavenumber", "power"),
    [
        ([1.0, 1.0], [1.0, 2.0]),
        ([1.0, 2.0], [1.0, 0.0]),
        ([1.0, np.nan], [1.0, 2.0]),
    ],
)
def test_tabulated_power_rejects_invalid_tables(wavenumber, power):
    """Reject a spectrum that cannot support a well-defined log interpolation.

    Parameters
    ----------
    wavenumber : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.
    power : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Duplicate k, nonpositive power for log interpolation, and NaN k each
    raise ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError):
        TabulatedPowerSpectrum(wavenumber, power)
