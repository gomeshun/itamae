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
    with pytest.raises(ValueError):
        TabulatedPowerSpectrum(wavenumber, power)
