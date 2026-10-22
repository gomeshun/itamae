"""Power-spectrum boundaries must reject lossy casts and malformed callbacks."""

import numpy as np
import pytest

from itamae.power import TabulatedPowerSpectrum, TransferModifiedPowerSpectrum


@pytest.mark.parametrize("bad", [np.array([1 + 1j, 2 + 1j]), ["1", "2"], [True, True]])
@pytest.mark.parametrize("field", ["wavenumber", "power", "query"])
def test_tabulated_power_rejects_non_real_numeric_values(bad, field):
    """Prevent complex or nonnumeric spectrum data from being silently cast to real power.

    Parameters
    ----------
    bad : object
        Deliberately invalid numerical value or array shape used to verify explicit
        rejection.
    field : str
        Physical field or provenance field deliberately varied in this test.

    Notes
    -----
    Pass criteria: Complex arrays, strings and Boolean values are rejected for table
    coordinates, power and queries.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    kwargs = dict(wavenumber=[1.0, 2.0], power=[1.0, 2.0])
    if field == "query":
        model = TabulatedPowerSpectrum(**kwargs)
        with pytest.raises(ValueError, match="real numeric"):
            model(bad)
    else:
        kwargs[field] = bad
        with pytest.raises(ValueError, match="real numeric"):
            TabulatedPowerSpectrum(**kwargs)


@pytest.mark.parametrize(
    "bad", [np.array([1 + 1j, 2 + 1j]), ["1", "2"], [True, True], np.ones((2, 1))]
)
def test_transfer_rejects_lossy_or_misaligned_ratio(bad):
    """Keep the suppression factor aligned with the queried wavenumbers.

    Parameters
    ----------
    bad : object
        Deliberately invalid numerical value or array shape used to verify explicit
        rejection.

    Notes
    -----
    Pass criteria: Complex, text, Boolean or wrong-shape power ratios raise ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    base = TabulatedPowerSpectrum([1.0, 2.0], [1.0, 2.0])
    model = TransferModifiedPowerSpectrum(base, lambda k: bad, ratio_identifier="invalid")
    with pytest.raises(ValueError, match="real numeric|shape"):
        model(np.array([1.0, 2.0]))


def test_transfer_rejects_nonfinite_product_without_warning():
    """Report spectral overflow as an invalid physical result.

    Notes
    -----
    Pass criteria: Overflow of finite base power times a finite ratio raises ValueError
    under strict floating-point handling.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    base = TabulatedPowerSpectrum([1.0, 2.0], [1e300, 1e300])
    model = TransferModifiedPowerSpectrum(base, lambda k: 1e300, ratio_identifier="overflow")
    with np.errstate(all="raise"), pytest.raises(ValueError, match="finite"):
        model([1.0, 2.0])


def test_scalar_transfer_keeps_query_shape():
    """Apply a scale-independent power suppression over any query shape.

    Notes
    -----
    Pass criteria: A scalar factor of one half preserves a two-dimensional query and returns
    half the base spectrum exactly.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    base = TabulatedPowerSpectrum([1.0, 2.0], [1.0, 2.0])
    model = TransferModifiedPowerSpectrum(base, lambda k: 0.5, ratio_identifier="half")
    query = np.array([[1.0, 2.0], [2.0, 1.0]])
    np.testing.assert_array_equal(model(query), 0.5 * base(query))
