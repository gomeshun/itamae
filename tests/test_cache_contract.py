"""Cache identities must follow actual inputs, even with human-readable labels."""

from dataclasses import replace

import numpy as np
import pytest

from itamae.power import SharpKWindow, TabulatedPowerSpectrum
from itamae.variance import (
    IntegratedVarianceModel,
    load_variance_cache,
    save_variance_cache,
    variance_cache_key,
)


def test_spectrum_label_cannot_alias_different_table_content():
    """Keep density variances tied to the actual input power spectrum.

    Notes
    -----
    Pass criteria: Changing a power sample changes the identifier despite an equal human
    label; numerically equal tables with different byte order retain one identifier.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    first = TabulatedPowerSpectrum([1.0, 2.0], [3.0, 4.0], identifier="same-source")
    second = TabulatedPowerSpectrum([1.0, 2.0], [3.0, 5.0], identifier="same-source")
    assert first.identifier != second.identifier
    copied = TabulatedPowerSpectrum(np.array([1, 2], dtype=">f8"), [3, 4], identifier="same-source")
    assert first.identifier == copied.identifier


def test_spectrum_identity_includes_extrapolation_policy():
    """Distinguish spectra that make different assumptions outside their tabulated range.

    Notes
    -----
    Pass criteria: The same power table with extrapolation enabled and disabled has
    different identities.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    first = TabulatedPowerSpectrum([1.0, 2.0], [3.0, 4.0])
    second = TabulatedPowerSpectrum([1.0, 2.0], [3.0, 4.0], extrapolate=True)
    assert first.identifier != second.identifier


def test_variance_identity_includes_derivative_resolution_but_not_batching():
    """Distinguish numerical controls that can change the EPS variance derivative.

    Notes
    -----
    Pass criteria: Derivative step and sharp-k quadrature order change the identity; chunk
    size does not.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    model = IntegratedVarianceModel(
        power=TabulatedPowerSpectrum([1.0, 2.0], [3.0, 4.0]),
        window=SharpKWindow(),
        rho_mean=1.0,
        k_min=1.0,
        k_max=2.0,
    )
    assert model.identifier != replace(model, derivative_step=2.0e-4).identifier
    assert model.identifier == replace(model, chunk_size=3).identifier
    assert model.identifier != replace(model, sharp_k_order=4).identifier


@pytest.mark.parametrize("field", ["mass", "variance"])
def test_cache_detects_finite_payload_corruption(tmp_path, field):
    """Detect an altered mass grid or variance table even when all entries remain physical.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated directory for the test catalogue, cache or provenance files.
    field : str
        Physical field or provenance field deliberately varied in this test.

    Notes
    -----
    Pass criteria: Changing one finite cached mass or variance value by ten percent causes
    the stored digest check to fail.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    mass = np.array([1.0, 2.0, 4.0])
    key = variance_cache_key("model", mass, backend_identifier="backend")
    path = tmp_path / "variance.npz"
    save_variance_cache(path, key=key, mass=mass, variance=mass**-1)
    with np.load(path, allow_pickle=False) as cache:
        payload = {name: cache[name] for name in cache.files}
    payload[field][1] *= 1.1  # Still finite, nonnegative and in mass order.
    np.savez(path, **payload)
    with pytest.raises(ValueError, match="digest|corrupt"):
        load_variance_cache(path, expected_key=key)


@pytest.mark.parametrize("mass", [[], [1.0]])
def test_cache_rejects_undersized_grids(tmp_path, mass):
    """Reject a variance cache that cannot define a usable tabulated mass dependence.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated directory for the test catalogue, cache or provenance files.
    mass : object
        Halo mass or deliberately invalid mass input for this controlled domain case.

    Notes
    -----
    Pass criteria: Empty and one-node mass grids are rejected with ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="at least two"):
        save_variance_cache(tmp_path / "variance.npz", key="0" * 64, mass=mass, variance=mass)


def test_cache_rejects_old_unverified_schema(tmp_path):
    """Prevent an unverifiable historical cache from being reused as trusted variance data.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated directory for the test catalogue, cache or provenance files.

    Notes
    -----
    Pass criteria: A version-1 archive without the current verified fields is rejected.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    path = tmp_path / "old.npz"
    np.savez(path, schema_version="1.0", key="0" * 64, mass=[1.0, 2.0], variance=[2.0, 1.0])
    with pytest.raises(ValueError, match="schema|field"):
        load_variance_cache(path, expected_key="0" * 64)
