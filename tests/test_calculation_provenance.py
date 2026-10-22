"""New calculation identities coexist with immutable historical provenance."""

import numpy as np
import pytest

from itamae.provenance import (
    CALCULATION_METADATA_KEYS,
    build_calculation_metadata,
    build_migration_metadata,
)
from itamae.types import WeightedSubhaloCatalog


@pytest.fixture
def identifiers(monkeypatch):
    """Supply independent toy model/backend identities with controlled exact source revisions.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the provider or environment under test, restored after the
        case.

    Returns
    -------
    dict
        Metadata input names that let tests isolate schema behavior from the current
        checkout.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    monkeypatch.setattr(
        "itamae.provenance._required_source_revision", lambda *args, **kwargs: "a" * 40
    )
    monkeypatch.setattr("itamae.provenance.package_version", lambda *args, **kwargs: "0.2.0rc1")
    return dict(
        variant="toy",
        distribution_name="toy",
        module_file=__file__,
        model_identifier="toy:model",
        backend_identifier="toy:backend",
        source_identifier="toy:source",
        variance_identifier="toy:variance",
        power_identifier="toy:power",
        solver_identifier="toy:solver",
    )


def test_calculation_schema_roundtrip_does_not_relabel_history(identifiers, tmp_path):
    """Keep old and current physical prescriptions distinguishable after saving a catalogue.

    Parameters
    ----------
    identifiers : dict
        Controlled model, backend, spectrum and solver identities used to exercise metadata
        without changing physical calculations.
    tmp_path : pathlib.Path
        Isolated directory for the test catalogue, cache or provenance files.

    Notes
    -----
    Pass criteria: Historical metadata retain physics_mode, current metadata retain
    calculation_specification, and both survive NPZ round trips exactly.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    historical = build_migration_metadata(**identifiers, physics_mode="legacy")
    current = build_calculation_metadata(**identifiers, calculation_specification="toy:v1")
    assert "physics_mode" not in current.as_mapping()
    assert set(CALCULATION_METADATA_KEYS) <= set(current.as_mapping())
    assert historical.extra["physics_mode"] == "legacy"
    assert "calculation_specification" not in historical.as_mapping()
    for name, metadata in (("old", historical), ("new", current)):
        catalog = WeightedSubhaloCatalog(
            columns={"mass": np.array([1.0, 2.0])},
            weights={"weight_base": np.array([0.5, 0.25])},
            metadata=metadata,
        )
        path = tmp_path / (name + ".npz")
        catalog.to_npz(path)
        restored = WeightedSubhaloCatalog.from_npz(path)
        assert dict(restored.metadata) == dict(metadata.as_mapping())


@pytest.mark.parametrize(
    "extra",
    [
        {"physics_mode": "legacy"},
        {"calculation_specification": "wrong:v1"},
        {"provenance_schema": "old"},
    ],
)
def test_new_provenance_rejects_conflicting_identity(identifiers, extra):
    """Prevent auxiliary metadata from overriding the declared calculation prescription.

    Parameters
    ----------
    identifiers : dict
        Controlled model, backend, spectrum and solver identities used to exercise metadata
        without changing physical calculations.
    extra : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Conflicting physics_mode, specification, or provenance-schema fields
    raise ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError):
        build_calculation_metadata(**identifiers, calculation_specification="toy:v1", extra=extra)


@pytest.mark.parametrize("specification", [None, "", "  ", 3])
def test_new_provenance_requires_a_specification(identifiers, specification):
    """Require an identifiable physical prescription for newly generated catalogues.

    Parameters
    ----------
    identifiers : dict
        Controlled model, backend, spectrum and solver identities used to exercise metadata
        without changing physical calculations.
    specification : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Missing, blank, or non-string calculation specifications raise
    ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError):
        build_calculation_metadata(**identifiers, calculation_specification=specification)
