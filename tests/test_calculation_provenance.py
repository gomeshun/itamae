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
    with pytest.raises(ValueError):
        build_calculation_metadata(**identifiers, calculation_specification="toy:v1", extra=extra)


@pytest.mark.parametrize("specification", [None, "", "  ", 3])
def test_new_provenance_requires_a_specification(identifiers, specification):
    with pytest.raises(ValueError):
        build_calculation_metadata(**identifiers, calculation_specification=specification)
