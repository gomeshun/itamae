"""Tests for standardized migration provenance metadata."""

from pathlib import Path

import pytest

from itamae.provenance import (
    MIGRATION_METADATA_KEYS,
    build_migration_metadata,
)


def test_build_migration_metadata_includes_standard_keys(monkeypatch) -> None:
    """The builder records every shared field and preserves variant extras."""
    monkeypatch.setenv("ITAMAE_SOURCE_REVISION", "i" * 40)
    monkeypatch.setenv("TEST_MODEL_SOURCE_REVISION", "s" * 40)

    metadata = build_migration_metadata(
        variant="test-model",
        distribution_name="test-model",
        module_file=str(Path(__file__)),
        model_identifier="test-model:legacy:v1",
        backend_identifier="array=numpy",
        source_identifier="test-model:adapter:v1",
        physics_mode="legacy",
        variance_identifier="test-model:variance:v1",
        power_identifier="test-model:power:v1",
        solver_identifier="test-model:solver:v1",
        extra={"variant_parameter": 3},
    )

    values = metadata.as_mapping()
    assert set(MIGRATION_METADATA_KEYS) <= set(values)
    assert values["itamae_source_revision"] == "i" * 40
    assert values["sashimi_source_revision"] == "s" * 40
    assert values["variant_parameter"] == 3
    assert values["catalog_schema_version"] == metadata.schema_version


def test_build_migration_metadata_rejects_standard_field_override() -> None:
    """Variant-specific extras may not silently redefine shared fields."""
    with pytest.raises(ValueError, match="conflicting standard fields"):
        build_migration_metadata(
            variant="test-model",
            distribution_name="test-model",
            module_file=str(Path(__file__)),
            model_identifier="test-model:consistent:v1",
            backend_identifier="array=numpy",
            source_identifier="test-model:adapter:v1",
            physics_mode="consistent",
            variance_identifier="test-model:variance:v1",
            power_identifier="test-model:power:v1",
            solver_identifier="test-model:solver:v1",
            extra={"physics_mode": "legacy"},
        )
