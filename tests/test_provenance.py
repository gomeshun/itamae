"""Tests for standardized migration provenance metadata."""

from pathlib import Path

import pytest

import itamae.provenance as provenance
from itamae.provenance import (
    MIGRATION_METADATA_KEYS,
    build_migration_metadata,
    source_revision,
)


def test_build_migration_metadata_includes_standard_keys(monkeypatch) -> None:
    """The builder records every shared field and preserves variant extras."""
    monkeypatch.delenv("ITAMAE_SOURCE_REVISION", raising=False)

    metadata = build_migration_metadata(
        variant="test-model",
        distribution_name="itamae",
        module_file=str(provenance.__file__),
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
    expected_revision = source_revision("itamae", module_file=str(provenance.__file__))
    assert len(expected_revision) == 40
    assert values["itamae_source_revision"] == expected_revision
    assert values["sashimi_source_revision"] == expected_revision
    assert values["variant_parameter"] == 3
    assert values["catalog_schema_version"] == metadata.schema_version


def test_build_migration_metadata_rejects_standard_field_override() -> None:
    """Variant-specific extras may not silently redefine shared fields."""
    with pytest.raises(ValueError, match="conflicting standard fields"):
        build_migration_metadata(
            variant="test-model",
            distribution_name="itamae",
            module_file=str(provenance.__file__),
            model_identifier="test-model:consistent:v1",
            backend_identifier="array=numpy",
            source_identifier="test-model:adapter:v1",
            physics_mode="consistent",
            variance_identifier="test-model:variance:v1",
            power_identifier="test-model:power:v1",
            solver_identifier="test-model:solver:v1",
            extra={"physics_mode": "legacy"},
        )


def test_embedded_revision_takes_precedence_over_ambient_environment(monkeypatch) -> None:
    """Installed-wheel metadata cannot be changed by an environment variable."""
    expected = "a" * 40
    monkeypatch.setenv("ITAMAE_SOURCE_REVISION", "b" * 40)
    monkeypatch.setattr(provenance, "_source_checkout_root", lambda module_file: None)
    monkeypatch.setattr(provenance, "_embedded_source_revision", lambda package_name: expected)
    assert source_revision("itamae", module_file="/outside/site-packages/provenance.py") == expected


def test_recognized_distribution_cannot_finish_with_unknown_revision(monkeypatch) -> None:
    """Durable provenance is mandatory for recognized family distributions."""
    monkeypatch.setattr(provenance, "_source_checkout_root", lambda module_file: None)
    monkeypatch.setattr(provenance, "_embedded_source_revision", lambda package_name: None)
    monkeypatch.setattr(provenance, "_direct_url_revision", lambda package_name: None)
    with pytest.raises(RuntimeError, match="No durable source revision"):
        source_revision("itamae")


def test_source_revision_does_not_walk_into_an_outer_repository_from_venv(
    tmp_path: Path, monkeypatch
) -> None:
    """A package under another repo's .venv must not inherit that repo's HEAD."""
    outer = tmp_path / "other-repository"
    module_file = outer / ".venv/lib/python3.11/site-packages/fake.py"
    module_file.parent.mkdir(parents=True)
    module_file.touch()
    (outer / ".git").mkdir()
    (outer / "pyproject.toml").write_text("[project]\nversion = '0.0.0'\n")

    def unexpected_git(*args, **kwargs):
        raise AssertionError("ambient repository Git lookup")

    monkeypatch.setattr(provenance.subprocess, "run", unexpected_git)
    monkeypatch.setattr(provenance, "_embedded_source_revision", lambda package_name: None)
    monkeypatch.setattr(provenance, "_direct_url_revision", lambda package_name: None)
    assert source_revision("missing-distribution", module_file=str(module_file)) == "unknown"
