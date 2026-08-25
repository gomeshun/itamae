"""Provenance helpers for migrated model catalogs."""

from __future__ import annotations

import json
import os
import subprocess
import tomllib
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Mapping

from . import __version__
from .types import CATALOG_SCHEMA_VERSION, CatalogMetadata
from .units import CANONICAL_UNIT_SCHEMA_VERSION

MIGRATION_METADATA_KEYS = (
    "sashimi_variant",
    "physics_mode",
    "itamae_version",
    "itamae_source_revision",
    "sashimi_version",
    "sashimi_source_revision",
    "canonical_unit_schema",
    "variance_identifier",
    "power_identifier",
    "solver_identifier",
    "catalog_schema_version",
)
UNKNOWN_SOURCE_REVISION = "unknown"


def _environment_variable(package_name: str) -> str:
    """Return the conventional source-revision environment variable name."""
    normalized = "".join(character if character.isalnum() else "_" for character in package_name)
    return f"{normalized.upper()}_SOURCE_REVISION"


def _project_version(module_file: str | None) -> str | None:
    """Read a source checkout's project version when package metadata is absent."""
    if module_file is None:
        return None
    module_path = Path(module_file).resolve()
    for directory in (module_path.parent, *module_path.parents):
        project_file = directory / "pyproject.toml"
        if not project_file.is_file():
            continue
        try:
            with project_file.open("rb") as stream:
                project = tomllib.load(stream)
        except (OSError, tomllib.TOMLDecodeError):
            continue
        version = project.get("project", {}).get("version")
        if isinstance(version, str) and version.strip():
            return version
    return None


def package_version(package_name: str, *, module_file: str | None = None) -> str:
    """Resolve a distribution version from installed or source-project metadata."""
    try:
        return importlib_metadata.version(package_name)
    except importlib_metadata.PackageNotFoundError:
        return _project_version(module_file) or UNKNOWN_SOURCE_REVISION


def _direct_url_revision(package_name: str) -> str | None:
    """Return a VCS commit recorded in a distribution's direct-url metadata."""
    try:
        direct_url = importlib_metadata.distribution(package_name).read_text("direct_url.json")
    except (OSError, importlib_metadata.PackageNotFoundError):
        return None
    if not direct_url:
        return None
    try:
        payload = json.loads(direct_url)
    except json.JSONDecodeError:
        return None
    vcs_info = payload.get("vcs_info")
    if not isinstance(vcs_info, dict):
        return None
    revision = vcs_info.get("commit_id")
    return revision.strip() if isinstance(revision, str) and revision.strip() else None


def _git_revision(module_file: str | None) -> str | None:
    """Return the current Git revision for a module in a source checkout."""
    if module_file is None:
        return None
    try:
        result = subprocess.run(
            ["git", "-C", str(Path(module_file).resolve().parent), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return None
    if result.returncode != 0:
        return None
    revision = result.stdout.strip()
    return revision or None


def source_revision(
    package_name: str,
    *,
    module_file: str | None = None,
    environment_variable: str | None = None,
) -> str:
    """Resolve a reproducible source revision without requiring Git at runtime.

    CI can provide an exact revision through the conventional environment
    variable. VCS-installed distributions and source checkouts are supported as
    fallbacks; wheels built without either source record return ``"unknown"``.
    """
    environment_names = (
        (environment_variable, _environment_variable(package_name))
        if environment_variable is not None
        else (_environment_variable(package_name),)
    )
    for name in environment_names:
        value = os.environ.get(name)
        if value and value.strip():
            return value.strip()

    revision = _git_revision(module_file)
    if revision is not None:
        return revision
    revision = _direct_url_revision(package_name)
    return revision or UNKNOWN_SOURCE_REVISION


def build_migration_metadata(
    *,
    variant: str,
    distribution_name: str,
    module_file: str,
    model_identifier: str,
    backend_identifier: str,
    source_identifier: str,
    physics_mode: str,
    variance_identifier: str,
    power_identifier: str,
    solver_identifier: str,
    canonical_unit_schema: str = CANONICAL_UNIT_SCHEMA_VERSION,
    extra: Mapping[str, Any] | None = None,
) -> CatalogMetadata:
    """Build a catalog metadata object with the shared migration vocabulary."""
    standard = {
        "sashimi_variant": variant,
        "physics_mode": physics_mode,
        "itamae_version": __version__,
        "itamae_source_revision": source_revision("itamae", module_file=__file__),
        "sashimi_version": package_version(distribution_name, module_file=module_file),
        "sashimi_source_revision": source_revision(
            distribution_name,
            module_file=module_file,
        ),
        "canonical_unit_schema": canonical_unit_schema,
        "variance_identifier": variance_identifier,
        "power_identifier": power_identifier,
        "solver_identifier": solver_identifier,
        "catalog_schema_version": CATALOG_SCHEMA_VERSION,
    }
    values = dict(extra or {})
    conflicting = {
        name: values[name] for name in standard if name in values and values[name] != standard[name]
    }
    if conflicting:
        raise ValueError(
            f"Migration metadata contains conflicting standard fields: {sorted(conflicting)}."
        )
    values.update(standard)
    return CatalogMetadata(
        model_identifier=model_identifier,
        backend_identifier=backend_identifier,
        source_identifier=source_identifier,
        schema_version=CATALOG_SCHEMA_VERSION,
        extra=values,
    )


__all__ = [
    "MIGRATION_METADATA_KEYS",
    "UNKNOWN_SOURCE_REVISION",
    "build_migration_metadata",
    "package_version",
    "source_revision",
]
