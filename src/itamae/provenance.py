"""Provenance helpers for migrated model catalogs."""

from __future__ import annotations

import importlib
import json
import re
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
SOURCE_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
_EMBEDDED_SOURCE_MODULES = {
    "itamae": "itamae._build_provenance",
    "sashimi-c": "_sashimi_c_build_provenance",
    "sashimi-si": "_sashimi_si_build_provenance",
    "sashimi-w": "_sashimi_w_build_provenance",
    "sashimi-f": "_sashimi_f_build_provenance",
}


def _valid_source_revision(value: Any) -> str | None:
    """Return a full lowercase commit SHA or ``None`` for invalid values."""
    if isinstance(value, str) and SOURCE_REVISION_PATTERN.fullmatch(value):
        return value
    return None


def _source_checkout_root(module_file: str | None) -> Path | None:
    """Return a repository root only for a recognized source layout."""
    if module_file is None:
        return None
    module_path = Path(module_file).resolve()
    for directory in (module_path.parent, *module_path.parents):
        project_file = directory / "pyproject.toml"
        if not project_file.is_file():
            continue
        relative_path = module_path.relative_to(directory)
        if not relative_path.parts:
            continue
        is_src_package = relative_path.parts[0] == "src" and (directory / "src" / "itamae").is_dir()
        is_top_level_module = len(relative_path.parts) == 1 and module_path.suffix == ".py"
        if not (is_src_package or is_top_level_module):
            continue
        try:
            result = subprocess.run(
                [
                    "git",
                    "-C",
                    str(directory),
                    "rev-parse",
                    "--show-toplevel",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
        except OSError:
            return None
        if result.returncode != 0:
            return None
        git_root = Path(result.stdout.strip()).resolve()
        return directory if git_root == directory.resolve() else None
    return None


def _project_version(module_file: str | None) -> str | None:
    """Read a source checkout's project version when package metadata is absent."""
    source_root = _source_checkout_root(module_file)
    if source_root is None:
        return None
    try:
        with (source_root / "pyproject.toml").open("rb") as stream:
            project = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError):
        return None
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
    return _valid_source_revision(vcs_info.get("commit_id"))


def _embedded_source_revision(package_name: str) -> str | None:
    """Read a source revision embedded in an installed distribution."""
    module_name = _EMBEDDED_SOURCE_MODULES.get(package_name)
    if module_name is None:
        return None
    try:
        module = importlib.import_module(module_name)
    except (ImportError, ModuleNotFoundError):
        return None
    return _valid_source_revision(getattr(module, "SOURCE_REVISION", None))


def _git_revision(module_file: str | None) -> str | None:
    """Return the current Git revision for a module in a source checkout."""
    source_root = _source_checkout_root(module_file)
    if source_root is None:
        return None
    try:
        result = subprocess.run(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return _valid_source_revision(result.stdout.strip())


def source_revision(
    package_name: str,
    *,
    module_file: str | None = None,
) -> str:
    """Resolve a reproducible source revision without using ambient Git state.

    Source checkouts use Git only when ``module_file`` belongs directly to the
    repository's recognized source layout. Installed wheels use their embedded
    build-provenance module, so a parent repository or environment variable
    cannot change their reported revision. VCS-installed distributions use
    ``direct_url.json`` as a fallback.
    """
    if _source_checkout_root(module_file) is not None:
        revision = _git_revision(module_file)
        if revision is not None:
            return revision
    revision = _embedded_source_revision(package_name)
    if revision is not None:
        return revision
    revision = _direct_url_revision(package_name)
    return revision or UNKNOWN_SOURCE_REVISION


def _required_source_revision(package_name: str, *, module_file: str) -> str:
    """Require an exact source revision for catalog metadata construction."""
    revision = source_revision(package_name, module_file=module_file)
    if revision == UNKNOWN_SOURCE_REVISION:
        raise RuntimeError(
            f"No exact source revision is available for {package_name!r}; "
            "build the distribution with embedded provenance."
        )
    return revision


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
        "itamae_source_revision": _required_source_revision("itamae", module_file=__file__),
        "sashimi_version": package_version(distribution_name, module_file=module_file),
        "sashimi_source_revision": _required_source_revision(
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
