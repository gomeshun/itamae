"""Tests for standardized migration provenance metadata."""

from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

import itamae.provenance as provenance
from itamae.provenance import (
    ITAMAE_DISTRIBUTION_NAME,
    MIGRATION_METADATA_KEYS,
    build_migration_metadata,
    package_version,
    source_revision,
)


@pytest.mark.parametrize(
    "module_path",
    ["src/itamae/provenance.py", "src/sashimi_c/__init__.py", "src/sashimi_f/_physics.py"],
)
def test_src_package_provenance_uses_its_own_repository(tmp_path, module_path):
    """Locate the source checkout that actually owns a core or variant calculation.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated directory for the test catalogue, cache or provenance files.
    module_path : str
        Representative package path inside a source checkout.

    Notes
    -----
    Pass criteria: Each supported src-layout module resolves to its own temporary project
    root.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    subprocess.run(["git", "init", "--quiet", str(tmp_path)], check=True)
    (tmp_path / "pyproject.toml").write_text('[project]\nname="example"\nversion="1.0"\n')
    module = tmp_path / module_path
    module.parent.mkdir(parents=True)
    module.write_text("# source module\n")
    assert provenance._source_checkout_root(str(module)) == tmp_path


def test_src_export_does_not_inherit_an_unrelated_parent_repository(tmp_path):
    """Avoid attributing an unpacked scientific package to an unrelated enclosing project.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated directory for the test catalogue, cache or provenance files.

    Notes
    -----
    Pass criteria: An exported src-layout package under another Git repository has no
    source-checkout root.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    subprocess.run(["git", "init", "--quiet", str(tmp_path)], check=True)
    exported = tmp_path / "unpacked"
    module = exported / "src" / "sashimi_c" / "__init__.py"
    module.parent.mkdir(parents=True)
    module.write_text("# exported module\n")
    (exported / "pyproject.toml").write_text('[project]\nname="example"\nversion="1.0"\n')
    assert provenance._source_checkout_root(str(module)) is None


def test_build_migration_metadata_includes_standard_keys(monkeypatch) -> None:
    """Record the sources and physical prescription needed to interpret a historical migration
    result.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.

    Notes
    -----
    Pass criteria: Metadata contain all shared keys, exact resolved source SHAs, schema
    identity and the supplied extra parameter.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
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
    """Prevent variant extras from changing the stated historical physical mode.

    Notes
    -----
    Pass criteria: An extra physics_mode conflicting with the declared mode raises
    ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
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


def test_logical_itamae_name_resolves_renamed_distribution(monkeypatch) -> None:
    """Resolve the core's logical name to the installed distribution used for provenance.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.

    Notes
    -----
    Pass criteria: Version lookup for itamae queries sashimi-itamae and returns the supplied
    version.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    seen = []

    def fake_version(name: str) -> str:
        """Record the distribution name requested by a controlled version lookup.

        Parameters
        ----------
        name : str
            Name of the quantity or component included in an error or provenance record so an
            invalid scientific input can be located.

        Returns
        -------
        str
            Known toy version, allowing the test to assert the exact lookup identity.

        Notes
        -----
        Resolve the core's logical name to the installed distribution used for provenance. Pass
        criteria: Version lookup for itamae queries sashimi-itamae and returns the supplied
        version. Failure occurs when the stated comparison or expected rejection is not met.
        Source attribution and format compatibility do not validate the physical prescription or
        numerical accuracy.
        """
        seen.append(name)
        return "9.9.9"

    monkeypatch.setattr(provenance.importlib_metadata, "version", fake_version)
    assert package_version("itamae") == "9.9.9"
    assert seen == [ITAMAE_DISTRIBUTION_NAME]


@pytest.mark.parametrize(
    "package_name",
    ["SASHIMI-ITAMAE", "sashimi_itamae", "sashimi.itamae"],
)
def test_distribution_alias_normalization_for_version_lookup(
    monkeypatch, package_name: str
) -> None:
    """Treat supported distribution-name spellings as one installed core version.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.
    package_name : str
        Supported logical or distribution-name spelling used for source/version attribution.

    Notes
    -----
    Pass criteria: Case, hyphen, underscore and dot spellings all query the canonical
    distribution and return one version.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    seen = []

    def fake_version(name: str) -> str:
        """Record the distribution name requested by a controlled version lookup.

        Parameters
        ----------
        name : str
            Name of the quantity or component included in an error or provenance record so an
            invalid scientific input can be located.

        Returns
        -------
        str
            Known toy version, allowing the test to assert the exact lookup identity.

        Notes
        -----
        Treat supported distribution-name spellings as one installed core version. Pass
        criteria: Case, hyphen, underscore and dot spellings all query the canonical
        distribution and return one version. Failure occurs when the stated comparison or
        expected rejection is not met. Source attribution and format compatibility do not
        validate the physical prescription or numerical accuracy.
        """
        seen.append(name)
        return "9.9.9"

    monkeypatch.setattr(provenance.importlib_metadata, "version", fake_version)
    assert package_version(package_name) == "9.9.9"
    assert seen == [ITAMAE_DISTRIBUTION_NAME]


def test_distribution_name_alias_uses_itamae_embedded_provenance(monkeypatch) -> None:
    """Use the same embedded scientific source for logical and distribution names.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.

    Notes
    -----
    Pass criteria: The public distribution alias resolves to the expected embedded SHA.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    expected = "a" * 40
    monkeypatch.setattr(provenance, "_source_checkout_root", lambda module_file: None)
    monkeypatch.setattr(provenance, "_embedded_source_revision", lambda package_name: expected)
    assert source_revision(ITAMAE_DISTRIBUTION_NAME) == expected


@pytest.mark.parametrize(
    "package_name",
    ["SASHIMI-ITAMAE", "sashimi_itamae", "sashimi.itamae"],
)
def test_distribution_alias_normalization_for_embedded_provenance(
    monkeypatch, package_name: str
) -> None:
    """Keep source identity stable under supported package-name spellings.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.
    package_name : str
        Supported logical or distribution-name spelling used for source/version attribution.

    Notes
    -----
    Pass criteria: All aliases load itamae._build_provenance and return the expected SHA.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    expected = "a" * 40
    imported = []

    def fake_import_module(module_name: str):
        """Record a controlled import of embedded scientific source attribution.

        Parameters
        ----------
        module_name : str
            Name of the embedded attribution module expected by the controlled import lookup.

        Returns
        -------
        types.SimpleNamespace
            Known SOURCE_REVISION descriptor.

        Notes
        -----
        Keep source identity stable under supported package-name spellings. Pass criteria: All
        aliases load itamae._build_provenance and return the expected SHA. Failure occurs when
        the stated comparison or expected rejection is not met. Source attribution and format
        compatibility do not validate the physical prescription or numerical accuracy.
        """
        imported.append(module_name)
        return SimpleNamespace(SOURCE_REVISION=expected)

    monkeypatch.setattr(provenance, "_source_checkout_root", lambda module_file: None)
    monkeypatch.setattr(provenance.importlib, "import_module", fake_import_module)
    assert source_revision(package_name) == expected
    assert imported == ["itamae._build_provenance"]


def test_embedded_revision_takes_precedence_over_ambient_environment(monkeypatch) -> None:
    """Protect installed scientific provenance from unrelated shell settings.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.

    Notes
    -----
    Pass criteria: A conflicting environment SHA cannot replace the wheel's embedded source
    SHA.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    expected = "a" * 40
    monkeypatch.setenv("ITAMAE_SOURCE_REVISION", "b" * 40)
    monkeypatch.setattr(provenance, "_source_checkout_root", lambda module_file: None)
    monkeypatch.setattr(provenance, "_embedded_source_revision", lambda package_name: expected)
    assert source_revision("itamae", module_file="/outside/site-packages/provenance.py") == expected


def test_recognized_distribution_cannot_finish_with_unknown_revision(monkeypatch) -> None:
    """Refuse to record an unidentified installed family calculation as reproducible.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.

    Notes
    -----
    Pass criteria: Recognized core aliases raise RuntimeError when no durable source
    revision is available.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    monkeypatch.setattr(provenance, "_source_checkout_root", lambda module_file: None)
    monkeypatch.setattr(provenance, "_embedded_source_revision", lambda package_name: None)
    monkeypatch.setattr(provenance, "_direct_url_revision", lambda package_name: None)
    with pytest.raises(RuntimeError, match="No durable source revision"):
        source_revision("itamae")
    with pytest.raises(RuntimeError, match="No durable source revision"):
        source_revision(ITAMAE_DISTRIBUTION_NAME)
    with pytest.raises(RuntimeError, match="No durable source revision"):
        source_revision("SASHIMI_ITAMAE")


def test_source_revision_does_not_walk_into_an_outer_repository_from_venv(
    tmp_path: Path, monkeypatch
) -> None:
    """Avoid assigning an application's revision to a dependency installed in its environment.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated directory for the test catalogue, cache or provenance files.
    monkeypatch : pytest.MonkeyPatch
        Scoped replacement of the particular physical provider or environment under test;
        restored after the case.

    Notes
    -----
    Pass criteria: An unidentified package inside another project's virtual environment
    returns unknown without invoking ambient Git.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: Source attribution and format compatibility do not validate the physical
    prescription or numerical accuracy.
    """
    outer = tmp_path / "other-repository"
    module_file = outer / ".venv/lib/python3.11/site-packages/fake.py"
    module_file.parent.mkdir(parents=True)
    module_file.touch()
    (outer / ".git").mkdir()
    (outer / "pyproject.toml").write_text("[project]\nversion = '0.0.0'\n")

    def unexpected_git(*args, **kwargs):
        """Fail if an installed package tries to borrow an enclosing repository revision.

        Parameters
        ----------
        *args : tuple
            Additional positional parameters forwarded unchanged to the supplied numerical or
            physical callback.
        **kwargs : dict
            Additional keyword options forwarded to the called routine; reserved scientific
            controls remain owned by that routine.

        Notes
        -----
        Avoid assigning an application's revision to a dependency installed in its environment.
        Pass criteria: An unidentified package inside another project's virtual environment
        returns unknown without invoking ambient Git. Failure occurs when the stated comparison
        or expected rejection is not met. Source attribution and format compatibility do not
        validate the physical prescription or numerical accuracy.
        """
        raise AssertionError("ambient repository Git lookup")

    monkeypatch.setattr(provenance.subprocess, "run", unexpected_git)
    monkeypatch.setattr(provenance, "_embedded_source_revision", lambda package_name: None)
    monkeypatch.setattr(provenance, "_direct_url_revision", lambda package_name: None)
    assert source_revision("missing-distribution", module_file=str(module_file)) == "unknown"
