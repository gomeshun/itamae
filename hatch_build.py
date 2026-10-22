"""Embed the ITAMAE source revision in build artifacts."""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

_SOURCE_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
_TARGET = "itamae/_build_provenance.py"


def _valid_revision(value: str | None) -> str | None:
    """Return a full lowercase source SHA or ``None``.

    Parameters
    ----------
    value : array_like
        Numerical input to be checked before use; shape and physical interpretation are set
        by the caller and retained by the validation helper.

    Returns
    -------
    str or None
        Full lowercase 40-character source revision if valid, otherwise None.

    Notes
    -----
    Build attribution must name this package source, even inside another Git repository. It
    supports tracing a scientific result from an installed wheel or source archive and does
    not measure physical accuracy.
    """
    if value is not None:
        value = value.strip()
    return value if value and _SOURCE_REVISION_PATTERN.fullmatch(value) else None


class CustomBuildHook(BuildHookInterface):
    """Add an exact source revision to each Hatchling artifact.

    Notes
    -----
    Build attribution must name this package source, even inside another Git repository. It
    supports tracing a scientific result from an installed wheel or source archive and does
    not measure physical accuracy.
    """

    _generated_path: Path | None = None

    def _revision_from_existing_artifact(self) -> str | None:
        """Preserve a revision when rebuilding from an exported source archive.

        Returns
        -------
        str or None
            Revision embedded in the package archive, or None if no valid descriptor exists.

        Notes
        -----
        Build attribution must name this package source, even inside another Git repository. It
        supports tracing a scientific result from an installed wheel or source archive and does
        not measure physical accuracy.
        """
        candidates = (
            Path(self.root) / "src" / _TARGET,
            Path(self.root) / _TARGET,
        )
        for path in candidates:
            if not path.is_file():
                continue
            match = re.search(r"SOURCE_REVISION\s*=\s*['\"]([0-9a-f]{40})['\"]", path.read_text())
            if match:
                return match.group(1)
        return None

    def _revision_from_git_checkout(self) -> str | None:
        """Use Git only when the package root is the repository root.

        Returns
        -------
        str or None
            Current package checkout revision, or None if this is not its own Git root.

        Notes
        -----
        Build attribution must name this package source, even inside another Git repository. It
        supports tracing a scientific result from an installed wheel or source archive and does
        not measure physical accuracy.
        """
        try:
            root_result = subprocess.run(
                ["git", "-C", self.root, "rev-parse", "--show-toplevel"],
                check=False,
                capture_output=True,
                text=True,
            )
            if root_result.returncode != 0:
                return None
            if Path(root_result.stdout.strip()).resolve() != Path(self.root).resolve():
                return None
            result = subprocess.run(
                ["git", "-C", self.root, "rev-parse", "HEAD"],
                check=False,
                capture_output=True,
                text=True,
            )
        except OSError:
            return None
        return _valid_revision(result.stdout if result.returncode == 0 else None)

    def _resolve_revision(self) -> str:
        """Preserve archive identity and reject conflicting build overrides.

        Returns
        -------
        str
            Exact revision supported by archive/checkout identity or an explicit non-conflicting
            override.

        Notes
        -----
        Build attribution must name this package source, even inside another Git repository. It
        supports tracing a scientific result from an installed wheel or source archive and does
        not measure physical accuracy.
        """
        supplied = os.environ.get("ITAMAE_SOURCE_REVISION")
        explicit = _valid_revision(supplied)
        if supplied is not None and explicit is None:
            raise RuntimeError(
                "ITAMAE_SOURCE_REVISION must contain a full 40-character lowercase SHA."
            )
        revision = self._revision_from_existing_artifact()
        if revision is None:
            revision = self._revision_from_git_checkout()
        if revision is not None and explicit is not None and revision != explicit:
            raise RuntimeError(
                "Explicit source revision conflicts with the source archive or checkout."
            )
        revision = revision or explicit
        if revision is None:
            raise RuntimeError(
                "ITAMAE artifact builds require an exact source revision from "
                "the source archive, ITAMAE_SOURCE_REVISION, or its own Git checkout."
            )
        return revision

    def initialize(self, version: str, build_data: dict[str, object]) -> None:
        """Generate and force-include the build provenance module.

        Parameters
        ----------
        version : str
            Build-system artifact version label; it does not replace the exact source revision.
        build_data : dict
            Build-system metadata updated to include the generated source-attribution module.

        Notes
        -----
        Build attribution must name this package source, even inside another Git repository. It
        supports tracing a scientific result from an installed wheel or source archive and does
        not measure physical accuracy.
        """
        revision = self._resolve_revision()
        descriptor_fd, descriptor_name = tempfile.mkstemp(
            prefix="itamae-build-provenance-", suffix=".py"
        )
        os.close(descriptor_fd)
        descriptor = Path(descriptor_name)
        descriptor.write_text(f"SOURCE_REVISION = {revision!r}\n", encoding="utf-8")
        self._generated_path = descriptor
        force_include = build_data.setdefault("force_include", {})
        force_include[str(descriptor)] = _TARGET

    def finalize(
        self,
        version: str,
        build_data: dict[str, object],
        artifact_path: str,
    ) -> None:
        """Remove the temporary source file after the artifact is written.

        Parameters
        ----------
        version : str
            Build-system artifact version label; it does not replace the exact source revision.
        build_data : dict
            Build-system metadata updated to include the generated source-attribution module.
        artifact_path : str
            Path of the completed distribution artifact passed by the build system.

        Notes
        -----
        Build attribution must name this package source, even inside another Git repository. It
        supports tracing a scientific result from an installed wheel or source archive and does
        not measure physical accuracy.
        """
        if self._generated_path is not None:
            self._generated_path.unlink(missing_ok=True)
            self._generated_path = None
