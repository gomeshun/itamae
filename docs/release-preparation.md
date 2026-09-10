# Release candidate and publication handoff

Candidate version: **0.2.0rc1**. Distribution: `sashimi-itamae`; import: `itamae`.
The family work-branch manifest records the reviewed source combination after
artifact verification. Candidate artifacts are supplied through a local wheelhouse;
publication on PyPI is not implied by successful local dependency resolution.

## Reproduce the candidate

1. Use a clean checkout of the exact reviewed SHA, on the migration branch or its
   topic branch. Record the Python/dependency environment.
2. Run Ruff, format, mypy, pytest/coverage and the Python 3.11–3.13 backend matrix.
3. Build wheel and sdist with `python -m build` or `uv build`. The build hook reads
   this repository's HEAD and embeds it. Record both artifact SHA-256 hashes.
4. Unpack the sdist outside a Git checkout and run `uv build --wheel` there without
   `ITAMAE_SOURCE_REVISION`. Install that rebuilt wheel in an empty environment;
   verify `itamae.provenance.source_revision("itamae")` equals the source SHA.
5. Install the reviewed family artifacts in clean environments outside all source
   trees. Run catalogs/observables, provenance and runtime-file collision checks.

## Separate publication procedure after peer review

Main integration, final artifact recreation, public repository settings, index
permissions and upload are subsequent operations requiring explicit user approval.
Rebuild from the actual final SHA and repeat affected regressions plus the exact
five-package family matrix. Do not reuse this candidate's artifact evidence for
changed source. Confirm name ownership and trusted publisher settings externally.

Publish this core first, verify a clean installation from the public index, then
publish C/SI/W/F with their versioned core dependency. Check metadata and hashes,
then perform the source-tree-free family smoke again against index artifacts.
TestPyPI is also publication and requires separate authorization.

`release-workflow.yml.example` is an inactive workflow template. It is deliberately
outside `.github/workflows`; preparing it does not register an upload trigger.
Before later enabling it, configure a protected `pypi` environment with required
human reviewers, OIDC trusted publishing and an approved final source SHA. Both
manual dispatch and environment approval are required. No release tags, scheduled
uploads or automatic publish jobs are part of preparation.
