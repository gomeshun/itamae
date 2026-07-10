# ITAMAE implementation plan

**ITAMAE: Integrated Toolkit for Analytical Merger-tree And Evolution**

This document defines the implementation policy for ITAMAE after reviewing the
current SASHIMI-C, SASHIMI-W, SASHIMI-SI, and SASHIMI-F implementations, the
FDM profile implementation in `dsph_fuzzy`, and the spatial development in the
`r-dependent` branch of SASHIMI-C.

ITAMAE is intended to be published on PyPI and installed and developed with
modern Python tooling, including `uv`. Scientific reproducibility, automated
verification, and clear documentation are first-class requirements rather than
post-release additions.

## 1. Executive decision

ITAMAE owns the common **computational language** of the SASHIMI family:

1. typed state and catalog objects;
2. backend-independent cosmology and unit interfaces;
3. spherical-halo and profile primitives;
4. grids, quadrature, interpolation, integration, caching, and batch execution;
5. generic evolution solvers operating on model-supplied equations;
6. deterministic weighted measures and stochastic realizations;
7. optional radial and phase-space infrastructure;
8. packaging, versioning, testing, documentation, and distribution machinery.

Each SASHIMI variant continues to own its physical prescriptions and default
composition:

- **SASHIMI-C:** CDM variance/concentration choices, Yang-type accretion model,
  calibrated average stripping law, structural response, disruption convention,
  and CDM observables;
- **SASHIMI-W:** WDM transfer function, filtering convention, WDM concentration,
  and WDM-specific calibration;
- **SASHIMI-SI:** SIDM cross section and gravothermal/profile evolution;
- **SASHIMI-F:** FDM transfer function, variance/filter choices, and FDM-specific
  structure assignment;
- **downstream dSph projects:** Jeans likelihoods, observational priors, and
  inference workflows.

> ITAMAE defines how components communicate and how a calculation is executed.  
> A SASHIMI variant defines which physical components constitute its model.

ITAMAE must not provide a monolithic default `SashimiModel` whose hidden choices
silently reproduce SASHIMI-C.

## 2. Common computational pipeline

All current SASHIMI variants follow the same conceptual sequence:

1. specify cosmology and a host halo at a reference redshift;
2. construct an accretion measure in mass and redshift;
3. evaluate an EPS-like accretion abundance;
4. integrate over host-history and concentration scatter;
5. assign internal structure at accretion;
6. evolve mass and structure after accretion;
7. evaluate survival, disruption, collapse, or validity conditions;
8. assemble a weighted catalog;
9. calculate observables, priors, or stochastic realizations.

The shared pipeline and its data contracts are the main objects ITAMAE should
support. Individual physical formulae remain explicit model components.

## 3. Architectural model

The basic object is a weighted measure transported through physical maps.

```text
x = (m_acc, z_acc, concentration_node, host_history_node, orbit_node, ...)
y = T_model(x)
expectation[O] = sum_i w_i O(y_i)
```

The existing SASHIMI catalog is a deterministic quadrature representation of a
population measure. ITAMAE preserves this as the canonical population
representation. Monte Carlo catalogs are generated as a separate operation.

The architecture has four layers:

1. **configuration and physical models** supplied by a SASHIMI package;
2. **initial measure construction** over mass, redshift, scatter, and optional
   orbital variables;
3. **state evolution** through mass, structure, orbit, and survival operators;
4. **catalog and observable operations** over the transported measure.

Backend selection is orthogonal to these layers. A physical model should not
need to be rewritten when switching cosmology or unit backends.

## 4. Backend policy

Backend choices must be explicit, immutable during a calculation, recorded in
metadata, and included in cache keys and regression records.

### 4.1 Cosmology backends

ITAMAE will support at least:

1. **native backend**
   - NumPy/SciPy implementation;
   - reproduces current SASHIMI formulae and conventions;
   - minimal dependencies and low overhead;
   - suitable for regression compatibility and large batch calculations.

2. **Colossus backend**
   - wraps `colossus.cosmology` and relevant halo utilities;
   - exposes established cosmologies and Colossus growth/power-spectrum tools;
   - supports independent validation and interoperability;
   - must not mutate global Colossus cosmology state invisibly.

```python
class CosmologyBackend(Protocol):
    @property
    def identifier(self) -> str: ...

    def H(self, z): ...
    def rho_crit(self, z): ...
    def rho_m(self, z): ...
    def omega_m(self, z): ...
    def growth_factor(self, z): ...
    def collapse_threshold(self, z): ...
    def cosmic_time(self, z): ...
    def lookback_time(self, z): ...
```

Proposed implementations:

```text
itamae.cosmology.NativeFlatLCDM
itamae.cosmology.ColossusCosmology
```

Adapters must state conventions for physical/comoving quantities, factors of
`h`, density definitions, growth normalization, scalar/array behavior, and
supported redshift ranges.

### 4.2 Unit backends

ITAMAE will support at least:

1. **native unit backend**
   - plain floating-point NumPy arrays;
   - one documented canonical internal unit system;
   - optimized for large quadrature catalogs and repeated evolution calls.

2. **Astropy unit backend**
   - accepts and returns `astropy.units.Quantity` where requested;
   - performs dimensional validation and explicit conversion;
   - supports user-facing, analysis-facing, and validation workflows.

```python
class UnitBackend(Protocol):
    @property
    def identifier(self) -> str: ...

    def to_internal(self, value, physical_type: str): ...
    def from_internal(self, value, unit): ...
    def validate(self, value, physical_type: str): ...
```

Proposed implementations:

```text
itamae.units.NativeUnits
itamae.units.AstropyUnits
```

Likely canonical internal units are:

```text
mass                    : Msun
length                  : Mpc
velocity                : km / s
time                    : Gyr
cross section per mass  : cm^2 / g
```

The exact schema must be finalized after regression fixtures exist for C, W,
SI, and F.

Astropy support is a real public backend. However, ITAMAE should not force
`Quantity` objects through every large internal batch. The preferred path is:

```text
Quantity input
  -> dimensional validation
  -> conversion to canonical floating arrays
  -> high-performance calculation
  -> optional Quantity output
```

### 4.3 Backend configuration

```python
config = BackendConfig(
    cosmology="native",  # or "colossus"
    units="native",      # or "astropy"
    array="numpy",
)
```

Objects retain backend identifiers in metadata. Existing objects and cached
results must not change when external global state changes.

### 4.4 Array and numerical backends

The initial numerical backend is NumPy/SciPy. JAX or other accelerated or
differentiable backends are future work. The first API should avoid gratuitous
NumPy-only assumptions, but backend generality must not delay a
regression-equivalent implementation.

## 5. Core data model

### 5.1 `HostState`

```python
@dataclass(frozen=True)
class HostState:
    redshift: Array
    time: Array
    m200: Array
    mvir: Array
    r200: Array
    rvir: Array
    concentration: Array
    metadata: Mapping[str, Any]
```

Host density, enclosed mass, potential, circular velocity, and local dynamical
time are exposed through a `HostPotential` interface.

### 5.2 `AccretionBatch`

```python
@dataclass(frozen=True)
class AccretionBatch:
    m200_acc: Array
    mvir_acc: Array
    z_acc: Array
    concentration_acc: Array
    weight_base: Array
    weight_concentration: Array
    metadata: Mapping[str, Any]
```

Optional fields include host-history and orbital-infall node identifiers. All
arrays share one leading batch shape.

### 5.3 `SubhaloState`

```python
@dataclass
class SubhaloState:
    m_bound: Array
    profile: ProfileParameters
    alive: Array
    flags: Array
    extra: Mapping[str, Array]
```

Profile parameters are named and extensible, allowing NFW, SIDM-cored, and
FDM-soliton structures.

### 5.4 `OrbitalState`

```python
@dataclass
class OrbitalState:
    energy: Array | None = None
    angular_momentum: Array | None = None
    radius: Array | None = None
    radial_velocity: Array | None = None
    tangential_velocity: Array | None = None
    pericenter: Array | None = None
    apocenter: Array | None = None
    phase: Array | None = None
```

Spatial information is separate from internal halo structure. A radial-PDF
model may populate radius nodes and spatial weights; an orbit-averaged model may
use `(E, L)`; an orbit-sampling model may populate instantaneous phase space.

### 5.5 `WeightedSubhaloCatalog`

The catalog retains independent weight factors:

```text
weight_base
weight_host_history
weight_concentration
weight_orbit
weight_survival
weight_final
```

It supports selection, concatenation, weighted reductions, serialization,
stochastic realization, optional Quantity export, and model-specific columns.

## 6. Planned package structure

```text
itamae/
  pyproject.toml
  README.md
  PLAN.md
  LICENSE
  CHANGELOG.md
  CITATION.cff
  uv.lock
  .github/
    workflows/
      test.yml
      release.yml
  src/
    itamae/
      __init__.py
      py.typed
      backends/
      units/
      cosmology/
      types/
      numerics/
      halo/
      protocols/
      evolution/
      measure/
      spatial/
      adapters/
      testing/
  tests/
    unit/
    integration/
    regression/
    property/
  docs/
```

Empty modules should not be created prematurely. Every added module must have a
clear purpose, tests where applicable, and module-level documentation.

## 7. Packaging, PyPI, and uv policy

ITAMAE is intended for publication on PyPI under the distribution name
`itamae`, subject to confirming name availability before the first release.
The import package is also `itamae`.

### 7.1 `pyproject.toml`

The repository must use a standards-compliant `pyproject.toml` as the single
source of packaging metadata, dependencies, development tooling, test settings,
lint settings, and build configuration. A lightweight PEP 517 backend such as
Hatchling is preferred initially.

Illustrative configuration:

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "itamae"
dynamic = ["version"]
description = "Integrated Toolkit for Analytical Merger-tree And Evolution"
readme = "README.md"
requires-python = ">=3.11"
license = { file = "LICENSE" }
authors = [{ name = "Shunichi Horigome" }]
dependencies = [
  "numpy>=1.26",
  "scipy>=1.11",
]

[project.optional-dependencies]
astropy = ["astropy>=6"]
colossus = ["colossus>=1.3"]
full = ["astropy>=6", "colossus>=1.3"]

[dependency-groups]
dev = [
  "pytest>=8",
  "pytest-cov>=5",
  "hypothesis>=6",
  "ruff>=0.5",
  "mypy>=1.10",
  "numpydoc>=1.7",
  "build>=1.2",
  "twine>=5",
]
docs = ["sphinx>=7", "myst-parser>=3", "furo>=2024.5.6"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "--strict-markers --strict-config"

[tool.hatch.version]
path = "src/itamae/__init__.py"

[tool.hatch.build.targets.wheel]
packages = ["src/itamae"]
```

Exact minimum dependency versions should be finalized after testing the current
SASHIMI repositories. Published metadata should use sensible ranges rather than
pinning the complete development lock file.

### 7.2 Installation with uv

End-user installation:

```bash
uv add itamae
uv add "itamae[astropy]"
uv add "itamae[colossus]"
uv add "itamae[full]"
```

Local development:

```bash
git clone https://github.com/gomeshun/itamae.git
cd itamae
uv sync --all-extras --dev
uv run pytest
```

The repository commits `uv.lock` for reproducible development and CI. The lock
file does not replace appropriate dependency ranges in `pyproject.toml`.

### 7.3 Build and release checks

Before every release:

```bash
uv sync --all-extras --dev
uv run ruff check .
uv run mypy src/itamae
uv run pytest --cov=itamae --cov-report=term-missing
uv run python -m numpydoc src/itamae
uv build
uv run twine check dist/*
```

CI must install the built wheel into a clean environment and run a smoke test.
Tests must cover the minimal installation and supported optional backends.

### 7.4 PyPI publication

- Publish initial candidates to TestPyPI first.
- Use GitHub Actions with PyPI trusted publishing where possible.
- Build source distributions and wheels from tagged commits.
- Publish only artifacts that passed the full release workflow.
- Verify installation of uploaded artifacts in a fresh uv environment.
- Update `CHANGELOG.md` for every public release.

## 8. Documentation and code-comment policy

Readable scientific software requires explanation of both software behavior and
physical conventions. Documentation is part of the implementation, not an
optional cleanup task.

### 8.1 NumPy-style docstrings

All public modules, classes, functions, methods, protocols, and dataclasses must
have sufficiently detailed NumPy-style docstrings. Important private helpers
should also be documented when their behavior, numerical assumptions, shapes,
or units are not obvious.

Docstrings should include the relevant NumPy documentation sections:

```text
Summary
Extended Summary, where useful
Parameters
Returns
Yields
Raises
Warns
Other Parameters
Attributes
Notes
References
Examples
See Also
```

Not every section is required for every object, but `Parameters`, `Returns`,
`Raises`, and `Notes` must be included whenever applicable.

Documentation must state:

- accepted scalar and array shapes;
- broadcasting rules;
- expected units and returned units;
- physical versus comoving conventions;
- factors of `h` in masses and lengths;
- mass-definition conventions;
- normalization conventions;
- valid parameter and redshift ranges;
- behavior outside interpolation domains;
- numerical algorithms and convergence assumptions;
- references for nontrivial physical formulae.

Example:

```python
def enclosed_mass(self, radius):
    """Evaluate the mass enclosed by the profile.

    Parameters
    ----------
    radius : float or numpy.ndarray or astropy.units.Quantity
        Physical radius. Plain floating-point input is interpreted in the
        canonical internal length unit. Array input follows NumPy broadcasting.

    Returns
    -------
    mass : float or numpy.ndarray or astropy.units.Quantity
        Enclosed mass with the same broadcast shape as ``radius``. Quantity
        output is returned when requested by the active unit backend.

    Raises
    ------
    ValueError
        If any radius is negative.
    astropy.units.UnitConversionError
        If Quantity input does not have dimensions of length.

    Notes
    -----
    The implementation uses the analytic NFW enclosed-mass function and its
    small-radius series expansion to avoid cancellation.
    """
```

### 8.2 Module and class documentation

Each module must begin with a module-level docstring explaining:

- the module's responsibility;
- what belongs and does not belong in the module;
- the main public objects;
- important numerical or physical conventions;
- relevant references when the module implements literature formulae.

Each class must document its responsibility, state invariants, mutability,
backend interactions, and thread/global-state assumptions. Dataclass fields
whose units or semantics are not obvious must be documented in the class
`Attributes` section.

### 8.3 Inline comments

Inline comments should explain **why** a non-obvious operation is necessary,
not merely restate the code. Detailed comments are required around:

- changes of variables and Jacobians;
- quadrature normalization;
- array reshaping and axis ordering;
- asymptotic expansions and numerical-stability branches;
- root-bracketing choices;
- cache-key construction;
- backend convention conversion;
- physical approximations and empirical calibration boundaries.

Comments must be kept synchronized with the implementation. Obsolete or
misleading comments are bugs and should be corrected in the same pull request.

### 8.4 Documentation enforcement

CI should run documentation-related checks, initially including:

- Ruff docstring rules selected for NumPy-style documentation;
- `numpydoc` validation for public API objects;
- Sphinx documentation build with warnings treated as errors once docs exist;
- doctests for stable executable examples where appropriate.

Documentation checks may begin with a documented temporary allow-list during
the initial migration, but new public APIs must not be added without compliant
docstrings.

## 9. Physical interfaces

Interfaces should use structural typing (`Protocol`) where practical.

```python
class HostHistoryModel(Protocol):
    def m200(self, host_reference, z, cosmology): ...
    def dmvir_dz(self, host_reference, z, cosmology): ...

class VarianceModel(Protocol):
    def sigma(self, mass, z=0.0): ...
    def variance(self, mass, z=0.0): ...
    def dvariance_dmass(self, mass, z=0.0): ...

class AccretionRateModel(Protocol):
    def differential_number(self, m_acc, z_acc, host, variance): ...

class ConcentrationModel(Protocol):
    def median(self, m200, z, cosmology): ...

class InitialStructureModel(Protocol):
    def assign(self, m200, z, concentration_nodes, context): ...

class MassLossLaw(Protocol):
    def rhs(self, state, host_state, orbital_state=None): ...

class ProfileEvolutionModel(Protocol):
    def evolve(self, initial_profile, mass_history, context): ...

class SurvivalModel(Protocol):
    def evaluate(self, state, context): ...
```

Numerical solvers and physical right-hand sides remain separate. The current
SASHIMI perturbative and Shanks methods belong to generic solver machinery;
calibrated stripping coefficients remain in the relevant SASHIMI package.

## 10. Spatial and orbital design

Spatial support is designed from the beginning, even though implementation
follows the non-spatial core.

### Level A: conditional radial measure

```python
class RadialMeasureModel(Protocol):
    def nodes(self, accretion_batch, host, target_redshift): ...
```

This supports `P(q | z_acc, ...)` without tracking individual orbits.

### Level B: local-environment evolution

A `LocalEnvironment` evaluates:

```text
rho_host(r, z)
M_host(<r, z)
Phi(r, z)
Vcirc(r, z)
tdyn_local(r, z)
tidal-radius or tidal-tensor inputs
```

Mass-loss and survival models may consume this context. Radius-dependent fitting
laws remain physical prescriptions in SASHIMI-C or a future spatial package.

### Level C: orbit-averaged transport

The research backend may represent the ensemble in `(E, L)` or action space.
ITAMAE provides turning points, radial periods, radial kernels, infall-variable
transformations, cached phase tables, and conservative transport machinery.
Infall distributions, dynamical-friction laws, Coulomb logarithms, and
disruption models remain model-supplied.

Every spatial implementation must satisfy normalization tests:

```text
integral dq P(q | x) = 1
integral dV p_V(r | E, L) = 1
sum child spatial weights = parent weight
integral dV n(r) = total surviving weight
```

A lone `radius` column is insufficient; its representation, epoch, and weight
semantics must be recorded.

## 11. Variant integration policy

### SASHIMI-C

First migration target. It retains its public API, model composition, calibrated
coefficients, disruption convention, observables, and scientific examples.
ITAMAE replaces only shared mechanisms.

### SASHIMI-W

Adapter-first migration because of its cgs units, WMAP7 setup, global file
loading, procedural initialization, sharp-k variance, and distinct
concentration calculation. Native and Astropy unit backends are especially
important for reproducing and validating this migration.

### SASHIMI-SI

Uses common accretion and tidal-history infrastructure but retains SIDM physics.
The catalog supports multiple named state views such as `cdm_reference` and
`sidm` with shared initial nodes and base weights.

### SASHIMI-F

Uses composable power-spectrum and variance components. FDM population
suppression and FDM core-halo structure remain separate model components. Its
existing Colossus usage is an initial integration target for the Colossus
backend.

## 12. Test strategy and bug prevention

Every implemented behavior must have tests proportional to its scientific and
software risk. A feature is not complete when it merely runs for one example;
it is complete when its contract and important failure modes are tested.

### 12.1 Test categories

The repository should distinguish:

1. **unit tests**
   - individual functions, classes, validation, and edge cases;
   - small and deterministic;
   - run on every pull request.

2. **integration tests**
   - interactions among backends, profiles, solvers, and catalogs;
   - package installation and import behavior;
   - run on every pull request where practical.

3. **golden regression tests**
   - compact reference outputs from SASHIMI-C/W/SI/F;
   - protect scientific results during extraction and refactoring;
   - tolerances and fixture provenance must be documented.

4. **invariant and property-based tests**
   - mathematical identities, monotonicity, normalization, and dimensional
     consistency;
   - use Hypothesis where random generation adds meaningful coverage.

5. **convergence tests**
   - verify stability under grid, quadrature, interpolation, and solver
     resolution changes;
   - expensive convergence tests may run on a scheduled workflow.

6. **failure-mode tests**
   - invalid units, negative masses/radii, unbracketed roots, unsupported
     redshifts, missing optional dependencies, invalid backend state, and
     out-of-domain interpolation.

7. **serialization and compatibility tests**
   - round-trip catalog serialization;
   - schema-version handling;
   - backward-compatibility adapters where promised.

### 12.2 Required invariants

Tests must cover, where applicable:

- finite, nonnegative quadrature and catalog weights;
- agreement between integrated weights and total expected abundance;
- normalization of concentration, host-history, and spatial measures;
- `m_bound <= m_acc` for monotonic stripping models;
- consistency between profile parameters and enclosed mass;
- invertibility of supported mass-definition transformations within tolerance;
- scalar/array and broadcasting equivalence;
- serial/chunked/parallel equivalence;
- native/Colossus agreement for shared cosmological quantities;
- native-float/Astropy-Quantity equivalence;
- recovery of the global model after normalized spatial marginalization;
- deterministic results for fixed random seeds;
- finite results or explicit flags rather than silent NaNs.

### 12.3 Tests accompanying changes

Every bug fix must include a regression test that fails before the fix and
passes after it. Every new public function or behavioral branch must include
unit tests. New numerical algorithms must include at least one comparison to an
independent calculation, analytic limit, or trusted legacy output.

Pull requests that intentionally omit tests must explain why the change is not
testable or why existing tests are sufficient. Such exceptions should be rare.

### 12.4 Coverage policy

Coverage percentage is a diagnostic, not a substitute for scientific test
quality. Nevertheless:

- coverage must be reported in CI;
- coverage must not decrease without an explicit explanation;
- critical modules such as units, backend conversion, weight construction,
  profile mass functions, and solver logic require branch coverage;
- generated code and defensive import guards may be excluded only explicitly.

An initial repository-wide threshold may start modestly during migration, then
increase as legacy code is replaced. New modules should target high line and
branch coverage from their first merge.

## 13. GitHub Actions continuous integration

Automated tests must be implemented in `.github/workflows/test.yml` and run on:

```yaml
on:
  pull_request:
  push:
    branches: [main]
  workflow_dispatch:
```

A scheduled workflow may run slower convergence and full regression suites.

### 13.1 Required CI jobs

The test workflow should contain separate or matrix-based jobs for:

- lint and formatting checks with Ruff;
- static typing with mypy;
- NumPy-style docstring validation;
- unit and integration tests;
- minimal dependency installation;
- Astropy extra installation;
- Colossus extra installation;
- full extra installation;
- supported Python-version matrix;
- wheel and sdist build;
- installation of the built wheel in a clean uv environment;
- import and minimal numerical smoke test.

Illustrative matrix:

```text
Python 3.11 / minimal
Python 3.12 / minimal
Python 3.13 / minimal, when dependencies support it
Python 3.12 / astropy
Python 3.12 / colossus
Python 3.12 / full
```

Optional-backend jobs may be marked separately only when an upstream dependency
is temporarily unavailable. Core tests must never be silently ignored.

### 13.2 Workflow commands

The primary commands should be reproducible locally:

```bash
uv sync --all-extras --dev
uv run ruff check .
uv run ruff format --check .
uv run mypy src/itamae
uv run pytest --cov=itamae --cov-branch --cov-report=term-missing
uv build
uv run twine check dist/*
```

CI-specific behavior should be minimal. A developer should be able to reproduce
a failure with the same uv commands.

### 13.3 Branch protection

Once the initial workflow is stable, the main branch should require successful
CI checks before merge. Recommended required checks are:

```text
lint
typecheck
tests-minimal
tests-astropy
tests-colossus
build-and-smoke-test
```

Direct release publication must depend on successful completion of the full test
workflow. A failing or skipped required check must prevent PyPI publication.

### 13.4 Scheduled and release workflows

A scheduled workflow should periodically run:

- full golden regressions;
- convergence tests;
- latest-compatible dependency resolution;
- optional lower-bound dependency tests;
- documentation build;
- cache and serialization compatibility checks.

The release workflow should:

1. verify that the tag matches the package version;
2. run or require the complete CI suite;
3. build wheel and sdist once;
4. test those exact artifacts;
5. publish the same artifacts through trusted publishing.

## 14. Caching and reproducibility

Cache keys include physical parameters, backend identifiers, cosmology,
canonical-unit schema version, power-spectrum source, numerical resolution,
package version, and serialization schema version.

Global backend state must not determine cached results. Randomized APIs must
accept an explicit `numpy.random.Generator` or seed and record sufficient
metadata for reproducibility.

## 15. Phased roadmap

### Phase 0: packaging, CI, documentation rules, and backend contracts

- add `pyproject.toml`, `src/itamae`, `tests`, and `.github/workflows`;
- configure uv and commit `uv.lock`;
- define supported Python/NumPy/SciPy versions;
- define `CosmologyBackend`, `UnitBackend`, and immutable `BackendConfig`;
- configure Ruff, mypy, pytest, coverage, and numpydoc;
- establish NumPy-style docstring and module-comment requirements;
- add initial CI matrix and clean-wheel smoke test;
- document canonical units and mass definitions;
- add golden-output scripts to C, W, SI, and F;
- make no scientific changes.

### Phase 1: native and Astropy unit support

- implement `NativeUnits` and `AstropyUnits`;
- add dimensional-validation and failure-mode tests;
- implement legacy-unit adapters for C/W/SI/F;
- retain plain internal arrays by default;
- document every public unit conversion with NumPy-style docstrings.

### Phase 2: native and Colossus cosmology support

- implement native flat-LCDM backend;
- implement Colossus adapter without hidden global-state changes;
- add backend-equivalence, convention, and global-state tests;
- pass cosmology explicitly into host-history and halo utilities.

### Phase 3: types, numerics, and halo primitives

- implement state/catalog dataclasses;
- implement grids, quadrature, integration, and interpolation helpers;
- implement robust NFW and mass-definition utilities;
- add analytic-limit, inverse-consistency, broadcasting, and property tests;
- provide unit-aware public wrappers.

### Phase 4: generic evolution solver

- extract ODE and perturbative runners;
- implement optional Shanks acceleration;
- support history grids and chunked execution;
- add comparisons against direct ODE integration and legacy results;
- keep calibrated mass-loss laws in SASHIMI packages.

### Phase 5: initial measure and common catalog builder

- construct accretion batches and independent weights;
- implement concentration quadrature;
- implement deterministic catalogs and stochastic realizations;
- test normalization, shape semantics, serialization, and reproducibility;
- migrate C first, then SI.

### Phase 6: power spectrum and variance protocols

- implement tabulated spectra, transfer functions, top-hat and sharp-k windows;
- implement variance integration, derivatives, interpolation, and safe caching;
- integrate W and F through adapters;
- test asymptotic behavior, cutoff behavior, derivatives, and legacy grids;
- permit native or Colossus-backed implementations where appropriate.

### Phase 7: spatial Level A and Level B

- implement radial measures and explicit spatial weights;
- implement host-potential and local-environment interfaces;
- implement normalized radial PDFs and radial observables;
- add normalization and global-model recovery tests;
- reimplement useful `r-dependent` functionality without embedding physical
  fitting laws in the catalog builder.

### Phase 8: SIDM/FDM structures and Level C research backend

- support multiple named state views;
- support cored/SIDM and soliton/FDM profile schemas;
- provide downstream dSph adapters;
- prototype orbit-averaged transport after radial validation;
- add conservation, normalization, and limiting-case tests.

### Phase 9: first public release

- freeze the initial public API;
- complete NumPy-style API documentation and examples;
- select and add the license;
- add citation and changelog files;
- require all branch-protection checks to pass;
- test wheel/sdist installation across supported Python versions;
- publish a release candidate to TestPyPI;
- verify `uv add itamae` and optional extras in clean environments;
- publish the first PyPI prerelease through trusted publishing.

## 16. Initial public API target

```python
from itamae.backends import BackendConfig
from itamae.cosmology import NativeFlatLCDM, ColossusCosmology
from itamae.units import NativeUnits, AstropyUnits
from itamae.types import AccretionBatch, WeightedSubhaloCatalog
from itamae.halo import NFWProfile, convert_mass_definition
from itamae.evolution import PerturbativeEvolutionSolver
```

```python
backend = BackendConfig(
    cosmology=ColossusCosmology(name="planck18"),
    units=AstropyUnits(),
)
```

A SASHIMI package assembles the physical model:

```python
from sashimi_c import SashimiCDM

model = SashimiCDM(backend=backend)
catalog = model.generate_catalog(host_mass=1.0e12, redshift=0.0)
```

or with Quantity input:

```python
import astropy.units as u

catalog = model.generate_catalog(
    host_mass=1.0e12 * u.Msun,
    redshift=0.0,
)
```

## 17. Immediate next tasks

The first implementation pull request should contain only:

1. `pyproject.toml`, `src` layout, uv configuration, and `uv.lock`;
2. `.github/workflows/test.yml` with lint, typing, documentation, test, backend,
   build, and clean-install jobs;
3. pytest, coverage, Ruff, mypy, Hypothesis, and numpydoc configuration;
4. package metadata and dynamic version setup;
5. backend protocols and immutable `BackendConfig`;
6. canonical-unit documentation;
7. `NativeUnits` and a minimal `AstropyUnits` adapter;
8. native flat-LCDM and a minimal Colossus adapter;
9. unit, failure-mode, and backend-equivalence tests;
10. NumPy-style module, class, and function docstrings;
11. wheel/sdist build and clean-install smoke tests.

It should not yet move EPS, concentration, stripping coefficients, or
SIDM/WDM/FDM physics into ITAMAE.

## 18. Definition of success

The initial refactor succeeds when:

- each SASHIMI variant still owns its physical assumptions;
- common numerical code is no longer copied among repositories;
- all variants remain regression-equivalent;
- every implemented feature has appropriate automated tests;
- every bug fix includes a regression test;
- required GitHub Actions checks run automatically on pull requests and `main`;
- failing tests prevent merge and release publication;
- public modules, classes, functions, and methods have detailed NumPy-style
  documentation;
- users can select native or Colossus cosmology explicitly;
- users can use native floats or Astropy Quantity interfaces explicitly;
- backend choices are reproducible and included in metadata and cache keys;
- weighted catalogs share one schema;
- WDM/FDM variance models are exchangeable without changing EPS integration;
- a radial measure can be added without rewriting the non-spatial pipeline;
- future orbit-averaged work reuses the same abstractions;
- PyPI artifacts can be installed with `uv add itamae`;
- minimal, Astropy, Colossus, and full installations pass clean-environment
  tests.
