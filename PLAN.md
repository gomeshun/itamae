# ITAMAE implementation plan

**ITAMAE: Integrated Toolkit for Analytical Merger-tree And Evolution**

This document defines the implementation policy for ITAMAE after reviewing the
current SASHIMI-C, SASHIMI-W, SASHIMI-SI, and SASHIMI-F implementations, the
FDM profile implementation in `dsph_fuzzy`, and the spatial development in the
`r-dependent` branch of SASHIMI-C.

The goal is to remove duplicated computational machinery without moving the
scientific identity of each SASHIMI variant into ITAMAE. ITAMAE is intended to
be distributed as a normal Python package through PyPI and installed and
developed with modern Python tooling, including `uv`.

## 1. Executive decision

ITAMAE owns the common **computational language** of the SASHIMI family:

1. typed state and catalog objects;
2. backend-independent cosmology and unit interfaces;
3. spherical-halo and profile primitives;
4. grids, quadrature, interpolation, integration, caching, and batch execution;
5. generic evolution solvers operating on model-supplied equations;
6. deterministic weighted measures and stochastic realizations;
7. optional radial and phase-space infrastructure;
8. stable packaging, versioning, and distribution infrastructure.

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
   - reproduces the current SASHIMI formulae and conventions;
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
  src/
    itamae/
      __init__.py
      py.typed
      backends/
        config.py
        registry.py
      units/
        base.py
        native.py
        astropy.py
      cosmology/
        base.py
        native.py
        colossus.py
      types/
        arrays.py
        state.py
        catalog.py
        flags.py
      numerics/
        grids.py
        quadrature.py
        integration.py
        interpolation.py
        root_finding.py
        batching.py
        cache.py
      halo/
        mass_definitions.py
        profiles.py
        nfw.py
        truncated_nfw.py
        potential.py
      protocols/
        variance.py
        host_history.py
        accretion.py
        concentration.py
        mass_loss.py
        profile_evolution.py
        survival.py
        infall.py
        orbit.py
        spatial.py
      evolution/
        ode.py
        perturbative.py
        operators.py
        runner.py
      measure/
        builder.py
        weights.py
        sampling.py
      spatial/
        radial.py
        phase_space.py
        orbit_averaging.py
        kernels.py
      adapters/
        legacy_c.py
        legacy_w.py
        legacy_si.py
        legacy_f.py
      testing/
        regression.py
        convergence.py
        backend_equivalence.py
  tests/
  docs/
```

This is a target layout; empty modules should not be created prematurely.

## 7. Packaging, PyPI, and uv policy

ITAMAE is intended for publication on PyPI under the distribution name
`itamae`, subject to confirming name availability before the first release.
The import package is also `itamae`.

### 7.1 `pyproject.toml`

The repository must use a standards-compliant `pyproject.toml` as the single
source of packaging metadata and dependency declarations. A lightweight PEP 517
backend such as Hatchling is preferred initially.

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
authors = [
  { name = "Shunichi Horigome" },
]
dependencies = [
  "numpy>=1.26",
  "scipy>=1.11",
]

[project.optional-dependencies]
astropy = ["astropy>=6"]
colossus = ["colossus>=1.3"]
full = [
  "astropy>=6",
  "colossus>=1.3",
]
dev = [
  "pytest>=8",
  "pytest-cov>=5",
  "ruff>=0.5",
  "mypy>=1.10",
  "build>=1.2",
  "twine>=5",
]
docs = [
  "sphinx>=7",
  "myst-parser>=3",
  "furo>=2024.5.6",
]

[project.urls]
Repository = "https://github.com/gomeshun/itamae"
Issues = "https://github.com/gomeshun/itamae/issues"

[tool.hatch.version]
path = "src/itamae/__init__.py"

[tool.hatch.build.targets.wheel]
packages = ["src/itamae"]
```

Exact minimum dependency versions should be finalized after testing the current
SASHIMI repositories. Loose but bounded compatibility is preferred over
unnecessarily strict pins in published package metadata.

### 7.2 Versioning

- Use semantic versioning where practical.
- Begin with development releases such as `0.1.0a1` or `0.1.0.devN`.
- Keep one authoritative version source.
- Include the package version and schema version in serialized metadata and
  cache keys.
- Tag release commits as `vX.Y.Z`.

### 7.3 Installation with uv

End-user installation from PyPI:

```bash
uv add itamae
```

With optional backends:

```bash
uv add "itamae[astropy]"
uv add "itamae[colossus]"
uv add "itamae[full]"
```

One-off isolated execution may use:

```bash
uvx --from itamae python -c "import itamae; print(itamae.__version__)"
```

For local development:

```bash
git clone https://github.com/gomeshun/itamae.git
cd itamae
uv sync --all-extras --dev
uv run pytest
```

The repository should commit `uv.lock` for reproducible development and CI.
The lock file is not a substitute for sensible dependency ranges in
`pyproject.toml`; PyPI users resolve from published metadata.

### 7.4 Dependency groups

Where supported by the selected uv/project configuration, development tools may
be declared through dependency groups rather than exposing all tooling as a
runtime-style optional extra. The intended groups are:

```text
dev     : pytest, coverage, ruff, mypy, build, twine
docs    : Sphinx/MyST/Furo
bench   : benchmark and profiling tools
```

Runtime optional extras remain user-facing:

```text
astropy
colossus
full
```

### 7.5 Build and release checks

Before every release:

```bash
uv sync --all-extras --dev
uv run pytest
uv run ruff check .
uv run mypy src/itamae
uv build
uv run twine check dist/*
```

CI should also install the built wheel into a clean environment and run a smoke
test. Tests must cover the minimal installation and each supported optional
backend combination.

Recommended test matrix:

```text
minimal: NumPy/SciPy only
astropy: minimal + Astropy
colossus: minimal + Colossus
full: all supported optional backends
```

### 7.6 PyPI publication

- Publish initial candidates to TestPyPI first.
- Use GitHub Actions with PyPI trusted publishing rather than storing a long-lived
  API token where possible.
- Build source distributions and wheels from tagged commits.
- Do not publish directly from an unclean local working tree.
- Verify installation of the uploaded artifact in a fresh environment.
- Add release notes and update `CHANGELOG.md` for each public release.

A future release workflow should be triggered by a GitHub release or version
tag, run the full test matrix, build artifacts once, and publish those exact
artifacts.

### 7.7 Public package quality requirements

Before the first PyPI release, ITAMAE should have:

- a selected open-source license;
- complete project metadata and classifiers;
- `README.md` installation and minimal usage examples;
- `CHANGELOG.md`;
- `CITATION.cff` and citation guidance;
- `py.typed` if public type annotations are supported;
- no package import-time data-file lookup relative to the current directory;
- packaged data declared explicitly;
- clear optional-dependency error messages;
- API documentation for all public objects.

## 8. Physical interfaces

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

## 9. Spatial and orbital design

Spatial support is designed from the beginning, even though implementation
follows the non-spatial core.

### Level A: conditional radial measure

```python
class RadialMeasureModel(Protocol):
    def nodes(self, accretion_batch, host, target_redshift):
        # q_nodes and normalized weight_orbit
        ...
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
laws remain physical prescriptions in SASHIMI-C or a future spatial model
package.

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

## 10. Variant integration policy

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

## 11. Testing, caching, and reproducibility

Required regression and invariant coverage includes:

- compact golden outputs for C, W, SI, F, and at least one radial case;
- finite and nonnegative quadrature weights;
- integrated-weight agreement with expected abundance;
- normalized concentration and spatial measures;
- profile mass consistency;
- scalar/batch and serial/parallel equivalence;
- recovery of the global model after normalized spatial marginalization;
- native cosmology versus Colossus comparisons;
- native floats versus Astropy Quantity inputs and outputs;
- dimensional-error tests;
- minimal and optional-dependency installation tests.

Cache keys include physical parameters, backend identifiers, cosmology,
canonical-unit schema version, power-spectrum source, numerical resolution,
package version, and serialization schema version.

## 12. Phased roadmap

### Phase 0: packaging, baseline, and backend contracts

- add `pyproject.toml`, `src/itamae`, tests, and CI;
- configure uv and commit `uv.lock`;
- define supported Python/NumPy/SciPy versions;
- define `CosmologyBackend`, `UnitBackend`, and immutable `BackendConfig`;
- document canonical units and mass definitions;
- add golden-output scripts to C, W, SI, and F;
- configure wheel and source-distribution builds;
- make no scientific changes.

### Phase 1: native and Astropy unit support

- implement `NativeUnits` and `AstropyUnits`;
- add dimensional-validation tests;
- implement legacy-unit adapters for C/W/SI/F;
- retain plain internal arrays by default.

### Phase 2: native and Colossus cosmology support

- implement native flat-LCDM backend;
- implement Colossus adapter without hidden global-state changes;
- add backend-equivalence tests;
- pass cosmology explicitly into host-history and halo utilities.

### Phase 3: types, numerics, and halo primitives

- implement state/catalog dataclasses;
- implement grids, quadrature, integration, and interpolation helpers;
- implement robust NFW and mass-definition utilities;
- provide unit-aware public wrappers.

### Phase 4: generic evolution solver

- extract ODE and perturbative runners;
- implement optional Shanks acceleration;
- support history grids and chunked execution;
- keep calibrated mass-loss laws in SASHIMI packages.

### Phase 5: initial measure and common catalog builder

- construct accretion batches and independent weights;
- implement concentration quadrature;
- implement deterministic catalogs and stochastic realizations;
- migrate C first, then SI.

### Phase 6: power spectrum and variance protocols

- implement tabulated spectra, transfer functions, top-hat and sharp-k windows;
- implement variance integration, derivatives, interpolation, and safe caching;
- integrate W and F through adapters;
- permit native or Colossus-backed implementations where appropriate.

### Phase 7: spatial Level A and Level B

- implement radial measures and explicit spatial weights;
- implement host-potential and local-environment interfaces;
- implement normalized radial PDFs and radial observables;
- reimplement useful `r-dependent` functionality without embedding physical
  fitting laws in the catalog builder.

### Phase 8: SIDM/FDM structures and Level C research backend

- support multiple named state views;
- support cored/SIDM and soliton/FDM profile schemas;
- provide downstream dSph adapters;
- prototype orbit-averaged transport after radial validation.

### Phase 9: first public release

- freeze the initial public API;
- complete documentation and examples;
- select and add the license;
- add citation and changelog files;
- test wheel/sdist installation across supported Python versions;
- publish a release candidate to TestPyPI;
- verify `uv add itamae` and optional extras in clean environments;
- publish the first PyPI prerelease through trusted publishing.

## 13. Initial public API target

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

## 14. Immediate next tasks

The first implementation pull request should contain only:

1. `pyproject.toml`, `src` layout, uv configuration, and CI;
2. package metadata and dynamic version setup;
3. backend protocols and immutable `BackendConfig`;
4. canonical-unit documentation;
5. `NativeUnits` and a minimal `AstropyUnits` adapter;
6. native flat-LCDM and a minimal Colossus adapter;
7. backend-equivalence tests;
8. wheel/sdist build and clean-install smoke tests.

It should not yet move EPS, concentration, stripping coefficients, or
SIDM/WDM/FDM physics into ITAMAE.

## 15. Definition of success

The initial refactor succeeds when:

- each SASHIMI variant still owns its physical assumptions;
- common numerical code is no longer copied among repositories;
- all variants remain regression-equivalent;
- users can select native or Colossus cosmology explicitly;
- users can use native floats or Astropy Quantity interfaces explicitly;
- backend choices are reproducible and included in metadata/cache keys;
- weighted catalogs share one schema;
- WDM/FDM variance models are exchangeable without changing EPS integration;
- a radial measure can be added without rewriting the non-spatial pipeline;
- future orbit-averaged work reuses the same abstractions;
- PyPI artifacts can be installed with `uv add itamae`;
- minimal, Astropy, Colossus, and full installations pass clean-environment
  tests.
