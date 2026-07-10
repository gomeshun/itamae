# ITAMAE implementation plan

**ITAMAE: Integrated Toolkit for Analytical Merger-tree And Evolution**

This document defines the implementation policy for ITAMAE after reviewing the
current SASHIMI-C, SASHIMI-W, SASHIMI-SI, and SASHIMI-F implementations, the
FDM profile implementation in `dsph_fuzzy`, and the spatial development in the
`r-dependent` branch of SASHIMI-C.

The goal is to remove duplicated computational machinery without moving the
scientific identity of each SASHIMI variant into ITAMAE.

## 1. Executive decision

ITAMAE owns the common **computational language** of the SASHIMI family:

1. typed state and catalog objects;
2. backend-independent cosmology and unit interfaces;
3. spherical-halo and profile primitives;
4. grids, quadrature, interpolation, integration, caching, and batch execution;
5. generic evolution solvers that operate on model-supplied equations;
6. deterministic weighted measures and stochastic realizations;
7. optional radial and phase-space infrastructure.

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

ITAMAE will support explicit, exchangeable backends. Backend choices must be
stored in configuration metadata and included in cache keys and regression
records.

### 4.1 Cosmology backends

ITAMAE will support at least two cosmology backends:

1. **native backend**
   - NumPy/SciPy implementation;
   - reproduces the current SASHIMI formulae and conventions;
   - minimal dependencies and low overhead;
   - suitable for regression compatibility and large batch calculations.

2. **Colossus backend**
   - wraps `colossus.cosmology` and relevant halo utilities;
   - allows established cosmologies and Colossus power-spectrum/growth tools;
   - useful for independent validation and interoperability;
   - must not mutate global Colossus cosmology state invisibly.

The public interface is backend-independent:

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

A backend adapter must define conventions explicitly, including:

- physical versus comoving distances;
- masses with or without factors of `h`;
- critical-density versus mean-density definitions;
- normalization of the growth factor;
- scalar and array behavior;
- supported redshift range.

Colossus is an optional dependency, exposed through an installation extra such
as:

```bash
pip install itamae[colossus]
```

The native and Colossus backends must be compared on a common test matrix. They
need not be bitwise identical, but differences must be understood and bounded.

### 4.2 Unit backends

ITAMAE will support at least two unit modes:

1. **native unit backend**
   - plain floating-point NumPy arrays;
   - one documented canonical internal unit system;
   - optimized for large quadrature catalogs and repeated evolution calls.

2. **Astropy unit backend**
   - accepts and returns `astropy.units.Quantity` where requested;
   - performs dimensional validation and explicit conversion;
   - supports user-facing, analysis-facing, and validation workflows.

The backend contract is conceptually:

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

The canonical internal units should be documented centrally. A likely initial
choice is:

```text
mass       : Msun
length     : Mpc
velocity   : km / s
time       : Gyr
cross section per mass : cm^2 / g
```

The exact choice should be finalized only after regression tests against C, W,
SI, and F are prepared.

Astropy support must be a real public backend, not merely a test helper.
Nevertheless, ITAMAE should not force `Quantity` objects through every large
internal batch. The recommended execution path is:

```text
Quantity input
  -> dimensional validation
  -> conversion to canonical internal floating arrays
  -> high-performance calculation
  -> optional Quantity output
```

This gives users safe unit-aware interfaces without imposing Quantity overhead
on every inner-loop operation.

Astropy is an optional dependency, exposed through an installation extra such
as:

```bash
pip install itamae[astropy]
```

An aggregate extra may also be provided:

```bash
pip install itamae[full]
```

### 4.3 Backend configuration

Backend selection must be explicit and immutable during a calculation:

```python
config = BackendConfig(
    cosmology="native",       # or "colossus"
    units="native",           # or "astropy"
    array="numpy",
)
```

Objects constructed under one backend configuration should retain the backend
identifier in metadata. Changing a global Colossus cosmology or changing unit
conventions must not silently alter an existing model or cached result.

### 4.4 Array and numerical backends

The initial array/numerical backend remains NumPy/SciPy. JAX or other
accelerated/differentiable backends are future work. The first API should avoid
unnecessary NumPy-only assumptions where a small abstraction is inexpensive,
but backend generality must not delay regression-equivalent implementation.

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

Metadata records the cosmology backend, unit backend, mass definitions, and
physical model identifiers.

Host density, enclosed mass, potential, circular velocity, and local dynamical
time are provided through a `HostPotential` interface.

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
arrays use a shared leading batch shape.

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
src/itamae/
  __init__.py

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
```

This is a target layout; empty modules should not be created prematurely.

## 7. Physical interfaces

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

## 8. Spatial and orbital design

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

## 9. Variant integration policy

### SASHIMI-C

First migration target. It retains its public API, model composition, calibrated
coefficients, disruption convention, observables, and scientific examples.
ITAMAE replaces only shared mechanisms.

### SASHIMI-W

Adapter-first migration because of its cgs units, WMAP7 setup, global file
loading, procedural initialization, sharp-k variance, and distinct
concentration calculation. The native and Astropy unit backends are especially
important for reproducing and validating this migration.

### SASHIMI-SI

Uses common accretion and tidal-history infrastructure but retains SIDM physics.
The catalog supports multiple named state views such as `cdm_reference` and
`sidm` with shared initial nodes and base weights.

### SASHIMI-F

Uses composable power-spectrum and variance components. FDM population
suppression and FDM core-halo structure remain separate model components.
Its existing Colossus usage provides an initial integration target for the
Colossus cosmology backend.

## 10. Testing policy

### Golden regression data

Each SASHIMI repository should produce compact reference outputs spanning:

- multiple host masses;
- zero and nonzero target redshifts;
- low and standard resolutions;
- concentration scatter on/off;
- evolved and unevolved profiles;
- representative WDM, SIDM, and FDM parameters;
- at least one radial configuration.

### Invariant tests

Required invariants include:

- finite, nonnegative quadrature weights;
- integrated-weight agreement with expected abundance;
- normalized concentration and spatial measures;
- profile mass consistency;
- scalar/batch equivalence;
- serial/parallel equivalence;
- recovery of the global model after normalized spatial marginalization.

### Backend-equivalence tests

The test suite must compare:

- native cosmology versus Colossus for shared cosmological quantities;
- native floats versus Astropy Quantity inputs and outputs;
- native and Astropy unit conversions for all public physical quantities;
- catalog results across backend combinations within documented tolerances;
- legacy C/W/SI/F outputs after conversion to one canonical unit system.

Tests must include deliberate unit mistakes and verify that the Astropy backend
raises clear dimensional errors.

## 11. Caching and reproducibility

Cache keys include:

- physical parameters;
- cosmology parameters and backend identifier;
- unit backend and canonical-unit schema version;
- power-spectrum source and version;
- numerical bounds and resolution;
- package/code version.

Global backend state must not determine cached results. In particular, Colossus
configuration should be isolated or explicitly restored, and existing objects
must retain their construction-time cosmology definition.

## 12. Phased roadmap

### Phase 0: baseline and backend contracts

- add packaging and CI;
- define supported Python/NumPy/SciPy versions;
- define `CosmologyBackend`, `UnitBackend`, and immutable `BackendConfig`;
- document canonical internal units and mass-definition conventions;
- add golden-output scripts to C, W, SI, and F;
- make no scientific changes.

### Phase 1: native and Astropy unit support

- implement `NativeUnits`;
- implement `AstropyUnits` with Quantity input/output conversion;
- add dimensional-validation tests;
- implement legacy-unit adapters for C/W/SI/F;
- ensure internal batch arrays remain plain floating arrays by default.

### Phase 2: native and Colossus cosmology support

- implement native flat-LCDM backend reproducing current formulae;
- implement the Colossus adapter without hidden global-state changes;
- add backend-equivalence tests;
- pass cosmology explicitly into host-history and halo utilities.

### Phase 3: types, numerics, and halo primitives

- implement state/catalog dataclasses;
- implement grids and Gauss-Hermite utilities;
- implement integration/interpolation wrappers;
- implement robust NFW mass and inverse-mass functions;
- implement mass-definition conversions;
- provide unit-aware public wrappers.

### Phase 4: generic evolution solver

- extract ODE and perturbative runners;
- implement optional Shanks acceleration;
- support history grids and chunked execution;
- keep calibrated mass-loss laws in SASHIMI packages.

### Phase 5: initial measure and common catalog builder

- implement accretion-batch construction;
- retain independent weight factors;
- implement concentration quadrature;
- implement deterministic catalogs and stochastic realizations;
- migrate C first, then SI.

### Phase 6: power spectrum and variance protocols

- implement tabulated spectra, transfer functions, top-hat and sharp-k windows;
- implement variance integration, derivatives, interpolation, and safe caching;
- integrate W and F through adapters;
- allow native or Colossus-backed variance implementations where appropriate.

### Phase 7: spatial Level A and Level B

- implement radial measures and explicit spatial weights;
- implement host-potential and local-environment interfaces;
- implement normalized radial PDFs and radial observables;
- reimplement useful `r-dependent` functionality without embedding its fitting
  laws in the catalog builder.

### Phase 8: SIDM/FDM structures and Level C research backend

- support multiple named state views;
- support cored/SIDM and soliton/FDM profile schemas;
- provide downstream dSph adapters;
- prototype orbit-averaged phase-space transport after radial validation.

## 13. Initial public API target

```python
from itamae.backends import BackendConfig
from itamae.cosmology import NativeFlatLCDM, ColossusCosmology
from itamae.units import NativeUnits, AstropyUnits
from itamae.types import AccretionBatch, WeightedSubhaloCatalog
from itamae.halo import NFWProfile, convert_mass_definition
from itamae.evolution import PerturbativeEvolutionSolver
```

Example backend selection:

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
catalog = model.generate_catalog(
    host_mass=1.0e12,  # floats use documented canonical units
    redshift=0.0,
)
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

1. repository packaging and CI;
2. backend protocols and immutable `BackendConfig`;
3. canonical-unit documentation;
4. `NativeUnits` and a minimal `AstropyUnits` adapter;
5. native flat-LCDM and a minimal Colossus adapter;
6. backend-equivalence tests for `H(z)`, `rho_crit(z)`, time, and basic unit
   conversions.

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
- future orbit-averaged work reuses the same host, state, measure, backend, and
  catalog abstractions.
