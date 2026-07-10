# ITAMAE implementation plan

**ITAMAE: Integrated Toolkit for Analytical Merger-tree And Evolution**

This document records the initial implementation policy for ITAMAE after
reviewing the current implementations of SASHIMI-C, SASHIMI-W, SASHIMI-SI,
SASHIMI-F, the FDM profile implementation in `dsph_fuzzy`, and the spatial
work in the `r-dependent` branch of SASHIMI-C.

The plan is deliberately conservative. The first goal is not to redesign every
SASHIMI model, but to establish a common computational foundation without
moving the scientific identity of each SASHIMI variant into ITAMAE.

## 1. Executive decision

ITAMAE will own the common **computational language** of the SASHIMI family:

1. typed state and catalog objects;
2. cosmological and halo-profile primitives;
3. grids, quadrature, interpolation, integration, caching, and batch execution;
4. interfaces connecting accretion, structure, evolution, survival, and spatial
   models;
5. generic solvers that operate on model-supplied equations;
6. deterministic weighted measures and stochastic realizations;
7. optional phase-space and radial-distribution infrastructure.

Each SASHIMI variant will continue to own its scientific prescription:

- **SASHIMI-C:** CDM variance/concentration choices, Yang-type accretion model,
  calibrated average mass-loss law, tidal structural response, disruption
  convention, and CDM observables;
- **SASHIMI-W:** WDM transfer function, filtering convention, WDM concentration
  calculation, and WDM-specific calibration;
- **SASHIMI-SI:** effective cross section and gravothermal/profile-evolution
  model;
- **SASHIMI-F:** FDM transfer function, variance/filter choices, and FDM-specific
  structure assignment;
- **downstream dSph projects:** Jeans likelihoods, observational priors, and
  inference workflows.

The boundary is summarized as:

> ITAMAE defines how components communicate and how a calculation is executed.  
> A SASHIMI variant defines the physical components that should be used.

ITAMAE must therefore not contain a default monolithic `SashimiModel` whose
hidden defaults reproduce SASHIMI-C. Such a class would make SASHIMI-C a thin
wrapper and obscure which assumptions define each scientific model.

## 2. Findings from the current implementations

### 2.1 Common pipeline

All four variants follow the same high-level sequence:

1. specify cosmology and a host halo at a reference redshift;
2. construct a grid or measure in accretion mass and accretion redshift;
3. evaluate an EPS-like accretion abundance;
4. integrate over host-history and concentration scatter;
5. assign internal halo structure at accretion;
6. evolve mass and structure after accretion;
7. evaluate survival, disruption, or model-validity conditions;
8. flatten quadrature nodes into a weighted catalog;
9. calculate mass functions, satellite counts, boost factors, priors, or a
   stochastic realization.

This common pipeline, rather than any one physical formula, is the main object
that ITAMAE should support.

### 2.2 Repeated implementations

SASHIMI-C, SASHIMI-SI, and SASHIMI-F contain closely related or copied
implementations of:

- units and constants;
- LCDM background functions;
- growth factor and collapse threshold;
- NFW auxiliary functions and mass-definition conversion;
- host mass-accretion history;
- Yang et al. conditional/accretion functions;
- Gauss-Hermite quadrature over concentration and host-history scatter;
- the calibrated average tidal mass-loss equation;
- ODE and perturbative solutions, including the Shanks-accelerated solution;
- conversion between mass fraction, `Vmax`, `rmax`, `rs`, and `rhos`;
- inversion of the NFW enclosed-mass function to obtain a truncation
  concentration;
- weight normalization and array flattening.

SASHIMI-W implements the same conceptual stages in an older, more procedural
form. It differs strongly in units, cosmological parameters, power-spectrum
input, sharp-k variance, and concentration calculation, but its output and
post-accretion pipeline are recognizably the same.

These differences argue for interfaces and adapters, not for copying the
SASHIMI-C implementations into ITAMAE and forcing all variants to inherit from
them.

### 2.3 Model-specific parts that remain outside the core

The following are not generic utilities even when they currently appear in
several files:

- the selected power-spectrum transfer function and window function;
- the selected host mass-accretion relation;
- the Yang-model variant and its normalization details;
- a concentration-mass-redshift relation;
- calibrated coefficients in the mean stripping law;
- the Peñarrubia-type structural response to bound-mass loss;
- a disruption threshold such as `ct > 0.77`;
- SIDM effective-cross-section and gravothermal evolution;
- FDM soliton/core-halo relations;
- WDM concentration and cutoff calibration.

ITAMAE may define protocols for these objects and may later host clearly named,
opt-in reference implementations. It must not silently choose them.

## 3. Architectural model

The fundamental object is a weighted measure transported through a sequence of
physical maps.

Let the initial variables be represented by `x`, for example

```text
x = (m_acc, z_acc, concentration_node, host_history_node, orbit_node, ...)
```

with quadrature or measure weight `w(x)`. A model evolves these variables into a
state `y` at the target epoch:

```text
y = T_model(x)
```

An expectation value is then

```text
sum_i w_i O(y_i)
```

or, in continuum notation,

```text
integral dx lambda(x) O(T_model(x)).
```

The existing SASHIMI weighted catalog is therefore a deterministic quadrature
representation of a population measure. ITAMAE will preserve this
representation as the canonical backend. Monte Carlo catalogs are generated
from it as a separate operation.

The architecture distinguishes four layers:

1. **configuration and physical models** supplied by a SASHIMI package;
2. **initial measure construction** over mass, redshift, scatter, and optional
   orbital variables;
3. **state evolution** through mass, structure, orbit, and survival operators;
4. **catalog and observable operations** over the transported weighted measure.

## 4. Core data model

The data model must be implemented before the physical modules are migrated.
The current long tuple outputs cannot accommodate SIDM fields, FDM core
properties, or spatial coordinates cleanly.

### 4.1 `HostState`

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
```

Host density, enclosed mass, potential, circular velocity, and local dynamical
time should be obtained through a `HostPotential` interface rather than stored
as unrelated helper functions.

### 4.2 `AccretionBatch`

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

Optional fields should include host-history node identifiers and orbital-infall
nodes. Arrays must use one shared leading batch shape.

### 4.3 `SubhaloState`

```python
@dataclass
class SubhaloState:
    m_bound: Array
    profile: ProfileParameters
    alive: Array
    flags: Array
    extra: Mapping[str, Array]
```

The profile parameters must not be restricted to NFW. The base interface should
support named parameter sets so that SIDM core radii and FDM soliton parameters
can coexist with the common state.

### 4.4 `OrbitalState`

Spatial information is represented separately from internal structure:

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

Not every spatial model must populate every field. A simple radial-PDF model may
only provide a radius node and radial weight. An orbit-averaged model may use
`(E, L)` and derive radial kernels. An orbit-sampling backend may populate
instantaneous phase-space coordinates.

### 4.5 `WeightedSubhaloCatalog`

The final catalog should retain separate weight factors:

```text
weight_base
weight_host_history
weight_concentration
weight_orbit
weight_survival
weight_final
```

This separation is required for diagnostics, reweighting, model comparisons,
and spatial marginalization. Boolean survival must not be irreversibly folded
into the only available weight.

The catalog should support:

- selection without losing metadata;
- concatenation and chunked construction;
- weighted sums and histograms;
- conversion to NumPy structured arrays and pandas DataFrames;
- NPZ/HDF5 or another stable serialization format;
- stochastic Poisson realization using an explicit random generator;
- model-specific extension columns without changing the base schema.

## 5. Planned package structure

```text
src/itamae/
  __init__.py

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

  cosmology/
    base.py
    flat_lcdm.py
    time.py

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
```

This layout is a target, not a requirement that every empty module be created
immediately.

## 6. Interfaces to implement

Interfaces should use structural typing (`Protocol`) where practical. Physical
classes should not be forced into a deep inheritance hierarchy.

### 6.1 Cosmology and host history

```python
class Cosmology(Protocol):
    def H(self, z): ...
    def rho_crit(self, z): ...
    def growth_factor(self, z): ...
    def cosmic_time(self, z): ...

class HostHistoryModel(Protocol):
    def m200(self, host_reference, z): ...
    def dmvir_dz(self, host_reference, z): ...
```

The current SASHIMI-C history may remain in SASHIMI-C while using ITAMAE's
cosmology and host-state containers.

### 6.2 Variance and power spectrum

```python
class VarianceModel(Protocol):
    def sigma(self, mass, z=0.0): ...
    def variance(self, mass, z=0.0): ...
    def dvariance_dmass(self, mass, z=0.0): ...
```

This is essential for unifying the variants:

- C and SI currently use an analytic CDM fit;
- W uses a tabulated spectrum, WDM transfer suppression, and sharp-k filtering;
- F uses a tabulated CDM spectrum, FDM transfer suppression, and selectable
  top-hat/sharp-k backends.

Power-spectrum evaluation, transfer functions, window functions, and variance
integration should be composable objects. File caching must be keyed by all
physical and numerical parameters, not only a particle mass.

### 6.3 Accretion measure

```python
class AccretionRateModel(Protocol):
    def differential_number(self, m_acc, z_acc, host, variance): ...
```

ITAMAE initially provides the quadrature engine, broadcasting rules, and
normalization checks. The concrete Yang-model formulas should remain in the
variant during the first migration. After regression equivalence is
established, identical formulas may be moved into an explicitly named optional
module such as `itamae.reference.yang2011`.

### 6.4 Concentration and initial structure

```python
class ConcentrationModel(Protocol):
    def median(self, m200, z): ...

class InitialStructureModel(Protocol):
    def assign(self, m200, z, concentration_nodes): ...
```

This separation is important because W and F may change the concentration model
independently of the EPS variance. ITAMAE must not assume that a dark-matter
model is fully specified by its transfer function.

### 6.5 Mass and profile evolution

```python
class MassLossLaw(Protocol):
    def rhs(self, state, host_state, orbital_state=None): ...

class ProfileEvolutionModel(Protocol):
    def evolve(self, initial_profile, mass_history, context): ...

class SurvivalModel(Protocol):
    def evaluate(self, state, context): ...
```

The numerical solver and the physical right-hand side must be separate. The
current average stripping law supplies `A`, `zeta`, and dynamical-time factors;
ITAMAE supplies ODE and perturbative integration. The Shanks transformation is
a solver option, not a CDM physical model.

The evolution context should include:

```text
host state
redshift/time grid
optional orbital state
local host density and enclosed mass
model metadata
```

This allows the same mass-loss interface to support both the original global
average law and future radius- or orbit-dependent laws.

## 7. Spatial and orbital design

Spatial support must be designed now even if it is implemented after the
non-spatial core.

### 7.1 Lessons from `r-dependent`

The current development branch includes two related but distinct ideas:

1. a conditional radial distribution `P(q | z_acc, ...)`, with
   `q = r / Rvir` and a minimum radius related to elapsed dynamical time;
2. a radius-dependent correction to the mean mass-loss rate.

It also contains a more general theoretical direction in which infalling
subhalos are injected into an `(E, L)` distribution and mapped to radial number
density by an orbit-averaged kernel. Dynamical friction and disruption become
drift and sink terms in integral-of-motion space.

These should not be collapsed into one `radius` option. ITAMAE should support
three levels.

#### Level A: conditional radial weighting

A spatial model returns quadrature nodes and weights in `q`:

```python
class RadialMeasureModel(Protocol):
    def nodes(self, accretion_batch, host, target_redshift):
        # returns q_nodes and normalized weight_orbit
        ...
```

This reproduces a radial-PDF approach and allows radial observables without
tracking individual orbits.

#### Level B: local-environment evolution

A `LocalEnvironment` object evaluates host quantities at the radius nodes:

```text
rho_host(r, z)
M_host(<r, z)
Phi(r, z)
Vcirc(r, z)
tdyn_local(r, z)
tidal tensor or tidal-radius inputs
```

A mass-loss or survival model may consume this context. This avoids embedding
functions such as `mdot_r_0`, `mdot_r_1`, or `mdot_r_2` in the catalog builder.
The radial correction remains a physical model owned by SASHIMI-C or a future
spatial SASHIMI package.

#### Level C: orbit-averaged phase-space transport

A more advanced backend represents the ensemble in `(E, L)` or another action
space. ITAMAE should provide reusable numerical pieces:

- spherical host potential interface;
- turning-point calculation;
- radial period;
- orbit-averaged radial kernel;
- mapping from infall velocity variables to `(E, L)`;
- quadrature/interpolation tables for phase kernels;
- conservative transport runner for source, drift, diffusion, and sink terms.

The physical infall distribution, Coulomb logarithm, dynamical-friction law,
and disruption prescription remain supplied by the model package.

### 7.2 Spatial normalization requirements

Every spatial backend must satisfy explicit normalization tests:

```text
integral dq P(q | x) = 1
integral dV p_V(r | E, L) = 1
sum spatial weights = parent node weight
integral dV n(r) = total surviving catalog weight
```

The implementation must keep clear whether a function is a density per `q`,
per `r`, per shell, or per volume. Variable names should encode this when
possible, for example `pdf_q`, `pdf_r`, and `number_density_3d`.

### 7.3 Do not store only instantaneous radius

A catalog column named `radius` alone is insufficient because it cannot state
whether the value is:

- a deterministic expected radius;
- a quadrature node in a radial PDF;
- a random orbital-phase realization;
- a pericenter or apocenter;
- a radius at accretion or at the target epoch.

Spatial values must be accompanied by a representation type and relevant
weights or orbital invariants.

## 8. Variant integration policy

### 8.1 SASHIMI-C

SASHIMI-C is the first migration target because it has the clearest current API
and the shared perturbative solver.

SASHIMI-C retains:

- its public `subhalo_properties` and observable API during migration;
- its standard physical model assembly;
- its calibrated coefficients and references;
- its default disruption and profile-response choices;
- its scientific examples and regression fixtures.

ITAMAE replaces only the internals listed in the phased roadmap below.

### 8.2 SASHIMI-W

SASHIMI-W requires an adapter-first migration because it uses a different unit
convention, WMAP7 input, global file loading, procedural initialization, and an
older SciPy/NumPy style.

The first objective is to expose:

```text
WDMVarianceModel
WDMConcentrationModel
WDMInitialStructureModel
```

behind common protocols while reproducing the existing output. Only after
regression tests pass should shared post-accretion and catalog machinery be
replaced.

### 8.3 SASHIMI-SI

SASHIMI-SI should use common CDM accretion and tidal-history infrastructure but
retain its own SIDM module. Its paired CDM/SIDM outputs indicate that the
catalog must support multiple named state views for one population node,
rather than duplicating the entire catalog.

A possible representation is:

```text
catalog.states["cdm_reference"]
catalog.states["sidm"]
```

with shared `m_acc`, `z_acc`, quadrature nodes, and base weights. SIDM validity,
survival, and collapse flags remain state-specific.

### 8.4 SASHIMI-F

SASHIMI-F provides the strongest motivation for composable power-spectrum and
variance objects. The FDM particle-mass setter currently triggers
power-spectrum and interpolation reconstruction. In ITAMAE this should become
an immutable, cacheable configuration object.

The FDM population suppression and FDM internal core-halo structure must be
separate components. The former enters the variance/accretion model; the latter
enters initial/profile structure and downstream priors.

## 9. Numerical and software policies

### 9.1 Units

ITAMAE should choose one documented internal unit convention for the first
release and provide explicit conversion at legacy adapters. Silent mixing of
the normalized units used by current C/SI/F and cgs units used by W is not
acceptable.

A lightweight internal unit convention is preferred over carrying Astropy
quantities through every large batch. Astropy may be used at boundaries and in
validation tests.

### 9.2 Array API and vectorization

All public numerical methods should accept scalar and NumPy-array inputs with
documented broadcasting rules. Batch shape must be separate from integration
axes. Hard-coded reshaping based on `N_ma`, `N_herm`, or `len(zdist)` should be
confined to legacy adapters.

Initial backend: NumPy/SciPy. JAX or other differentiable backends are future
work and should not shape the first API unless the abstraction is cost-free.

### 9.3 Caching

Cache keys must include:

- physical model parameters;
- cosmology;
- power-spectrum source and version;
- integration bounds and resolution;
- code version or schema version.

The current practice of excluding `self` from a pickle-cache key is unsafe when
instances differ physically.

### 9.4 Parallelism

Parallel execution belongs in a runner, not in physical model functions.
Models should expose pure batch-in/batch-out methods. The runner may chunk over
accretion redshift, mass, concentration, or orbital nodes. Serial execution
must remain the reference backend.

### 9.5 Flags and failure handling

Do not globally suppress `RuntimeWarning`. Numerical issues should be encoded
using explicit flags such as:

```text
INVALID_VARIANCE
OUTSIDE_INTERPOLATION_RANGE
ROOT_NOT_BRACKETED
PROFILE_UNPHYSICAL
DISRUPTED
COLLAPSED
ORBIT_UNBOUND
```

Model packages decide whether flagged nodes are removed, assigned zero weight,
or reported to the user.

## 10. Testing strategy

### 10.1 Golden regression data

Before extraction, each SASHIMI repository should generate compact reference
outputs for a small matrix of configurations:

- multiple host masses;
- redshift 0 and at least one nonzero target redshift;
- low and standard grid resolutions;
- concentration scatter on/off;
- evolved/unevolved profile options;
- representative WDM, SIDM, and FDM parameters;
- at least one spatial/radial configuration.

Store summary arrays and selected raw nodes, not only plotted results.

### 10.2 Invariant tests

Required invariants include:

- nonnegative finite quadrature weights;
- agreement between integrated weights and total expected abundance;
- normalization of concentration and spatial quadratures;
- `m_bound <= m_acc` for models that assume monotonic stripping;
- profile mass consistency at the truncation radius;
- scalar/batch equivalence;
- convergence with grid and interpolation resolution;
- serial/parallel equivalence;
- recovery of the original non-spatial result after spatial marginalization
  when the radial correction is normalized to unity.

The last condition is particularly important for the `r-dependent` extension.

## 11. Phased roadmap

### Phase 0: baseline and repository skeleton

- add `pyproject.toml`, `src/itamae`, and a minimal test configuration;
- define supported Python/NumPy/SciPy versions;
- add golden-output scripts to C, W, SI, and F;
- document unit and naming conventions;
- make no scientific changes.

### Phase 1: types, numerics, and NFW primitives

Implement:

- state and catalog dataclasses;
- grid and Gauss-Hermite helpers;
- integration/interpolation wrappers;
- safe NFW mass and inverse-mass functions;
- mass-definition conversion utilities;
- flags, serialization, and regression helpers.

Adopt these first in SASHIMI-C while preserving its public tuple-returning API
through a compatibility wrapper.

### Phase 2: generic evolution solver

Extract and test:

- ODE runner;
- perturbative coefficients and solution assembly;
- optional Shanks acceleration;
- history-grid handling;
- batch/chunk execution.

The SASHIMI-C mass-loss law remains in SASHIMI-C and is supplied to the solver.

### Phase 3: initial measure and common catalog builder

Implement:

- accretion-batch construction;
- independent weight factors;
- concentration quadrature;
- model-component runner;
- deterministic weighted catalog and Monte Carlo realization.

Migrate C first, then SI. Do not migrate W/F until their variance interfaces
have regression coverage.

### Phase 4: power spectrum and variance protocols

Implement composable:

- tabulated power spectrum;
- transfer function;
- top-hat and sharp-k windows;
- variance integration and derivative;
- interpolation/cache objects.

Integrate W and F through adapters. Compare against their original grids and
mass functions before optimizing.

### Phase 5: spatial Level A and Level B

Implement:

- radial-node measure and explicit spatial weights;
- host-potential/local-environment interface;
- normalized conditional radial PDFs;
- context-aware mass-loss and survival protocols;
- shell counts and radial number-density observables.

Reimplement the useful parts of `r-dependent` using these interfaces. The
existing fitting functions remain in SASHIMI-C until their physical status is
settled.

### Phase 6: SIDM and FDM structural extensions

- support multiple named state views;
- add profile-parameter schemas for cored/SIDM and soliton/FDM profiles;
- provide adapters to dSph profile/prior packages;
- keep Jeans and inference code outside ITAMAE.

### Phase 7: orbit-averaged Level C research backend

After the radial implementation is validated:

- implement spherical turning points and radial periods;
- cache orbit kernels in dimensionless host coordinates where possible;
- implement infall-to-`(E,L)` transformations;
- prototype conservative transport with source/drift/sink operators;
- compare radial distributions and mass segregation with simulations and
  simpler radial-PDF models.

This phase is research work and should not block a stable non-spatial ITAMAE
release.

## 12. Initial public API target

ITAMAE should initially expose low-level, stable objects rather than a universal
model constructor:

```python
from itamae.types import AccretionBatch, WeightedSubhaloCatalog
from itamae.numerics import gauss_hermite_lognormal
from itamae.halo import NFWProfile, convert_mass_definition
from itamae.evolution import PerturbativeEvolutionSolver
```

A SASHIMI package assembles them:

```python
from sashimi_c import SashimiCDM

model = SashimiCDM()
catalog = model.generate_catalog(host_mass=1.0e12, redshift=0.0)
```

Later, advanced users may construct a custom pipeline explicitly, but ITAMAE
must not advertise a default physical configuration as the canonical SASHIMI
calculation.

## 13. Immediate next tasks

The first implementation pull request should be limited to:

1. repository packaging and CI;
2. `HostState`, `AccretionBatch`, `SubhaloState`, `OrbitalState`, and
   `WeightedSubhaloCatalog` drafts;
3. grid and Gauss-Hermite utilities;
4. NFW enclosed mass and robust inversion;
5. unit tests and one small SASHIMI-C compatibility example.

It should not yet move EPS, concentration, stripping coefficients, or
SIDM/WDM/FDM physics into ITAMAE.

After that pull request, the second milestone is extraction of the generic
perturbative solver with regression tests against SASHIMI-C and SASHIMI-SI.

## 14. Definition of success

The initial ITAMAE refactor is successful when:

- each SASHIMI variant still clearly owns its physical assumptions;
- common numerical code is no longer copied among repositories;
- all variants produce regression-equivalent results;
- weighted catalogs use a shared schema;
- WDM/FDM variance implementations can be exchanged without changing EPS
  integration code;
- SIDM/FDM-specific profile fields do not require new long tuple signatures;
- a radial measure can be added without rewriting the non-spatial pipeline;
- the original global result can be recovered by spatial marginalization;
- future orbit-averaged work can reuse the same host, state, measure, and
  catalog abstractions.
