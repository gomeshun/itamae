# ITAMAE Implementation Plan

**Integrated Toolkit for Analytical Merger-tree And Evolution**

This document defines the initial architecture and migration plan for ITAMAE after comparing the current implementations of:

- [SASHIMI-C](https://github.com/gomeshun/sashimi-c)
- [SASHIMI-W](https://github.com/gomeshun/sashimi-w)
- [SASHIMI-SI](https://github.com/gomeshun/sashimi-si)
- [SASHIMI-F](https://github.com/gomeshun/sashimi-f)
- the FDM halo/profile implementation in [dsph_fuzzy](https://github.com/gomeshun/dsph_fuzzy)

The plan is deliberately conservative: ITAMAE should remove duplicated machinery without taking ownership of the scientific identity of each SASHIMI variant.

## 1. Architectural principle

ITAMAE owns reusable **mechanisms**. Each SASHIMI package owns its concrete **physical prescriptions and default composition**.

> ITAMAE knows how to perform a SASHIMI calculation.  
> A SASHIMI variant specifies which calculation represents its dark-matter model.

SASHIMI-C must therefore remain a scientifically meaningful package rather than becoming an almost empty wrapper.

### ITAMAE owns

- common data structures and array conventions;
- numerical integration, interpolation, root finding, and ODE infrastructure;
- quadrature, sampling, batching, parallel execution, and caching;
- configurable cosmology and spherical-halo primitives;
- interfaces for power spectra, variance, collapse barriers, concentration, accretion, mass loss, structural response, disruption, and spatial evolution;
- generic weighted-catalog operations;
- state/history containers, serialization, validation, and regression utilities.

### Each SASHIMI variant owns

- its reference cosmology where scientifically required;
- its concrete power-spectrum and transfer-function choices;
- collapse barrier and EPS prescriptions;
- concentration calibration;
- calibrated tidal-mass-loss and structural-response laws;
- model-specific profile evolution and survival criteria;
- published observables, reference configurations, and compatibility APIs.

Concrete named models may eventually live in optional `itamae.models.*` modules, but the ITAMAE core must never select them silently.

## 2. Findings from the current family

### 2.1 Strongly duplicated components

SASHIMI-C, SASHIMI-SI, and SASHIMI-F contain closely related implementations of:

- units and constants;
- background cosmology, linear growth, critical density, and virial overdensity;
- NFW functions and conversions among `M200`, `Mvir`, `rs`, `rhos`, `Vmax`, and `rmax`;
- host mass-accretion histories (`Mzi`, `Mzzi`, and `dMdz`);
- Yang et al.-type EPS accretion calculations;
- Gauss-Hermite treatment of concentration and host-history scatter;
- post-accretion tidal mass-loss equations;
- zeroth- through third-order perturbative solutions and Shanks acceleration;
- construction and flattening of weighted subhalo catalogs.

These are the first extraction targets.

### 2.2 SASHIMI-W requires genuinely model-independent interfaces

SASHIMI-W uses its own unit convention and reference cosmology, reads a tabulated matter power spectrum, applies a WDM transfer function, computes sharp-k and top-hat variances, uses a mass-dependent WDM collapse barrier, and calculates a WDM concentration relation. These quantities then enter an EPS and tidal-evolution pipeline similar in broad structure to SASHIMI-C.

Therefore the following must be independent interfaces:

- `PowerSpectrum`;
- `TransferFunction`;
- `WindowFunction`;
- `VarianceModel`;
- `CollapseBarrier`;
- `ConcentrationModel`.

In particular, ITAMAE must not define `sigma(M,z)` as a CDM-only fitting function.

### 2.3 SASHIMI-SI requires complete evolutionary histories

SASHIMI-SI first builds a CDM-like mass, `Vmax`, and `rmax` history, then integrates the SIDM parametric/gravothermal response along that history. It also distinguishes CDM and SIDM survival and weights and stores quantities such as `tt_ratio` and core-collapse state.

ITAMAE must therefore support:

- trajectories on common or ragged redshift/time grids;
- model components that consume another component's trajectory;
- multiple named states for one initial node, such as `cdm_reference` and `sidm`;
- explicit model-dependent masks, probabilities, and weights.

Evolution cannot be represented only as a map from an accretion state to one final state.

### 2.4 SASHIMI-F separates population suppression from halo structure

SASHIMI-F modifies the power spectrum and variance entering EPS and supports top-hat, sharp-k, and external/Colossus-style backends. The FDM models in `dsph_fuzzy` separately implement a soliton core joined to an outer NFW-like profile, smooth/continuous transition conditions, enclosed masses, truncation, `r200`, and core mass.

ITAMAE should distinguish:

1. **population physics**: the abundance and initial distribution of accreted halos;
2. **internal-structure physics**: profile assignment and structural evolution;
3. **inference adapters**: conversion of a catalog/profile into satellite priors or Jeans-model parameters.

Jeans likelihoods, MCMC, WBIC, and project-specific dSph analyses remain outside ITAMAE.

### 2.5 The current monolithic classes mix independent concerns

The present scripts commonly mix:

1. physical formulas;
2. numerical solution methods;
3. multidimensional grid construction;
4. output formatting and observables.

ITAMAE should separate these concerns rather than merely moving duplicated functions into another monolithic file.

## 3. Core mathematical model

ITAMAE will treat a deterministic SASHIMI catalog as a quadrature representation of an initial measure followed by one or more evolution maps.

Let

\[
  x=(m_{\rm acc},z_{\rm acc},c_{\rm acc},\eta,\ldots),
\]

with initial measure

\[
  dN=\lambda(x)\,dx.
\]

A model-specific evolution map produces a state or trajectory

\[
  y=T_\theta(x),
\]

and an observable is

\[
  \langle O\rangle=\int \lambda(x)O[T_\theta(x)]\,dx.
\]

The current weighted catalogs are quadrature rules for this integral. ITAMAE will preserve the deterministic weighted representation as the canonical internal form. Monte Carlo and quasi-Monte Carlo catalogs will be derived representations.

For future spatial calculations, the state becomes

\[
  y=(m,\mathrm{profile},\boldsymbol{x},\boldsymbol{v},\mathrm{orbit},\mathrm{flags},\ldots),
\]

without changing the measure/catalog abstraction.

## 4. Proposed package structure

```text
src/itamae/
  __init__.py
  constants.py
  typing.py
  exceptions.py

  cosmology/
    base.py
    flat_lcdm.py
    time.py

  numerics/
    grids.py
    quadrature.py
    integration.py
    interpolation.py
    roots.py
    ode.py
    perturbation.py
    random.py

  power/
    base.py
    tabulated.py
    transfer.py
    windows.py
    variance.py

  halos/
    mass_definition.py
    concentration.py
    profiles/
      base.py
      nfw.py
      truncated_nfw.py
      soliton.py
      soliton_nfw.py

  histories/
    state.py
    trajectory.py
    host.py

  models/
    protocols.py
    composition.py

  measure/
    nodes.py
    weights.py
    eps_engine.py

  evolution/
    executor.py
    mass.py
    structure.py
    disruption.py
    spatial.py

  catalog/
    schema.py
    catalog.py
    views.py
    sampling.py
    io.py

  observables/
    generic.py

  adapters/
    legacy.py
    pandas.py
    xarray.py

  validation/
    checks.py
    regression.py
```

Directory names may evolve, but these conceptual boundaries should remain.

## 5. Stable core data model

### 5.1 Array convention

The unflattened canonical quadrature shape should be

```text
(z_acc, concentration_node, mass_node, optional_extra_node)
```

Core calculations must preserve named dimensions or explicit dimension metadata. Flattening should occur only in catalog views or legacy adapters.

The first implementation can use NumPy arrays plus dimension metadata. `xarray` should be an optional adapter rather than a mandatory dependency.

### 5.2 `InitialNodeBatch`

Required fields:

```text
m200_acc
mvir_acc
z_acc
concentration_node
weight_base
weight_quadrature
```

Optional fields:

```text
formation_redshift
host_history_node
model_parameters
spatial_initial_conditions
```

### 5.3 `HaloStateBatch`

Common fields:

```text
redshift
cosmic_time
bound_mass
m200
mvir
rs
rhos
rmax
vmax
truncation_radius
truncation_concentration
survive
status_flags
```

Model-specific fields should be namespaced or registered, for example:

```text
sidm.rc
sidm.tt_ratio
sidm.core_collapsed
fdm.rc
fdm.rt
fdm.rho_c
```

### 5.4 `HaloTrajectoryBatch`

This stores states along redshift or cosmic time and must support:

- dense common grids for vectorized calculations;
- ragged trajectories through offsets/index arrays when required;
- interpolation onto another component's grid;
- derived rates and event flags;
- optional orbital positions and velocities.

This object is required by SIDM and is also the natural foundation for future spatial/orbital evolution.

### 5.5 `WeightedSubhaloCatalog`

Weights must not be overwritten destructively. Store at least:

```text
weight_base
weight_concentration
weight_host_history
weight_orbit
weight_survival
weight_selection
weight_final
```

Masks and probabilities remain separate:

```text
is_accreted
is_disrupted
is_numerically_valid
survival_probability
selection_probability
```

This is necessary to compare CDM, WDM, SIDM, and FDM on a common initial measure and to diagnose the origin of suppression.

## 6. Model protocols

Composition should replace the current deep inheritance pattern.

```python
class VarianceModel(Protocol):
    def sigma(self, mass, redshift): ...
    def variance(self, mass, redshift): ...
    def dvariance_dmass(self, mass, redshift): ...

class CollapseBarrier(Protocol):
    def delta_c(self, mass, redshift): ...

class HostHistoryModel(Protocol):
    def mass(self, host_mass_at_reference, redshift, reference_redshift): ...
    def dmass_dredshift(self, ...): ...

class AccretionMeasureModel(Protocol):
    def differential_number(self, nodes, host): ...

class ConcentrationModel(Protocol):
    def median(self, mass, redshift): ...
    def quadrature(self, mass, redshift, scatter, order): ...

class MassLossLaw(Protocol):
    def rhs(self, bound_mass, redshift, host_state): ...

class StructuralResponseModel(Protocol):
    def evolve(self, initial_state, mass_history, context): ...

class DisruptionModel(Protocol):
    def evaluate(self, state, context): ...

class SpatialModel(Protocol):
    def initialize(self, nodes, host): ...
    def evolve(self, initial_phase_space, histories, host): ...
```

Protocols express requirements without forcing every variant into one class hierarchy.

## 7. What belongs in ITAMAE immediately

### 7.1 Numerical utilities

- logarithmic and mixed grids;
- Gauss-Hermite nodes and weights;
- weighted Simpson/trapezoid integration;
- stable interpolation and monotonic inverse interpolation;
- reusable root-solving utilities;
- ODE wrappers;
- generic polynomial perturbation-series machinery;
- Shanks transformation with diagnostics;
- chunked mapping and optional parallel backends;
- deterministic random-number handling.

The perturbative solver is shared, while coefficients generated from a calibrated physical mass-loss law remain model-owned.

### 7.2 Cosmology and halo primitives

- configurable flat-LCDM background;
- `H(z)`, `rho_crit(z)`, growth factor, growth derivative, lookback time, and cosmic age;
- spherical-overdensity definitions;
- `M200`/`Mvir` conversion infrastructure;
- NFW density, enclosed mass, circular velocity, `Vmax`, `rmax`, and truncation inversion;
- unit conversion helpers at package boundaries.

ITAMAE must not hard-code Planck18 or WMAP7 as the only cosmology. A SASHIMI variant selects its reference cosmology.

### 7.3 Power-spectrum and variance infrastructure

- tabulated power-spectrum loader with units and provenance;
- transfer-function composition;
- top-hat and sharp-k windows;
- stable small-argument spherical-Bessel evaluation;
- numerical `sigma(M,z)` and `dS/dM` calculation;
- interpolation/cache layer;
- k-range convergence diagnostics;
- optional external-backend adapters.

Concrete WDM and FDM transfer functions may initially remain in their SASHIMI repositories while using ITAMAE interfaces. They can move later to optional named model modules after regression equivalence and ownership are agreed.

### 7.4 Generic EPS engine

ITAMAE should contain tensor/grid machinery and integration infrastructure for conditional accretion calculations:

- host-history scatter quadrature;
- broadcasting over mass and accretion-redshift nodes;
- normalization integrals;
- base catalog-weight construction;
- validation of positive variance differences and finite kernels.

The Yang model-1/2/3 kernels and the WDM barrier are scientific prescriptions. During the first migration they should remain in the corresponding SASHIMI package and be passed into the engine.

### 7.5 Evolution execution

ITAMAE should provide:

- scalar and vectorized mass-evolution execution;
- final-state and full-history modes;
- common redshift/time-grid construction;
- trajectory interpolation;
- component pipelines in which one evolution model consumes another model's reference history;
- event/status propagation;
- chunking and parallel execution.

Calibrated `A(M,z)`, `zeta(M,z)`, Peñarrubia-like response, SIDM gravothermal response, and FDM core relations remain model-owned.

### 7.6 Catalog and generic observables

ITAMAE should provide only generic operations:

- weighted histograms and cumulative distributions;
- weighted sums, expectations, and moments;
- filtering and named views;
- Poisson realization from a weighted expectation catalog;
- reproducible resampling;
- optional NPZ/HDF5/Parquet-style serialization;
- conversion hooks for legacy tuples.

Published SASHIMI observables, such as specific annihilation-boost or satellite-formation prescriptions, remain in the scientific package.

### 7.7 Profile primitives needed by FDM downstream work

ITAMAE can host reusable mathematical profile primitives:

- soliton density and enclosed mass;
- NFW outer envelope;
- continuous and smooth soliton--NFW matching;
- truncation and total mass;
- derived `r200`, core mass, and circular velocity.

The selected core--halo relation and its calibration remain in SASHIMI-F or the downstream FDM package.

## 8. What remains in each SASHIMI package

### SASHIMI-C

- canonical CDM composition and defaults;
- selected CDM variance/power-spectrum configuration;
- the published host-history and concentration prescriptions;
- calibrated mass-loss law;
- CDM structural response and disruption definition;
- prompt-cusp extension and associated observables;
- published-reference scripts and public API;
- legacy `subhalo_properties_calc` compatibility.

### SASHIMI-W

- WDM particle-mass parameterization;
- WDM transfer function and reference cosmology/data selection;
- WDM collapse barrier;
- WDM concentration calibration;
- WDM-specific mass limits and validity checks;
- published WDM observables and reference configurations.

### SASHIMI-SI

- scattering cross-section models;
- effective-cross-section calculation;
- collapse-time prescription;
- SIDM parametric/gravothermal evolution;
- SIDM survival and core-collapse definitions;
- dual CDM-reference/SIDM output semantics;
- published SIDM observables.

### SASHIMI-F

- FDM particle-mass convention;
- FDM transfer function and selected window/backend defaults;
- FDM population-suppression model;
- FDM concentration and core--halo prescriptions;
- FDM-specific validity domains;
- published FDM observables and dSph adapters.

## 9. Spatial and orbital extension

Position information must be included in the architecture now, even if the first release does not implement a calibrated spatial model.

### 9.1 Do not append only a final radius column

A final `r` value alone is insufficient for applications involving:

- radial selection effects;
- disk shocking and baryonic disruption;
- orbit-dependent stripping;
- lensing projections;
- gamma-ray and stellar-stream observables;
- correlations between accretion time, stripping, and present position.

Spatial evolution should be an optional phase-space trajectory coupled to mass and structural histories.

### 9.2 Spatial state

Reserve fields such as:

```text
position_cartesian       # (..., 3)
velocity_cartesian       # (..., 3)
host_centric_radius
radial_velocity
tangential_velocity
specific_energy
specific_angular_momentum
pericenter
apocenter
infall_direction
orbit_status
```

Not every spatial model must populate every field.

### 9.3 Separate spatial ingredients

Define independent protocols for:

1. `InfallPhaseSpaceModel` — phase-space distribution at accretion;
2. `HostPotentialModel` — spherical, axisymmetric, time-dependent, or externally supplied potential;
3. `OrbitEvolutionModel` — orbit integration or empirical radial mapping;
4. `OrbitDependentMassLossModel` — optional feedback from orbit to stripping;
5. `SpatialSelectionModel` — observer, survey, lensing, or aperture selection;
6. `ProjectionModel` — conversion from 3D coordinates to projected observables.

A specific radial distribution must not be baked into the catalog generator.

### 9.4 Explicit coupling levels

Support three levels:

- **Level 0: post-processing position model**  
  Draw or weight present-day positions conditionally on existing halo properties, with no feedback to mass loss.

- **Level 1: orbit-aware tagging**  
  Generate an orbital history and derived pericenter/apocenter while retaining the existing calibrated orbit-averaged mass-loss law.

- **Level 2: coupled orbital evolution**  
  Integrate orbit, mass loss, structural response, and baryonic effects together.

Current SASHIMI calculations are approximately orbit-averaged. ITAMAE must preserve that fast path and impose no spatial overhead when spatial evolution is disabled.

### 9.5 Spatial measure and weights

Spatial variables may be represented by additional quadrature nodes or conditional sampling:

\[
 dN=\lambda(m_{\rm acc},z_{\rm acc},c_{\rm acc})
 p(\Gamma_{\rm acc}\mid m_{\rm acc},z_{\rm acc},c_{\rm acc})
 \,d\Gamma_{\rm acc}\,dm\,dz\,dc,
\]

where `Gamma_acc` is the phase-space state at accretion.

The catalog should distinguish:

```text
weight_population
weight_orbit
weight_survival
weight_selection
```

This permits deterministic radial distributions, Monte Carlo orbit realizations, and hybrid methods without redesigning the pipeline.

## 10. Public API direction

The model-facing API should be small and composition-based.

```python
from itamae import CatalogRequest
from itamae.measure import InitialMeasureEngine
from itamae.evolution import EvolutionExecutor

request = CatalogRequest(
    host_mass=1.0e12,
    evaluation_redshift=0.0,
    accretion_redshift_grid=...,
    accretion_mass_grid=...,
    concentration_quadrature=...,
    retain_history=False,
    retain_phase_space=False,
)

nodes = measure_engine.build_nodes(request, model=cdm_model)
catalog = evolution_executor.run(nodes, model=cdm_model, request=request)
```

End users should normally call the SASHIMI package:

```python
from sashimi_c import SashimiCDM

catalog = SashimiCDM().generate_catalog(...)
```

Direct low-level composition is intended for model development and cross-model studies.

## 11. Legacy compatibility

Each SASHIMI package should provide an adapter reproducing its current tuple output exactly.

```python
legacy_tuple = catalog.to_legacy("sashimi-c-v1")
```

Regression tests should compare:

- node-by-node catalog quantities before survival filtering;
- final filtered catalogs;
- integrated total weights;
- mass and velocity functions;
- representative profile quantities;
- published benchmark observables.

Generic adapter infrastructure belongs in ITAMAE; model-specific field order belongs in the relevant SASHIMI package.

## 12. Testing strategy

### Unit tests

- cosmology identities and dimensions;
- NFW analytic relations;
- mass-definition round trips;
- quadrature normalization;
- variance derivatives against finite differences;
- perturbative solutions against direct ODE integration;
- catalog weight bookkeeping;
- serialization round trips.

### Cross-repository regression tests

Freeze small benchmark outputs for C, W, SI, and F. Use small grids for CI and larger optional scientific benchmarks.

Suggested benchmark matrix:

```text
host mass:       1e9, 1e12, 1e15 Msun where valid
redshift:        0, 1
mass nodes:      16 or 32
z nodes:         coarse fixed grid
concentration:   1 and 5 Hermite nodes
model points:    one CDM, two WDM, two SIDM, two FDM settings
```

### Invariants

- base weights are finite and non-negative;
- total weight is invariant under flatten/unflatten operations;
- survival masks never alter base weights;
- final mass does not exceed accretion mass unless explicitly permitted;
- profile mass and stored bound mass agree within tolerance;
- trajectory redshift/time axes are monotonic;
- spatial states satisfy their declared coordinate convention.

## 13. Performance strategy

Optimization follows numerical equivalence.

1. Vectorize over mass and concentration nodes.
2. Chunk over accretion redshift or node index when memory-limited.
3. Cache expensive power-spectrum and variance tables with complete configuration hashes.
4. Support optional multiprocessing/joblib execution.
5. Avoid mandatory JAX, Numba, or NumExpr dependencies in the core.
6. Add accelerated backends only behind identical protocols.
7. Benchmark final-state and full-history modes separately.

Cache keys must include cosmology, model parameters, grid ranges, backend version, and tolerances. Ignoring instance state in cache keys is unsafe when that state changes the result.

## 14. Migration phases

### Phase 0 — Baseline and scaffolding

- add `pyproject.toml`, `src/itamae`, tests, linting, and CI;
- define supported Python/NumPy/SciPy versions;
- freeze reference outputs from all four SASHIMI variants;
- document units and array conventions;
- add architecture decision records for ownership boundaries.

### Phase 1 — Pure numerical and halo primitives

Extract without changing scientific formulas:

- grids and quadrature;
- configurable cosmology base;
- NFW/truncated-NFW utilities;
- mass-definition conversions;
- integration/interpolation/root helpers;
- perturbation and Shanks utilities.

Migrate SASHIMI-C to these pieces while preserving its API.

### Phase 2 — Catalog schema and weight bookkeeping

- implement `InitialNodeBatch`, `HaloStateBatch`, and `WeightedSubhaloCatalog`;
- add named weights and masks;
- implement flattening and the SASHIMI-C legacy adapter;
- migrate generic weighted histograms and Poisson realization;
- verify agreement with SASHIMI-C.

### Phase 3 — Generic EPS and host-history engine

- define variance, barrier, host-history, concentration, and accretion protocols;
- implement tensor/quadrature execution;
- keep concrete Yang/CDM and WDM prescriptions in their packages initially;
- migrate C, then W and F population generation;
- compare pre-evolution measures before moving further.

### Phase 4 — Evolution histories

- implement final-state and full-history execution;
- separate mass-loss law from numerical solver;
- migrate SASHIMI-C tidal evolution;
- add structural-response interfaces;
- migrate SASHIMI-SI using a CDM reference trajectory;
- add explicit event and numerical-failure flags.

### Phase 5 — WDM/FDM power infrastructure

- implement tabulated spectra, transfer composition, windows, variance, and derivatives;
- migrate SASHIMI-W and SASHIMI-F backends;
- add k-range convergence warnings and cache provenance;
- compare top-hat and sharp-k results against current codes.

### Phase 6 — FDM profiles and inference adapters

- implement soliton and soliton--NFW primitives;
- port tested enclosed-mass and matching utilities from `dsph_fuzzy`;
- define a core--halo-relation protocol;
- create adapters without moving Jeans likelihoods into ITAMAE.

### Phase 7 — Spatial foundation

- add optional spatial fields and phase-space support to `HaloTrajectoryBatch`;
- define infall, potential, orbit, projection, and selection protocols;
- migrate the experimental SASHIMI-C position-dependent implementation first as an external model using these interfaces;
- validate radial distributions and normalization;
- preserve a zero-overhead orbit-averaged default path.

### Phase 8 — Coupled spatial evolution and baryons

- add optional orbit-dependent mass loss;
- add time-dependent host potential and disk/baryonic components;
- support pericenter events and disruption channels;
- quantify differences from the orbit-averaged calibration before adopting any new default.

## 15. First usable milestone

The first usable ITAMAE release should satisfy all of the following:

- SASHIMI-C uses ITAMAE numerical, halo, and catalog infrastructure;
- the old SASHIMI-C public API remains available;
- reference CDM catalog quantities agree within documented tolerances;
- base weights, masks, probabilities, and final weights remain separate;
- full-history storage exists even if disabled by default;
- protocols required by W, SI, F, and spatial extensions are present;
- no WDM, SIDM, FDM, or spatial physical prescription is silently hard-coded into the core.

This milestone is deliberately smaller than a complete migration of every repository.

## 16. Immediate implementation tasks

1. Add packaging and CI scaffolding.
2. Write `UNITS.md` defining canonical internal units and boundary conversions.
3. Write `SCHEMA.md` defining dimensions, required fields, weights, masks, and status flags.
4. Implement NumPy-based quadrature and grid utilities.
5. Implement configurable flat-LCDM and NFW primitives.
6. Implement mass-definition conversions with regression tests.
7. Implement `WeightedSubhaloCatalog` and legacy-view hooks.
8. Implement generic perturbation/Shanks helpers separated from the CDM mass-loss law.
9. Add frozen small-grid outputs from SASHIMI-C.
10. Refactor SASHIMI-C on a migration branch to consume these components.
11. Add variance and barrier protocols before attempting WDM/FDM migration.
12. Add trajectory and optional spatial schemas before freezing the public API.

## 17. Open decisions

The following should be resolved through small prototypes:

- canonical internal units: current SASHIMI conventions versus explicit astrophysical units;
- dataclasses plus NumPy versus an optional `xarray`-first representation;
- dense versus ragged storage for histories;
- whether named scientific implementations live in `itamae.models` or only in SASHIMI repositories;
- primary on-disk format;
- minimum spatial representation for the first spatial milestone;
- regression tolerances for legacy interpolation and ODE behavior.

The default choice should minimize dependencies and maximize reproducibility.

## 18. Definition of success

ITAMAE succeeds when:

- a physical formula is implemented once unless the models genuinely differ;
- differences among C, W, SI, and F are visible as explicit components and configurations;
- SASHIMI-C remains a scientifically meaningful package rather than an empty wrapper;
- weighted expectation catalogs and stochastic realizations share one schema;
- SIDM and future spatial models can consume full histories without redesigning the package;
- new dark-matter or baryonic models are added by implementing protocols rather than copying an entire SASHIMI script;
- published SASHIMI results remain reproducible and attributable to the correct scientific package.
