# ITAMAE

**Integrated Toolkit for Analytical Merger-tree And Evolution**

ITAMAE is a shared computational toolkit for the
[SASHIMI](https://github.com/gomeshun/sashimi-c) family of semi-analytical
subhalo models.

The name follows the culinary theme of SASHIMI: an *itamae* is the chef who
prepares and assembles the ingredients. In the same way, ITAMAE provides the
common numerical machinery and data structures used to construct SASHIMI
subhalo catalogs.

> **Project status:** early alpha migration.
> The package foundation, generic evolution, weighted-catalog, power/variance,
> and spherical-orbit mechanisms are implemented. The public API remains
> provisional until all SASHIMI variants complete their golden regressions.

## Motivation

Several SASHIMI variants share substantial computational infrastructure,
including:

* cosmological background calculations;
* halo mass-definition and NFW-profile utilities;
* numerical integration, interpolation, and quadrature;
* solvers for post-accretion evolution;
* weighted subhalo-catalog data structures;
* sampling, caching, parallel execution, and input/output utilities.

Historically, parts of this machinery have been implemented independently in
SASHIMI-C, SASHIMI-SIDM, SASHIMI-WDM, and FDM-related projects. ITAMAE aims to
remove unnecessary duplication while keeping the scientific assumptions of
each SASHIMI variant explicit and independently testable.

## Design principle

ITAMAE provides the **computational mechanisms** required by SASHIMI models.
Each SASHIMI package remains responsible for selecting and defining its
**physical model**.

In particular:

* **ITAMAE** owns reusable interfaces, numerical solvers, catalog schemas,
  profile primitives, and execution infrastructure.
* **SASHIMI-C** owns the standard CDM model choices, including its accretion,
  concentration, tidal-evolution, disruption, and observable prescriptions.
* **SASHIMI-SIDM, SASHIMI-WDM, and SASHIMI-FDM** own their corresponding
  dark-matter-specific extensions and calibrations.

A concise summary is:

> **ITAMAE knows how to perform a SASHIMI calculation.
> Each SASHIMI variant knows which physical calculation should be performed.**

ITAMAE is therefore not intended to replace SASHIMI-C or any other SASHIMI
variant with a single monolithic implementation.

## Planned scope

The initial development will focus on the smallest set of components that are
genuinely shared across the SASHIMI family.

### Core numerical infrastructure

* logarithmic mass and redshift grids;
* Gauss-Hermite and other quadrature utilities;
* numerical integration and interpolation;
* generic ordinary differential-equation solvers;
* perturbative evolution solvers;
* chunked and parallel batch execution;
* reproducible random-number handling;
* caching and serialization.

### Cosmology and halo utilities

* cosmological background quantities;
* critical density, growth factor, and cosmic time;
* halo mass-definition conversions;
* NFW and truncated-NFW profile primitives;
* conversions among
  (M_{200}), (M_{\mathrm{vir}}), (r_s), (\rho_s),
  (V_{\max}), and (r_{\max}).

### Common interfaces

ITAMAE will define interfaces for model components such as:

* variance and power-spectrum models;
* halo mass-accretion histories;
* subhalo accretion models;
* concentration models;
* post-accretion mass-loss models;
* density-profile evolution models;
* survival and disruption prescriptions.

Concrete scientific prescriptions may live in a SASHIMI variant or in an
optional ITAMAE model module, but they will not be silently selected by the
ITAMAE core.

### Weighted catalogs

A central object will be a structured weighted subhalo catalog containing, for
example:

* accretion mass and redshift;
* structural parameters at accretion;
* evolved bound mass and structural parameters;
* base quadrature weight;
* concentration-scatter weight;
* model-dependent survival weight;
* status and validity flags;
* model-specific optional fields.

The catalog layer will support both deterministic expectation-value
calculations and stochastic realizations.

## Non-goals

The initial ITAMAE package will not attempt to own:

* the canonical physical definition of SASHIMI-C;
* SIDM cross-section or gravothermal-evolution models;
* WDM- or FDM-specific calibrations;
* dwarf-galaxy Jeans likelihoods;
* MCMC, nested-sampling, or WBIC analysis pipelines;
* project-specific plotting scripts and notebooks.

These belong in the corresponding scientific packages. ITAMAE may provide
interfaces or adapters that allow them to exchange data cleanly.

## Relationship to the SASHIMI family

The intended dependency structure is:

```text
                         +------------------+
                         |      ITAMAE      |
                         | shared numerical |
                         | infrastructure   |
                         +---------+--------+
                                   |
             +---------------------+---------------------+
             |                     |                     |
      +------v------+       +------v------+       +------v------+
      |  SASHIMI-C  |       | SASHIMI-SI  |       | SASHIMI-W/F |
      | CDM physics |       | SIDM physics |       | WDM/FDM     |
      +-------------+       +-------------+       +-------------+
```

Downstream projects, including dwarf-spheroidal analyses, should normally
depend on the relevant SASHIMI variant rather than reconstructing the complete
physical pipeline directly from low-level ITAMAE components.

## Illustrative API

The API has not yet been finalized. The following example only illustrates the
intended separation between reusable infrastructure and model-specific
physics.

```python
from itamae.catalog import WeightedSubhaloCatalog
from itamae.integration import GaussHermiteQuadrature
from itamae.evolution import PerturbativeEvolutionSolver

# Model-specific prescriptions are supplied by a SASHIMI package.
from sashimi_c import SashimiCDM

model = SashimiCDM()
catalog: WeightedSubhaloCatalog = model.generate_catalog(
    host_mass=1.0e12,
    redshift=0.0,
)
```

## Development strategy

Development will proceed incrementally:

1. Record reference outputs from the existing SASHIMI implementations.
2. Extract only clearly identical numerical and structural components.
3. Add regression tests before changing any scientific implementation.
4. Migrate SASHIMI-C first while preserving its public results.
5. Integrate SIDM, WDM, and FDM extensions through explicit interfaces.
6. Optimize vectorization, batching, caching, and parallel execution only
   after numerical equivalence has been established.

Backward compatibility and scientific reproducibility are higher priorities
than rapid API expansion.

## Installation

ITAMAE is not yet released on PyPI. The public repository can be installed for
development with `uv`:

```bash
git clone https://github.com/gomeshun/itamae.git
cd itamae
uv sync --all-extras --group dev
uv run pytest
```

The canonical internal units and required reproducibility metadata are
documented in [`docs/canonical-units.md`](docs/canonical-units.md). Foundation
branch-consolidation decisions are recorded in
[`docs/foundation-integration.md`](docs/foundation-integration.md).

WDM and FDM packages retain their physical transfer functions and compose them
with ITAMAE explicitly:

```python
from itamae.power import SharpKWindow, TabulatedPowerSpectrum
from itamae.variance import IntegratedVarianceModel

power = TabulatedPowerSpectrum(
    k,
    p_modified,
    identifier="sashimi-variant:documented-transfer:v1",
)
variance = IntegratedVarianceModel(
    power=power,
    window=SharpKWindow(),
    rho_mean=rho_mean,
    k_min=k.min(),
    k_max=k.max(),
    filter_scale=calibrated_c,
)
```

The identifier, backend, mass grid, and numerical settings can be combined
with `variance_cache_key`; cache loading rejects a mismatched key rather than
silently reusing data from another physical model.

## Contributing

Issues and pull requests are welcome once the initial package structure has
been established.

Contributions should, where applicable, include:

* unit tests for numerical utilities;
* regression tests against an existing SASHIMI implementation;
* explicit documentation of physical assumptions;
* numerical-convergence or performance checks;
* clear separation between generic infrastructure and model-specific physics.

## Citation

A dedicated ITAMAE citation will be added when the package and its scientific
scope are stabilized.

When using a SASHIMI model, please cite the publications associated with the
specific SASHIMI variant and with the physical prescriptions used in the
calculation.

## License

The license has not yet been selected.
