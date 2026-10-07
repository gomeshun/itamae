# Population execution contract

`PopulationPipeline` owns stage order, validation, weight transport and catalog
assembly. `PopulationComponents` supplies named variant-owned implementations of
those same stages. No physical prescription is selected by the executor.

Every batch uses one flat population axis. All metadata keys and nested values
must match exactly across batches, including source revisions, units, model
specification and backend choices. Mapping insertion order is immaterial; numeric
metadata is compared without a scientific tolerance. Batch-local indices belong
in `contexts`; node-local coordinates belong in aligned arrays. Optional weight
factors must be present in every batch or absent from every batch. These checks
run before any callback so incompatible inputs never partially execute.

Stages return mappings of named finite real or boolean arrays with exactly the
batch shape. The columns stage must provide at least one physical column; its
names, and the survival view names, must match across batches. ITAMAE does not
impose NFW columns or a variant's mass definition. Survival masks contain booleans
or binary 0/1 values; arbitrary truthy numbers are invalid. Callback and stage
validation failures report the batch index, stage and callback, preserving the
original exception as the cause.

An empty input iterable is an error because no output schema exists. A zero-node
batch is valid and invokes the stages on empty arrays. A zero-survivor population
retains every physical node with zero survival weights. No clipping, normalization
or removal of invalid physical values occurs at this boundary.

Execution follows input order. To partition a population, split all node/weight
arrays identically and repeat the appropriate physical context. Pointwise stages
must yield the same ordered columns, masks and weights as unpartitioned execution.
Components using batch-wide normalization or randomness must supply that state
explicitly; changing the partition must not change the physical model. The
independent toy integration test checks exact transport equality and seeded
catalog realizations for a partition that includes a zero-node chunk.

`AccretionBatch` and count catalogs require nonnegative independent weights.
Historical signed measures belong in independent reference calculations, not
production count catalogs. The final count weights multiply the independent
batch factors and the selected survival mask exactly once.

## Validation record (2026-09-10)

The new contract tests failed in 18 cases before the fix (3 passed). With the
fix, the complete ITAMAE suite passes 80 tests on Python 3.11.15 using the locked
full-extra dependencies; branch coverage is 79% versus 78% before the change.
Ruff lint/format and mypy pass. C's 53 and SI's 14 existing regressions pass with
their batch-local redshift indices moved into contexts. This is execution
validation, not a new claim about the physical prescriptions.
