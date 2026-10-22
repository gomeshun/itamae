# Scientific test review: itamae — 2026-10-14

Reviewed **145 test definitions** individually. This count is not the number of parametrized cases. The review proposes future changes; all test bodies, assertions and tolerances are retained.

Executable-code baseline before this documentation change: `40aa3bd225e0de614650e04c21d92f49b1999436`.
ITAMAE has no author-upstream baseline; all its new library definitions are in scope.

## Decisions

| Decision | Definitions |
| --- | ---: |
| Retain | 115 |
| Consolidate without losing conditions | 19 |
| Relax only incidental assertions | 9 |
| Retire or replace, conditionally | 2 |

A consolidation retains every independent scenario and positive/negative control. A relaxation removes an implementation constraint only after preserving the physical or provenance invariant. A retirement is conditional on the replacement/surviving coverage stated below; no candidate is an instruction to delete a test.

Exact equality remains appropriate for lossless serialization, unchanged input identity, deterministic weight factorization and algebraic contracts. A near-roundoff comparison with a frozen fixture protects reproduction of a specified calculation; it is not an independent physical accuracy requirement. Analytic solutions, independent quadrature, unit identities and controlled refinement have a separate scientific role.

## Validation of this documentation change

- NumPy docstrings parsed for 476 in-scope definitions; all explicit parameters are documented.
- All 74 tracked Python files have identical executable ASTs after removing only module/class/function docstrings. The test bodies, decorators and assertions are unchanged.
- Preserved reference and historical analysis files remain byte-for-byte unchanged.
- Python 3.11.15: 333 existing pytest cases passed, 0 failures/errors and 0 skips.
- Ruff checks pass; ITAMAE additionally passes its full format and mypy checks.
- The saved numerical sweeps/catalogues retain their original calculation revisions. This audit does not relabel them as runs of a later documentation commit.

## Coverage gap found while documenting units

The current `ColossusCosmology.rho_crit` and `rho_m` multiply Colossus density only by `1e9`. Colossus documents its values in `Msun*h^2/kpc^3`; conversion to the canonical `Msun/Mpc^3` also needs multiplication by the selected `h**2`. The existing cross-provider test checks only H and growth. This is a missing density-unit check, not a reason to remove that smoke test. The implementation and assertions are unchanged in this documentation review. A follow-up should check density against the selected H/G convention and explicitly cover h factors.

## Candidate index

| Test | Proposal | Reason and retained coverage |
| --- | --- | --- |
| [test_backend_config_validates_and_exports_metadata](../tests/test_backends.py#L8) | Relax only incidental assertions | Keep backend selection and metadata coverage, but the exact identifier prefix is a formatting choice that can be relaxed without weakening physical provenance. |
| [test_archive_rejects_conflicting_explicit_revision](../tests/test_build_source_identity.py#L212) | Consolidate without losing conditions | Consolidate with matching-override and malformed-override cases into one parametrized provenance test; preserve all conflict scenarios. |
| [test_exact_checkout_rejects_conflicting_explicit_revision](../tests/test_build_source_identity.py#L241) | Consolidate without losing conditions | Consolidate override cases while retaining both checkout and archive origins; they protect different identity sources. |
| [test_malformed_explicit_revision_does_not_silently_fall_back](../tests/test_build_source_identity.py#L269) | Consolidate without losing conditions | Consolidate malformed values into the common override matrix; rejecting ambiguous provenance remains required. |
| [test_matching_explicit_archive_revision_is_allowed](../tests/test_build_source_identity.py#L297) | Consolidate without losing conditions | Merge this positive control with conflicting-override coverage so the rejection tests cannot pass by rejecting every build. |
| [test_variance_identity_includes_derivative_resolution_but_not_batching](../tests/test_cache_contract.py#L57) | Relax only incidental assertions | Keep sensitivity to accuracy controls; consolidate the chunk-size identity assertion with chunk-independence tests rather than requiring an identifier design. |
| [test_execution_reuses_canonical_variance_protocol](../tests/test_execution.py#L41) | Retire or replace, conditionally | Retirement candidate: this is object identity only; test_execution_composition plus real variance-adapter tests exercise usable behavior. Retain only if this exact import identity is a supported API promise. |
| [test_population_components_delegate_to_canonical_pipeline](../tests/test_execution_composition.py#L242) | Retire or replace, conditionally | Retirement candidate: the adjacent callback-versus-component execution test verifies the behavior that matters. Requiring this concrete class can forbid equivalent implementations without improving scientific reliability. |
| [test_population_components_match_callback_execution_and_catalogs](../tests/test_execution_composition.py#L265) | Relax only incidental assertions | Retain behavioral equivalence; exact column iteration order is incidental and may be relaxed if field-to-node correspondence remains verified. |
| [test_errors_identify_batch_stage_and_original_cause](../tests/test_execution_contract.py#L206) | Relax only incidental assertions | Retain failure propagation and context; exact message wording can be relaxed or consolidated with callback-exception coverage. |
| [test_callback_exception_reports_component_and_retains_cause](../tests/test_execution_contract.py#L235) | Consolidate without losing conditions | Consolidate with stage-failure diagnostics; preserving the cause matters, but requiring the helper's exact function name is incidental. |
| [test_enclosed_mass_retains_small_radius_precision](../tests/test_nfw_contract.py#L37) | Relax only incidental assertions | Keep cancellation coverage at scientifically used radii; the 1e-150 stress point and near-machine relative accuracy are candidates for a separate numerical stress suite, not a universal physics acceptance threshold. |
| [test_inverse_preserves_small_positive_radii](../tests/test_nfw_contract.py#L56) | Relax only incidental assertions | Keep small-radius positivity and inversion accuracy; extreme subnormal-scale coverage can be separated from model-domain validation with the supported radius range documented. |
| [test_lognormal_invalid_inputs](../tests/test_numerics.py#L35) | Consolidate without losing conditions | Retain as part of a consolidated quadrature-domain test; the supported scatter and order validation should remain covered. |
| [test_success_message_remains_a_diagnostic](../tests/test_odeint_controls.py#L269) | Consolidate without losing conditions | Consolidate with solver diagnostics or remove exact warning wording; numerical success and failure propagation are already tested separately. |
| [test_constant_grid_still_validates_solver_options](../tests/test_odeint_controls.py#L431) | Consolidate without losing conditions | Consolidate into invalid-controls coverage as a zero-duration parameter case; preserve the bypass regression. |
| [test_cache_reuse_and_stale_table_rejection_after_provider_change](../tests/test_picard.py#L167) | Relax only incidental assertions | Keep stale-state rejection and changed physics; exact cache-object reuse can be relaxed if equivalent safe caching is introduced. |
| [test_mass_unit_rescaling_preserves_fractional_evolution](../tests/test_picard.py#L343) | Relax only incidental assertions | Retain dimensional invariance; the extreme scale and near-machine tolerance may move to a numerical stress tier once realistic unit conversions have explicit coverage. |
| [test_nonreal_profile_and_window_coordinates_fail](../tests/test_profile_window_inputs.py#L25) | Consolidate without losing conditions | Retain complex-value protection; consolidate the shared input matrix across entry points without eliminating domain-specific checks. |
| [test_profile_parameters_do_not_accept_nonreal_values](../tests/test_profile_window_inputs.py#L55) | Consolidate without losing conditions | Can merge with test_profile_requires_finite_positive_parameters as one parameter-domain matrix while retaining all distinct invalid cases. |
| [test_logical_itamae_name_resolves_renamed_distribution](../tests/test_provenance.py#L157) | Consolidate without losing conditions | Consolidate with distribution-alias version lookup as one spelling matrix; the duplicated body adds no independent scientific coverage. |
| [test_distribution_alias_normalization_for_version_lookup](../tests/test_provenance.py#L213) | Consolidate without losing conditions | Consolidate with logical-name lookup while retaining all declared aliases. |
| [test_distribution_name_alias_uses_itamae_embedded_provenance](../tests/test_provenance.py#L269) | Consolidate without losing conditions | Consolidate with normalized-alias embedded-provenance coverage; retain a case for the canonical spelling. |
| [test_distribution_alias_normalization_for_embedded_provenance](../tests/test_provenance.py#L298) | Consolidate without losing conditions | Consolidate alias cases; exact module-load call assertions can be relaxed if durable source identity stays verified. |
| [test_npz_rejects_nonmapping_manifest](../tests/test_public_input_contracts.py#L206) | Consolidate without losing conditions | Consolidate with schema and duplicate-field tests as malformed-archive cases; keep the load-time boundary. |
| [test_background_density_overflow_is_explicit](../tests/test_public_input_contracts.py#L398) | Relax only incidental assertions | Keep overflow detection; move fantastically large epochs to a numerical stress tier if routine scientific CI is simplified. |
| [test_subnormal_nfw_target_has_a_finite_warning_free_bracket](../tests/test_public_input_contracts.py#L416) | Consolidate without losing conditions | Candidate for a separate stress suite or consolidation with small-radius inversion; it does not constrain the resolved subhalo domain by itself. |
| [test_shell_probabilities_normalize](../tests/test_spatial.py#L33) | Consolidate without losing conditions | Consolidate with radial-PDF normalization if one parametrized test still covers both point-density and finite-shell representations. |
| [test_native_units_reject_nonfinite](../tests/test_units.py#L7) | Consolidate without losing conditions | Merge the nonfinite case with public unit-boundary coverage while retaining unknown-unit and dimension cases. |
| [test_integrator_preserves_strict_tabulated_endpoints](../tests/test_variance.py#L169) | Consolidate without losing conditions | Consolidate this smoke case with the analytic sharp-k integral test while preserving strict endpoint enforcement. |

## Individual assessment of every test

### tests/test_backends.py

#### 1. `test_backend_config_validates_and_exports_metadata`

[Source](../tests/test_backends.py#L8) · **Relax only incidental assertions** · Role: guard.

**Purpose.** Record the cosmology and unit implementation used for a physical calculation.
**Pass.** Supported NumPy/native backends produce consistent metadata; selecting the unsupported JAX array backend raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep backend selection and metadata coverage, but the exact identifier prefix is a formatting choice that can be relaxed without weakening physical provenance.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 2. `test_backend_config_rejects_incomplete_objects`

[Source](../tests/test_backends.py#L29) · **Retain** · Role: guard.

**Purpose.** Fail before a calculation uses an object lacking cosmology or unit operations.
**Pass.** Objects missing the CosmologyBackend or UnitBackend interfaces raise TypeError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A missing background or unit operation must not become a partially computed catalogue.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_build_source_identity.py

#### 3. `test_archive_in_unrelated_repository_keeps_embedded_revision`

[Source](../tests/test_build_source_identity.py#L140) · **Retain** · Role: reproducibility.

**Purpose.** Keep a published catalogue attributable to the code inside its source archive.
**Pass.** An exported package inside an unrelated checkout resolves to its embedded source SHA rather than the enclosing repository SHA.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain one integration case per distinct build hook; otherwise a valid numerical result may cite the enclosing project's revision.
**Limit.** These checks establish source attribution, not the numerical or physical correctness of the packaged model.

#### 4. `test_unidentified_export_rejects_enclosing_repository`

[Source](../tests/test_build_source_identity.py#L164) · **Retain** · Role: reproducibility.

**Purpose.** Reject an export whose scientific source revision cannot be established.
**Pass.** A source export lacking embedded identity raises RuntimeError instead of using the surrounding repository's SHA.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This prevents false reproducibility claims; it is independent of physical accuracy.
**Limit.** These checks establish source attribution, not the numerical or physical correctness of the packaged model.

#### 5. `test_exact_source_checkout_uses_own_revision`

[Source](../tests/test_build_source_identity.py#L190) · **Retain** · Role: reproducibility.

**Purpose.** Identify the actual source checkout used to produce a catalogue.
**Pass.** A package at its own Git root resolves to that checkout's exact revision.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain an end-to-end provenance case for each packaging system.
**Limit.** These checks establish source attribution, not the numerical or physical correctness of the packaged model.

#### 6. `test_archive_rejects_conflicting_explicit_revision`

[Source](../tests/test_build_source_identity.py#L212) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Prevent a build override from relabelling a preserved scientific source archive.
**Pass.** A supplied revision that conflicts with the archive's embedded SHA raises RuntimeError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with matching-override and malformed-override cases into one parametrized provenance test; preserve all conflict scenarios.
**Limit.** These checks establish source attribution, not the numerical or physical correctness of the packaged model.

#### 7. `test_exact_checkout_rejects_conflicting_explicit_revision`

[Source](../tests/test_build_source_identity.py#L241) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Prevent an environment override from relabelling the checked-out calculation.
**Pass.** A supplied revision different from the checkout's SHA raises RuntimeError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate override cases while retaining both checkout and archive origins; they protect different identity sources.
**Limit.** These checks establish source attribution, not the numerical or physical correctness of the packaged model.

#### 8. `test_malformed_explicit_revision_does_not_silently_fall_back`

[Source](../tests/test_build_source_identity.py#L269) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Reject an unusable source identifier before a scientific package is built.
**Pass.** Malformed explicit revision values raise RuntimeError even when another source identity is available.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate malformed values into the common override matrix; rejecting ambiguous provenance remains required.
**Limit.** These checks establish source attribution, not the numerical or physical correctness of the packaged model.

#### 9. `test_matching_explicit_archive_revision_is_allowed`

[Source](../tests/test_build_source_identity.py#L297) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Permit an explicitly confirmed archive identity without changing its provenance.
**Pass.** An explicit SHA equal to the embedded archive SHA is accepted and recorded unchanged.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Merge this positive control with conflicting-override coverage so the rejection tests cannot pass by rejecting every build.
**Limit.** These checks establish source attribution, not the numerical or physical correctness of the packaged model.

### tests/test_cache_contract.py

#### 10. `test_spectrum_label_cannot_alias_different_table_content`

[Source](../tests/test_cache_contract.py#L17) · **Retain** · Role: reproducibility.

**Purpose.** Keep density variances tied to the actual input power spectrum.
**Pass.** Changing a power sample changes the identifier despite an equal human label; numerically equal tables with different byte order retain one identifier.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A stale variance table can change halo abundances silently, so content identity is scientifically necessary.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 11. `test_spectrum_identity_includes_extrapolation_policy`

[Source](../tests/test_cache_contract.py#L38) · **Retain** · Role: reproducibility.

**Purpose.** Distinguish spectra that make different assumptions outside their tabulated range.
**Pass.** The same power table with extrapolation enabled and disabled has different identities.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Outside-range power affects integrated small-scale variance and must not reuse an incompatible cache.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 12. `test_variance_identity_includes_derivative_resolution_but_not_batching`

[Source](../tests/test_cache_contract.py#L57) · **Relax only incidental assertions** · Role: reproducibility.

**Purpose.** Distinguish numerical controls that can change the EPS variance derivative.
**Pass.** Derivative step and sharp-k quadrature order change the identity; chunk size does not.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep sensitivity to accuracy controls; consolidate the chunk-size identity assertion with chunk-independence tests rather than requiring an identifier design.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 13. `test_cache_detects_finite_payload_corruption`

[Source](../tests/test_cache_contract.py#L84) · **Retain** · Role: reproducibility.

**Purpose.** Detect an altered mass grid or variance table even when all entries remain physical.
**Pass.** Changing one finite cached mass or variance value by ten percent causes the stored digest check to fail.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Range checks alone cannot detect a corrupted variance that biases accretion weights.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 14. `test_cache_rejects_undersized_grids`

[Source](../tests/test_cache_contract.py#L118) · **Retain** · Role: reproducibility.

**Purpose.** Reject a variance cache that cannot define a usable tabulated mass dependence.
**Pass.** Empty and one-node mass grids are rejected with ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A single point cannot support the interpolation or derivative needed by the population model.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 15. `test_cache_rejects_old_unverified_schema`

[Source](../tests/test_cache_contract.py#L142) · **Retain** · Role: reproducibility.

**Purpose.** Prevent an unverifiable historical cache from being reused as trusted variance data.
**Pass.** A version-1 archive without the current verified fields is rejected.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain while old cache files can be encountered; retire only with an explicit cache-format support policy, not because the test is nonphysical.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_calculation_provenance.py

#### 16. `test_calculation_schema_roundtrip_does_not_relabel_history`

[Source](../tests/test_calculation_provenance.py#L52) · **Retain** · Role: reproducibility.

**Purpose.** Keep old and current physical prescriptions distinguishable after saving a catalogue.
**Pass.** Historical metadata retain physics_mode, current metadata retain calculation_specification, and both survive NPZ round trips exactly.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Relabelling a historical result would make comparisons scientifically misleading.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 17. `test_new_provenance_rejects_conflicting_identity`

[Source](../tests/test_calculation_provenance.py#L100) · **Retain** · Role: reproducibility.

**Purpose.** Prevent auxiliary metadata from overriding the declared calculation prescription.
**Pass.** Conflicting physics_mode, specification, or provenance-schema fields raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A stored catalogue must not carry two incompatible claims about its physical calculation.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 18. `test_new_provenance_requires_a_specification`

[Source](../tests/test_calculation_provenance.py#L128) · **Retain** · Role: reproducibility.

**Purpose.** Require an identifiable physical prescription for newly generated catalogues.
**Pass.** Missing, blank, or non-string calculation specifications raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Preserve an unambiguous link between results and their model assumptions.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_catalog.py

#### 19. `test_catalog_factorized_weights_and_selection`

[Source](../tests/test_catalog.py#L52) · **Retain** · Role: guard.

**Purpose.** Interpret population weights as the product of abundance and survival factors.
**Pass.** The toy catalogue has final weights [1,1,0], weighted mass 3, and the requested two-row selection.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This checks the measure used by mass functions, not just array storage.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 20. `test_catalog_rejects_invalid_schema_and_weights`

[Source](../tests/test_catalog.py#L72) · **Retain** · Role: guard.

**Purpose.** Prevent malformed population measures from entering observable calculations.
**Pass.** Missing base weights, mismatched node shapes, negative weights, missing identities, and reserved-field overrides raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Misaligned or negative expected counts can silently corrupt observables; retain these input boundaries.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 21. `test_catalog_histogram_realization_and_concatenation`

[Source](../tests/test_catalog.py#L118) · **Retain** · Role: guard.

**Purpose.** Preserve expected subhalo counts when binning, realizing, or joining populations.
**Pass.** Toy histogram counts are [1,1]; a seeded realization has the expected fields and one dimension; joining duplicates preserves six nodes and their weights.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain histogram and concatenation checks; the realization part is only a shape smoke test and is not a Poisson-distribution validation.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 22. `test_poisson_realization_requires_explicit_generator`

[Source](../tests/test_catalog.py#L146) · **Retain** · Role: guard.

**Purpose.** Make stochastic subhalo realizations use a caller-controlled random stream.
**Pass.** Passing no NumPy Generator raises TypeError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Explicit random-state ownership is required to reproduce Monte Carlo realizations.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 23. `test_catalog_npz_round_trip_without_pickle`

[Source](../tests/test_catalog.py#L163) · **Retain** · Role: guard.

**Purpose.** Preserve masses, weights, and scientific metadata exactly in a saved catalogue.
**Pass.** All columns, weight factors, names, and metadata are identical after NPZ save and reload.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Exact equality is appropriate for lossless storage; it is not a demand for machine-precision physical predictions.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 24. `test_catalog_npz_rejects_wrong_archive_schema`

[Source](../tests/test_catalog.py#L196) · **Retain** · Role: guard.

**Purpose.** Reject saved data whose fields cannot be interpreted as a supported catalogue.
**Pass.** An archive declaring an unknown schema raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Without a known schema, units and population meaning cannot be trusted.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_cosmology.py

#### 25. `test_native_cosmology_basic_relations`

[Source](../tests/test_cosmology.py#L7) · **Retain** · Role: numerical.

**Purpose.** Check the present-day normalization and sign of the native cosmological clock.
**Pass.** H(0)=67.4 km/s/Mpc, D(0)=1, cosmic times at z=0 and 1 are positive, and the z=1 lookback time is positive.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These are basic background identities; they complement derivative and limiting-cosmology tests.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 26. `test_colossus_backend_matches_order_of_magnitude`

[Source](../tests/test_cosmology.py#L28) · **Retain** · Role: numerical.

**Purpose.** Detect gross background mismatches against the optional Colossus implementation.
**Pass.** At z=0,1,3 the two Planck-like backgrounds agree within 3 percent in H and 4 percent in growth.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain as a broad independent smoke check; parameter and growth-prescription differences preclude interpreting these loose tolerances as precision calibration.
**Limit.** Only H and growth are compared, not critical/matter density. The current Colossus wrapper multiplies its Msun*h^2/kpc^3 density by 1e9 without applying h^2; this test cannot detect that canonical-unit mismatch.

### tests/test_evolution.py

#### 27. `test_ode_exponential_decay`

[Source](../tests/test_evolution.py#L7) · **Retain** · Role: numerical.

**Purpose.** Recover a known loss history with the generic evolution integrator.
**Pass.** For dy/dt=-2y and y(0)=1, 21 outputs over t=0..2 agree with exp(-2t) to relative tolerance 1e-6.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** An independent analytic solution detects integration-sign and output-grid errors.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 28. `test_shanks_geometric_sequence`

[Source](../tests/test_evolution.py#L26) · **Retain** · Role: numerical.

**Purpose.** Check the sequence-acceleration formula used by the historical stripping approximation.
**Pass.** Applying Shanks to 1,1.5,1.75 recovers the geometric-series limit 2.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain while Shanks remains a supported solver; this algebraic check does not validate its use for a particular tidal history.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 29. `test_perturbative_solver_constant_rhs_and_exponential_sequence`

[Source](../tests/test_evolution.py#L42) · **Retain** · Role: numerical.

**Purpose.** Recover analytic constant-loss and exponential histories through Picard iteration.
**Pass.** The constant RHS yields 2-t; order-eight iteration follows exp(-t) with rtol=1e-6 and atol=1e-7 on 401 nodes.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These independent solutions validate the numerical integration component rather than a frozen implementation output.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 30. `test_perturbative_solver_validates_grid_and_rhs`

[Source](../tests/test_evolution.py#L66) · **Retain** · Role: numerical.

**Purpose.** Stop a numerical evolution whose time order or derivative is undefined.
**Pass.** A nonmonotonic time grid or NaN right-hand side raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Returning a plausible array after such failures would contaminate subsequent structure calculations.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_execution.py

#### 31. `test_execution_reuses_canonical_variance_protocol`

[Source](../tests/test_execution.py#L41) · **Retire or replace, conditionally** · Role: guard.

**Purpose.** Track the compatibility alias for the variance interface used by population execution.
**Pass.** The re-exported VarianceModel is the same Python object as the canonical protocol.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retirement candidate: this is object identity only; test_execution_composition plus real variance-adapter tests exercise usable behavior. Retain only if this exact import identity is a supported API promise.
**Limit.** This is an API or output-compatibility check; it does not test a physical equation or numerical accuracy.

#### 32. `test_population_pipeline_transports_stages_and_factorized_weights`

[Source](../tests/test_execution.py#L58) · **Retain** · Role: guard.

**Purpose.** Keep each accretion node aligned with its evolved mass and abundance factors.
**Pass.** The two toy slices yield bound masses [1.1,2.2,4.2,6.4], unchanged base factors, and final CDM weights [0.14,0.24,0.14,0.24].
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Node or factor misalignment directly changes mass functions and must be detected independently of a full halo model.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 33. `test_population_pipeline_rejects_misaligned_contexts`

[Source](../tests/test_execution.py#L95) · **Retain** · Role: guard.

**Purpose.** Prevent redshift or host context from being assigned to the wrong accretion slice.
**Pass.** A context list with a different length from the batch list raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A mismatched context can evolve an otherwise valid node with the wrong host history.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_execution_composition.py

#### 34. `test_population_components_delegate_to_canonical_pipeline`

[Source](../tests/test_execution_composition.py#L242) · **Retire or replace, conditionally** · Role: guard.

**Purpose.** Check the concrete pipeline class returned by component composition.
**Pass.** The composed object is an instance of PopulationPipeline.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retirement candidate: the adjacent callback-versus-component execution test verifies the behavior that matters. Requiring this concrete class can forbid equivalent implementations without improving scientific reliability.
**Limit.** This is an API or output-compatibility check; it does not test a physical equation or numerical accuracy.

#### 35. `test_population_components_match_callback_execution_and_catalogs`

[Source](../tests/test_execution_composition.py#L265) · **Relax only incidental assertions** · Role: guard.

**Purpose.** Preserve the physical population when the same stages are supplied as objects or callbacks.
**Pass.** Both routes return exactly the same toy columns, survival views, weight factors, and strict-view catalogue.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain behavioral equivalence; exact column iteration order is incidental and may be relaxed if field-to-node correspondence remains verified.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_execution_contract.py

#### 36. `test_mixed_metadata_rejected_before_callbacks`

[Source](../tests/test_execution_contract.py#L87) · **Retain** · Role: guard.

**Purpose.** Prevent merging populations with different cosmologies, units, revisions, or prescriptions.
**Pass.** Each incompatible pair raises ValueError before any initialization callback runs.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Mixing individually valid populations can produce an uninterpretable mass function; reject before physical evolution.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 37. `test_concatenation_checks_nested_metadata_independent_of_mapping_order`

[Source](../tests/test_execution_contract.py#L121) · **Retain** · Role: guard.

**Purpose.** Distinguish physical metadata content from dictionary insertion order.
**Pass.** Equivalent nested cosmology metadata concatenate to eight nodes; changed or absent cosmology data are rejected.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain content comparison to prevent both unsafe mixing and unnecessary rejection of equivalent descriptions.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 38. `test_survival_rejects_truthy_non_masks`

[Source](../tests/test_execution_contract.py#L148) · **Retain** · Role: guard.

**Purpose.** Keep a disruption decision distinct from a fractional population weight.
**Pass.** Fractional, negative, or text-valued survival masks raise an error.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Implicit truth conversion would count disrupted or fractional-weight nodes incorrectly.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 39. `test_boolean_and_binary_masks_retain_selection`

[Source](../tests/test_execution_contract.py#L172) · **Retain** · Role: guard.

**Purpose.** Preserve the same disrupted and surviving nodes for supported mask encodings.
**Pass.** Boolean and exact 0/1 masks are transported without changing selected nodes.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This positive control guards against over-strict validation that rejects physically equivalent masks.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 40. `test_errors_identify_batch_stage_and_original_cause`

[Source](../tests/test_execution_contract.py#L206) · **Relax only incidental assertions** · Role: guard.

**Purpose.** Make a failed structure or survival stage traceable to its physical slice.
**Pass.** Malformed stage outputs raise ValueError naming the batch and stage while retaining the original exception.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain failure propagation and context; exact message wording can be relaxed or consolidated with callback-exception coverage.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 41. `test_callback_exception_reports_component_and_retains_cause`

[Source](../tests/test_execution_contract.py#L235) · **Consolidate without losing conditions** · Role: guard.

**Purpose.** Preserve a model-domain failure rather than returning a partial population.
**Pass.** An ArithmeticError from the evolver becomes a contextual error and remains the chained cause.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with stage-failure diagnostics; preserving the cause matters, but requiring the helper's exact function name is incidental.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 42. `test_column_order_is_not_a_schema_difference`

[Source](../tests/test_execution_contract.py#L274) · **Retain** · Role: guard.

**Purpose.** Join physically identical catalogues despite different field insertion order.
**Pass.** Two toy slices with reversed dictionary order retain the expected accreted-mass sequence.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This protects usability without relaxing field-name or node-alignment checks.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 43. `test_optional_weight_mismatch_precedes_callbacks`

[Source](../tests/test_execution_contract.py#L330) · **Retain** · Role: guard.

**Purpose.** Reject a population mix where only some slices include host-history weighting.
**Pass.** Presence mismatch of weight_host_history raises ValueError before callbacks run.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Otherwise the combined effective number of subhaloes has inconsistent normalization.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 44. `test_scalar_batch_has_an_explicit_population_axis_error`

[Source](../tests/test_execution_contract.py#L352) · **Retain** · Role: guard.

**Purpose.** Reject a batch that lacks a node axis needed for population alignment.
**Pass.** Replacing every one-node array by a scalar raises a population-axis error.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep this guard until scalar batches have a defined, tested broadcast convention.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 45. `test_empty_batches_zero_survivors_and_deterministic_partition`

[Source](../tests/test_execution_contract.py#L385) · **Retain** · Role: guard.

**Purpose.** Keep population results independent of how accretion nodes are partitioned.
**Pass.** Empty slices and zero survivors are handled, no-batch input fails, and partitioned versus whole populations agree exactly including a fixed-seed realization.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This is deterministic transport and random-stream reproducibility, not statistical calibration of Poisson sampling.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_halo.py

#### 46. `test_nfw_inverse_roundtrip`

[Source](../tests/test_halo.py#L9) · **Retain** · Role: numerical.

**Purpose.** Recover the dimensionless truncation radius from its NFW enclosed mass.
**Pass.** For generated radii 0..1e4, inversion recovers the radius with relative and absolute tolerance 1e-9.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain broad radius coverage, but interpret this as an internal inverse check; the Decimal-reference tests provide the independent equation check.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 47. `test_nfw_profile_mass_is_monotonic`

[Source](../tests/test_halo.py#L32) · **Retain** · Role: numerical.

**Purpose.** Keep an NFW halo's enclosed mass increasing with physical radius.
**Pass.** The enclosed mass increases strictly over 100 radii from 1e-5 to 1 Mpc for the fixed positive profile.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Monotonic enclosed mass is a physical invariant distinct from matching a reference implementation.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_measure.py

#### 48. `test_build_accretion_batch_broadcasts_and_preserves_weights`

[Source](../tests/test_measure.py#L6) · **Retain** · Role: numerical.

**Purpose.** Broadcast a shared epoch and concentration without changing accretion-node weights.
**Pass.** Two masses retain their virial masses while scalar concentration and host-history factors broadcast to the two-node shape.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A population measure must keep each physical state paired with its quadrature weight.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_nfw_contract.py

#### 49. `test_enclosed_mass_retains_small_radius_precision`

[Source](../tests/test_nfw_contract.py#L37) · **Relax only incidental assertions** · Role: numerical.

**Purpose.** Resolve the NFW central mass without cancellation in log(1+x)-x/(1+x).
**Pass.** The dimensionless mass agrees with 340-digit Decimal evaluation at radii 1e-150..1e6 with rtol=5e-14 and zero absolute tolerance.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep cancellation coverage at scientifically used radii; the 1e-150 stress point and near-machine relative accuracy are candidates for a separate numerical stress suite, not a universal physics acceptance threshold.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 50. `test_inverse_preserves_small_positive_radii`

[Source](../tests/test_nfw_contract.py#L56) · **Relax only incidental assertions** · Role: numerical.

**Purpose.** Avoid turning a small positive bound mass into a zero truncation radius.
**Pass.** Inverting independent Decimal masses over radii 1e-150..1e6 recovers radii to rtol=5e-12; zero and empty inputs keep their meanings.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep small-radius positivity and inversion accuracy; extreme subnormal-scale coverage can be separated from model-domain validation with the supported radius range documented.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 51. `test_nfw_mass_rejects_nonfinite_or_negative_inputs`

[Source](../tests/test_nfw_contract.py#L79) · **Retain** · Role: numerical.

**Purpose.** Reject a radius or enclosed-mass target that has no physical NFW interpretation.
**Pass.** NaN, infinities and negative values raise ValueError in both mass and inverse functions.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Silent clipping would conceal a broken evolved halo state.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 52. `test_profile_requires_finite_positive_parameters`

[Source](../tests/test_nfw_contract.py#L106) · **Retain** · Role: numerical.

**Purpose.** Require a physical NFW scale radius and scale density.
**Pass.** Zero, negative or nonfinite r_s or rho_s raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These parameters determine masses, velocities and annihilation luminosity, so invalid values must fail early.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 53. `test_central_potential_uses_its_limit_without_evaluating_zero_over_zero`

[Source](../tests/test_nfw_contract.py#L130) · **Retain** · Role: numerical.

**Purpose.** Evaluate the finite NFW central potential by its analytic limit.
**Pass.** At r=0 the potential equals -4*pi*G*rho_s*r_s**2 to rtol=2e-15 without RuntimeWarning; negative radius fails.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** The central-limit identity is exact for the stated profile and avoids NaNs in orbit calculations; retain it.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_numerics.py

#### 54. `test_lognormal_quadrature_is_positive_and_normalized`

[Source](../tests/test_numerics.py#L9) · **Retain** · Role: numerical.

**Purpose.** Keep concentration-scatter quadrature a normalized positive population measure.
**Pass.** For generated positive medians and scatter 0..0.5, seven-point Hermite nodes are positive and weights sum to one.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Normalization is required so scatter averaging does not create or destroy expected subhalo counts.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 55. `test_lognormal_invalid_inputs`

[Source](../tests/test_numerics.py#L35) · **Consolidate without losing conditions** · Role: numerical.

**Purpose.** Reject a lognormal concentration distribution with nonpositive median.
**Pass.** A zero median raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain as part of a consolidated quadrature-domain test; the supported scatter and order validation should remain covered.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 56. `test_population_grids_are_bounded_and_explicit`

[Source](../tests/test_numerics.py#L52) · **Retain** · Role: numerical.

**Purpose.** Keep accretion integration grids within the requested physical support.
**Pass.** The logarithmic mass grid and redshift cells match the declared endpoints; optional final-cell omission works; reversed or nonpositive bounds fail.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Grid endpoints control which accretion epochs and masses contribute to abundance, so preserve this coverage.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_ode_contract.py

#### 57. `test_odeint_preserves_scipy_defaults_and_time_order`

[Source](../tests/test_ode_contract.py#L11) · **Retain** · Role: numerical.

**Purpose.** Preserve historical odeint histories in both time and redshift directions.
**Pass.** Forward and reversed grids match direct SciPy odeint exactly and exponential analytic solutions to rtol=1e-7.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep the analytic and sign checks; bitwise parity is a wrapper-compatibility promise, not a scientific accuracy requirement for alternative integrators.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 58. `test_invalid_rhs_is_rejected`

[Source](../tests/test_ode_contract.py#L47) · **Retain** · Role: numerical.

**Purpose.** Stop an evolution with undefined or dimensionally inconsistent derivatives.
**Pass.** RK45 and odeint reject NaN derivatives or a derivative vector of the wrong size.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A successful return with corrupted derivatives would invalidate all downstream mass and structure values.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 59. `test_invalid_inputs_are_rejected`

[Source](../tests/test_ode_contract.py#L85) · **Retain** · Role: numerical.

**Purpose.** Reject evolution inputs that cannot define a finite real initial-value problem.
**Pass.** Both solver paths reject nonfinite or empty initial states, invalid time grids, negative/NaN tolerances and complex states.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These failures prevent silent data loss or misinterpreted numerical accuracy settings.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_odeint_controls.py

#### 60. `test_selected_controls_preserve_direct_odeint_values`

[Source](../tests/test_odeint_controls.py#L10) · **Retain** · Role: numerical.

**Purpose.** Apply the caller's time-step and work limits to the actual physical integration.
**Pass.** Forward and reversed forced-decay solutions exactly match direct SciPy with the same hmax, h0, mxstep and tolerances, leaving the options unchanged.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain control forwarding; exact SciPy equality is suitable for a wrapper but must not be demanded of different numerical methods.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 61. `test_stiff_jacobian_uses_canonical_time_first_signature`

[Source](../tests/test_odeint_controls.py#L79) · **Retain** · Role: numerical.

**Purpose.** Use the supplied Jacobian for a system with widely separated decay time-scales.
**Pass.** Decay rates 1 and 1000 with a time-first Jacobian reproduce direct odeint exactly.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This detects argument-order mistakes that can corrupt stiff evolution even when nonstiff examples pass.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 62. `test_invalid_controls_fail_before_rhs`

[Source](../tests/test_odeint_controls.py#L188) · **Retain** · Role: numerical.

**Purpose.** Reject ambiguous integration settings before invoking the physical model.
**Pass.** Invalid step bounds, orders, work counts, Jacobian settings and unsupported overrides raise ValueError before the RHS runs.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain explicit rejection; unsupported controls must not be silently ignored.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 63. `test_odeint_controls_are_not_silently_used_by_another_solver`

[Source](../tests/test_odeint_controls.py#L232) · **Retain** · Role: numerical.

**Purpose.** Prevent a user from believing LSODA-specific accuracy controls affected another solver.
**Pass.** Passing odeint_options without selecting odeint raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This guards the interpretation of numerical convergence studies.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 64. `test_failed_internal_step_limit_is_not_reported_as_success`

[Source](../tests/test_odeint_controls.py#L249) · **Retain** · Role: numerical.

**Purpose.** Do not return an uncompleted physical trajectory as a valid solution.
**Pass.** For mxstep=1 on a nontrivial decay interval, the wrapper raises RuntimeError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Solver failure is scientifically different from a finite-looking partial result.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 65. `test_success_message_remains_a_diagnostic`

[Source](../tests/test_odeint_controls.py#L269) · **Consolidate without losing conditions** · Role: numerical.

**Purpose.** Distinguish an informational LSODA success warning from failed integration.
**Pass.** printmessg emits the expected success warning while the final decay value agrees with exp(-1) to rtol=1e-7.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with solver diagnostics or remove exact warning wording; numerical success and failure propagation are already tested separately.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 66. `test_banded_jacobian_shape_and_orientation`

[Source](../tests/test_odeint_controls.py#L291) · **Retain** · Role: numerical.

**Purpose.** Interpret row- and column-oriented banded Jacobians consistently in stiff evolution.
**Pass.** Both diagonal-band orientations reproduce direct odeint for decay rates 1 and 1000.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain while the public solver accepts banded Jacobians; shape mistakes otherwise alter derivatives silently.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 67. `test_invalid_jacobian_is_not_returned_as_a_solution`

[Source](../tests/test_odeint_controls.py#L373) · **Retain** · Role: numerical.

**Purpose.** Reject an invalid linearization of the physical evolution equations.
**Pass.** Wrong-shape, NaN or complex Jacobians raise ValueError rather than returning a trajectory.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This is a necessary boundary for user-supplied stiff-system derivatives.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 68. `test_explicit_repeated_outputs_preserve_scipy_shape_and_values`

[Source](../tests/test_odeint_controls.py#L404) · **Retain** · Role: numerical.

**Purpose.** Allow repeated observation epochs without dropping or moving physical samples.
**Pass.** Descending, constant and singleton time grids match SciPy shapes and values when repeats are explicitly enabled.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Repeated epochs occur in stored SIDM histories; retain their correspondence to requested outputs.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 69. `test_constant_grid_still_validates_solver_options`

[Source](../tests/test_odeint_controls.py#L431) · **Consolidate without losing conditions** · Role: numerical.

**Purpose.** Apply the same numerical-control validation to zero-duration evolution.
**Pass.** A negative hmax is rejected even when every requested epoch is identical.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate into invalid-controls coverage as a zero-duration parameter case; preserve the bypass regression.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_picard.py

#### 70. `test_endpoint_and_history_match_independent_analytic_solution`

[Source](../tests/test_picard.py#L87) · **Retain** · Role: numerical.

**Purpose.** Recover the analytic tidal mass history for constant host and loss coefficients.
**Pass.** For exponents 0,0.04,-0.01, Picard masses agree within rtol=5e-4 and DOP853 within 5e-10, including repeated outputs and an omitted accretion-time output.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This independent solution checks integration direction, mass scaling and history interpolation without using a frozen Picard output.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 71. `test_grid_refinement_reduces_time_dependent_coefficient_error`

[Source](../tests/test_picard.py#L118) · **Retain** · Role: numerical.

**Purpose.** Separate redshift-quadrature error from nonlinear Picard convergence.
**Pass.** For rate=0.1*z**2 and zero exponent, refinement from 33 to 129 nodes reduces analytic endpoint error by more than tenfold and below 2e-4.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain this controlled convergence test; it does not imply every physical model improves monotonically at every queried point.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 72. `test_endpoint_broadcasting_and_history_shape`

[Source](../tests/test_picard.py#L142) · **Retain** · Role: numerical.

**Purpose.** Keep accretion masses and requested epochs aligned in vectorized stripping queries.
**Pass.** Endpoint broadcasting and time-leading history shapes are correct; zero-duration histories preserve mass; empty inputs preserve empty axes.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A shape error can mix halo populations even when each individual computed mass is finite.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 73. `test_cache_reuse_and_stale_table_rejection_after_provider_change`

[Source](../tests/test_picard.py#L167) · **Relax only incidental assertions** · Role: numerical.

**Purpose.** Recompute stripping when the physical mass-loss rate changes.
**Pass.** Unchanged state reuses a table; doubling the rate rejects the old table, builds a replacement and lowers final mass.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep stale-state rejection and changed physics; exact cache-object reuse can be relaxed if equivalent safe caching is introduced.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 74. `test_domain_fallback_is_explicit_and_recorded`

[Source](../tests/test_picard.py#L195) · **Retain** · Role: numerical.

**Purpose.** Use a direct physical integration outside the caller-validated Picard domain.
**Pass.** A z_acc=8 query outside the declared domain warns, matches the analytic solution to 5e-10, and records one fallback covering two masses in JSON-safe diagnostics.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This prevents unvalidated table extrapolation from being presented as a verified approximation.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 75. `test_unconverged_iteration_uses_direct_reference`

[Source](../tests/test_picard.py#L226) · **Retain** · Role: numerical.

**Purpose.** Avoid accepting an unconverged nonlinear mass-loss approximation.
**Pass.** A deliberately difficult one-iteration case warns, returns the analytic reference to 5e-10, and counts the fallback mass.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** The fallback is a numerical safeguard; it is not evidence that the fixed Picard settings converge everywhere.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 76. `test_invalid_endpoint_options_fail`

[Source](../tests/test_picard.py#L262) · **Retain** · Role: numerical.

**Purpose.** Reject resolution or interpolation settings that cannot define a valid stripping table.
**Pass.** Zero iterations, fractional node counts, Boolean counts, zero tolerance and unknown interpolation each raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Unsupported numerical settings must not produce an apparently successful scientific catalogue.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 77. `test_bad_inputs_fail_without_clipping_or_extrapolation`

[Source](../tests/test_picard.py#L287) · **Retain** · Role: numerical.

**Purpose.** Keep invalid masses and reversed cosmic evolution out of the stripping result.
**Pass.** Nonpositive mass, final redshift above accretion, NaN epochs, increasing history grids, and direct table queries beyond bounds are rejected.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Do not silently replace invalid physical states with clipped masses or extrapolated table values.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 78. `test_invalid_physical_callback_fails_in_every_path`

[Source](../tests/test_picard.py#L315) · **Retain** · Role: numerical.

**Purpose.** Apply the same physical coefficient checks to Picard and direct integration.
**Pass.** Negative or NaN stripping rates raise ValueError for endpoint, history and direct paths.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Otherwise choosing a fallback could silently change which physical inputs are accepted.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 79. `test_mass_unit_rescaling_preserves_fractional_evolution`

[Source](../tests/test_picard.py#L343) · **Relax only incidental assertions** · Role: numerical.

**Purpose.** Make stripping depend on the subhalo-to-host mass ratio rather than a chosen mass unit.
**Pass.** Scaling host and subhalo masses together by 1e-30 preserves final masses after rescaling to rtol=5e-14.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain dimensional invariance; the extreme scale and near-machine tolerance may move to a numerical stress tier once realistic unit conversions have explicit coverage.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_power.py

#### 80. `test_tabulated_power_spectrum_log_interpolation_and_provenance`

[Source](../tests/test_power.py#L13) · **Retain** · Role: numerical.

**Purpose.** Recover a power-law spectrum without silently extending its support.
**Pass.** A 7*k**-2 table reproduces queried powers, exposes its source label and protocol, and rejects a point below k_min.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain analytic interpolation and support checks; source-label text alone is not provenance completeness.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 81. `test_transfer_modified_spectrum_takes_explicit_power_ratio`

[Source](../tests/test_power.py#L39) · **Retain** · Role: numerical.

**Purpose.** Apply a power suppression once, without confusing transfer amplitude with its square.
**Pass.** Multiplying a unit base spectrum by exp(-k**2) yields that power ratio and retains its identifier.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This distinguishes P_modified/P_base from the transfer amplitude T and prevents a systematic suppression error.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 82. `test_linear_power_interpolation_preserves_physical_zero`

[Source](../tests/test_power.py#L65) · **Retain** · Role: numerical.

**Purpose.** Represent a genuine spectral node without imposing a positive floor.
**Pass.** Linear interpolation preserves a zero tabulated power; log-log interpolation rejects that zero.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain for oscillatory FDM spectra, where replacing zero power changes the variance.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 83. `test_window_limits_and_protocols`

[Source](../tests/test_power.py#L90) · **Retain** · Role: numerical.

**Purpose.** Use correctly normalized smoothing windows at long wavelength.
**Pass.** Top-hat tends to one at x=0 and 1e-6; sharp-k gives [1,1,1,0] at x=[0,1e-6,1,2].
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** The filter normalization sets density variance and halo abundance; retain analytic boundary checks.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 84. `test_tabulated_power_rejects_invalid_tables`

[Source](../tests/test_power.py#L122) · **Retain** · Role: numerical.

**Purpose.** Reject a spectrum that cannot support a well-defined log interpolation.
**Pass.** Duplicate k, nonpositive power for log interpolation, and NaN k each raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These guards prevent undefined or misleading variance integrals.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_power_value_contract.py

#### 85. `test_tabulated_power_rejects_non_real_numeric_values`

[Source](../tests/test_power_value_contract.py#L11) · **Retain** · Role: guard.

**Purpose.** Prevent complex or nonnumeric spectrum data from being silently cast to real power.
**Pass.** Complex arrays, strings and Boolean values are rejected for table coordinates, power and queries.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain rejection of lossy casts; Boolean/string convenience acceptance would require an explicit input policy.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 86. `test_transfer_rejects_lossy_or_misaligned_ratio`

[Source](../tests/test_power_value_contract.py#L47) · **Retain** · Role: guard.

**Purpose.** Keep the suppression factor aligned with the queried wavenumbers.
**Pass.** Complex, text, Boolean or wrong-shape power ratios raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Broadcasting the wrong physical axis can produce a finite but incorrect spectrum.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 87. `test_transfer_rejects_nonfinite_product_without_warning`

[Source](../tests/test_power_value_contract.py#L72) · **Retain** · Role: guard.

**Purpose.** Report spectral overflow as an invalid physical result.
**Pass.** Overflow of finite base power times a finite ratio raises ValueError under strict floating-point handling.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A nonfinite power product must not enter density-variance quadrature.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 88. `test_scalar_transfer_keeps_query_shape`

[Source](../tests/test_power_value_contract.py#L92) · **Retain** · Role: guard.

**Purpose.** Apply a scale-independent power suppression over any query shape.
**Pass.** A scalar factor of one half preserves a two-dimensional query and returns half the base spectrum exactly.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This positive control avoids rejecting legitimate constant suppression when enforcing callback shape rules.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_profile_window_inputs.py

#### 89. `test_nonreal_profile_and_window_coordinates_fail`

[Source](../tests/test_profile_window_inputs.py#L25) · **Consolidate without losing conditions** · Role: guard.

**Purpose.** Prevent lossy conversion of physical radii and filter arguments.
**Pass.** Complex, text and Boolean coordinates raise ValueError without warnings in NFW mass, inverse, density, potential and both windows.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain complex-value protection; consolidate the shared input matrix across entry points without eliminating domain-specific checks.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 90. `test_profile_parameters_do_not_accept_nonreal_values`

[Source](../tests/test_profile_window_inputs.py#L55) · **Consolidate without losing conditions** · Role: guard.

**Purpose.** Reject non-real NFW scale radius and density before model construction.
**Pass.** Complex, text and Boolean r_s or rho_s raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Can merge with test_profile_requires_finite_positive_parameters as one parameter-domain matrix while retaining all distinct invalid cases.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_provenance.py

#### 91. `test_src_package_provenance_uses_its_own_repository`

[Source](../tests/test_provenance.py#L23) · **Retain** · Role: reproducibility.

**Purpose.** Locate the source checkout that actually owns a core or variant calculation.
**Pass.** Each supported src-layout module resolves to its own temporary project root.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain runtime source attribution; this exercises a different path from build-hook archive provenance.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 92. `test_src_export_does_not_inherit_an_unrelated_parent_repository`

[Source](../tests/test_provenance.py#L52) · **Retain** · Role: reproducibility.

**Purpose.** Avoid attributing an unpacked scientific package to an unrelated enclosing project.
**Pass.** An exported src-layout package under another Git repository has no source-checkout root.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Wrong repository attribution can make a result irreproducible despite valid numerical arrays.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 93. `test_build_migration_metadata_includes_standard_keys`

[Source](../tests/test_provenance.py#L80) · **Retain** · Role: reproducibility.

**Purpose.** Record the sources and physical prescription needed to interpret a historical migration result.
**Pass.** Metadata contain all shared keys, exact resolved source SHAs, schema identity and the supplied extra parameter.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain as a compatibility check for reading/reproducing historical results.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 94. `test_build_migration_metadata_rejects_standard_field_override`

[Source](../tests/test_provenance.py#L127) · **Retain** · Role: reproducibility.

**Purpose.** Prevent variant extras from changing the stated historical physical mode.
**Pass.** An extra physics_mode conflicting with the declared mode raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Preserve the boundary between result interpretation and arbitrary annotations.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 95. `test_logical_itamae_name_resolves_renamed_distribution`

[Source](../tests/test_provenance.py#L157) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Resolve the core's logical name to the installed distribution used for provenance.
**Pass.** Version lookup for itamae queries sashimi-itamae and returns the supplied version.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with distribution-alias version lookup as one spelling matrix; the duplicated body adds no independent scientific coverage.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 96. `test_distribution_alias_normalization_for_version_lookup`

[Source](../tests/test_provenance.py#L213) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Treat supported distribution-name spellings as one installed core version.
**Pass.** Case, hyphen, underscore and dot spellings all query the canonical distribution and return one version.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with logical-name lookup while retaining all declared aliases.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 97. `test_distribution_name_alias_uses_itamae_embedded_provenance`

[Source](../tests/test_provenance.py#L269) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Use the same embedded scientific source for logical and distribution names.
**Pass.** The public distribution alias resolves to the expected embedded SHA.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with normalized-alias embedded-provenance coverage; retain a case for the canonical spelling.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 98. `test_distribution_alias_normalization_for_embedded_provenance`

[Source](../tests/test_provenance.py#L298) · **Consolidate without losing conditions** · Role: reproducibility.

**Purpose.** Keep source identity stable under supported package-name spellings.
**Pass.** All aliases load itamae._build_provenance and return the expected SHA.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate alias cases; exact module-load call assertions can be relaxed if durable source identity stays verified.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 99. `test_embedded_revision_takes_precedence_over_ambient_environment`

[Source](../tests/test_provenance.py#L353) · **Retain** · Role: reproducibility.

**Purpose.** Protect installed scientific provenance from unrelated shell settings.
**Pass.** A conflicting environment SHA cannot replace the wheel's embedded source SHA.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain: published result provenance must refer to installed code, not a mutable environment label.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 100. `test_recognized_distribution_cannot_finish_with_unknown_revision`

[Source](../tests/test_provenance.py#L380) · **Retain** · Role: reproducibility.

**Purpose.** Refuse to record an unidentified installed family calculation as reproducible.
**Pass.** Recognized core aliases raise RuntimeError when no durable source revision is available.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A successful numerical calculation and a reproducible one are different; retain this provenance boundary.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

#### 101. `test_source_revision_does_not_walk_into_an_outer_repository_from_venv`

[Source](../tests/test_provenance.py#L411) · **Retain** · Role: reproducibility.

**Purpose.** Avoid assigning an application's revision to a dependency installed in its environment.
**Pass.** An unidentified package inside another project's virtual environment returns unknown without invoking ambient Git.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain because this common installation layout differs from a directly unpacked archive.
**Limit.** Source attribution and format compatibility do not validate the physical prescription or numerical accuracy.

### tests/test_public_input_contracts.py

#### 102. `test_direct_catalog_rejects_nonreal_or_nonfinite_columns`

[Source](../tests/test_public_input_contracts.py#L47) · **Retain** · Role: guard.

**Purpose.** Protect physical catalogue columns even when callers bypass population execution.
**Pass.** Direct construction rejects NaN, infinity, complex and text-valued columns with ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Public catalogue construction is an independent entry point; pipeline validation alone is insufficient.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 103. `test_direct_catalog_rejects_complex_weights_before_cast`

[Source](../tests/test_public_input_contracts.py#L71) · **Retain** · Role: guard.

**Purpose.** Prevent complex expected counts from being silently converted to real weights.
**Pass.** Direct construction with complex weight_base raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Loss of the imaginary component could hide a failed upstream calculation.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 104. `test_factor_product_overflow_is_an_error`

[Source](../tests/test_public_input_contracts.py#L88) · **Retain** · Role: guard.

**Purpose.** Reject a nonfinite final abundance even when each weight factor is finite.
**Pass.** Multiplying two factors of 1e200 raises ValueError instead of storing infinite weight.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Finite individual quadrature factors do not guarantee a finite physical measure.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 105. `test_weighted_sum_rejects_invalid_values`

[Source](../tests/test_public_input_contracts.py#L109) · **Retain** · Role: guard.

**Purpose.** Reject undefined observables before applying a population measure.
**Pass.** NaN, infinity or complex weighted values raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A single invalid contribution must not silently corrupt an integrated observable.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 106. `test_weighted_sum_overflow_is_an_error`

[Source](../tests/test_public_input_contracts.py#L132) · **Retain** · Role: guard.

**Purpose.** Report overflow of an integrated observable explicitly.
**Pass.** Summing two 1e308 contributions raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain failure detection; the extreme values are controlled stress inputs, not a claimed physical regime.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 107. `test_empty_density_histogram_has_no_silent_nan`

[Source](../tests/test_public_input_contracts.py#L149) · **Retain** · Role: guard.

**Purpose.** Distinguish zero subhaloes from an undefined normalized probability density.
**Pass.** An empty population gives zero count bins; requesting density normalization raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** An empty mass function is physical, but normalizing it to unit probability is undefined.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 108. `test_npz_rejects_duplicate_names_in_manifest`

[Source](../tests/test_public_input_contracts.py#L170) · **Retain** · Role: guard.

**Purpose.** Prevent saved catalogue fields from ambiguously overwriting one another.
**Pass.** An NPZ manifest naming the same physical column twice raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Duplicate field names could silently replace a mass or velocity array.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 109. `test_npz_rejects_nonmapping_manifest`

[Source](../tests/test_public_input_contracts.py#L206) · **Consolidate without losing conditions** · Role: guard.

**Purpose.** Reject saved data lacking an interpretable catalogue description.
**Pass.** A JSON list in place of the archive manifest raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with schema and duplicate-field tests as malformed-archive cases; keep the load-time boundary.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 110. `test_npz_checks_physical_arrays_on_load`

[Source](../tests/test_public_input_contracts.py#L230) · **Retain** · Role: guard.

**Purpose.** Validate physical values after reading persisted population data.
**Pass.** Replacing a saved physical array entry by NaN causes load to raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Serialization must not bypass the physical-value checks used for direct construction.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 111. `test_native_units_reject_unrepresentable_inputs_and_outputs`

[Source](../tests/test_public_input_contracts.py#L261) · **Retain** · Role: guard.

**Purpose.** Keep canonical unit conversion from discarding complex or nonfinite values.
**Pass.** Both conversion directions reject complex, infinite or NaN inputs.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These checks complement scale-conversion tests and prevent lossy casting.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 112. `test_astropy_units_reject_complex_before_cast`

[Source](../tests/test_public_input_contracts.py#L285) · **Retain** · Role: guard.

**Purpose.** Apply the same real-valued physical-domain rule to quantity-aware conversion.
**Pass.** Complex solar-mass quantities and complex internal values raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain optional-backend parity; silently casting Astropy inputs would bypass native safeguards.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 113. `test_native_cosmology_rejects_invalid_parameters`

[Source](../tests/test_public_input_contracts.py#L320) · **Retain** · Role: guard.

**Purpose.** Restrict the implemented flat-LambdaCDM background to its supported parameter domain.
**Pass.** Matter density outside 0<Omega_m0<=1 and nonpositive/nonfinite h are rejected.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** The restriction describes this implementation, not all possible cosmologies; retain until broader physics is explicitly supported.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 114. `test_native_cosmology_rejects_invalid_redshift`

[Source](../tests/test_public_input_contracts.py#L347) · **Retain** · Role: guard.

**Purpose.** Reject epochs that cannot define the implemented real cosmological background.
**Pass.** All public background functions reject z<=-1, NaN, infinity and complex redshift.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A shared background-domain guard prevents inconsistent treatment of the same epoch across growth, density and time.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 115. `test_native_cosmology_identifies_distinct_float_parameters`

[Source](../tests/test_public_input_contracts.py#L382) · **Retain** · Role: guard.

**Purpose.** Avoid cache aliasing between numerically distinct cosmological parameters.
**Pass.** h=0.674 and h=0.67400000001 produce different identifiers.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Exact input provenance is intentional even for observationally negligible differences; consolidating with general identity tests is possible, but coarsening identities needs a declared tolerance policy.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 116. `test_background_density_overflow_is_explicit`

[Source](../tests/test_public_input_contracts.py#L398) · **Relax only incidental assertions** · Role: guard.

**Purpose.** Do not return an infinite cosmological density as a usable background.
**Pass.** At the stress redshift 1e100, mean and critical density raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Keep overflow detection; move fantastically large epochs to a numerical stress tier if routine scientific CI is simplified.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 117. `test_subnormal_nfw_target_has_a_finite_warning_free_bracket`

[Source](../tests/test_public_input_contracts.py#L416) · **Consolidate without losing conditions** · Role: guard.

**Purpose.** Keep the NFW mass inverse defined at the representable positive-mass boundary.
**Pass.** For a dimensionless target 1e-310, inversion is warning-free and agrees with sqrt(2m) to rtol=1e-10.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Candidate for a separate stress suite or consolidation with small-radius inversion; it does not constrain the resolved subhalo domain by itself.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 118. `test_einstein_de_sitter_limits_and_future_redshift`

[Source](../tests/test_public_input_contracts.py#L438) · **Retain** · Role: guard.

**Purpose.** Recover the analytic matter-only cosmological background.
**Pass.** At z=-0.5,0,1,3, Omega_m=1, D=1/(1+z), and cosmic time follows 2/(3H0)*(1+z)**(-3/2) within the stated roundoff/quadrature tolerances.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain an independent physical limit, including supported future epochs.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_sharp_k_integral.py

#### 119. `test_cutoff_plateau_is_monotonic_without_projection`

[Source](../tests/test_sharp_k_integral.py#L87) · **Retain** · Role: numerical.

**Purpose.** Resolve a suppressed spectrum's variance plateau without imposing monotonicity afterward.
**Pass.** For P(k)=exp(-k), variance decreases with mass and matches the incomplete-gamma integral at rtol=3e-12, atol=3e-16.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This independent integral detects the small-scale cancellation that can corrupt EPS weights.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 120. `test_moving_boundary_derivative_agrees_with_independent_integral`

[Source](../tests/test_sharp_k_integral.py#L111) · **Retain** · Role: numerical.

**Purpose.** Make the sharp-k derivative describe the same continuous variance integral.
**Pass.** A finite difference agrees with the moving-boundary derivative to rtol=4e-8; chunking preserves values exactly; doubling k-grid resolution agrees to 2e-13.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain derivative consistency and convergence; the analytic spectrum makes the tight tolerances meaningful for this controlled case.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 121. `test_tabulated_spectrum_knots_have_exact_piecewise_polynomial_integrals`

[Source](../tests/test_sharp_k_integral.py#L139) · **Retain** · Role: numerical.

**Purpose.** Integrate across input-spectrum knots without missing changes of interpolation slope.
**Pass.** Variance at six off-knot cutoffs matches direct piecewise polynomial integration to rtol=3e-13.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This catches quadrature errors that smooth toy spectra cannot expose; it is a mathematical identity for linear interpolation.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_spatial.py

#### 122. `test_radial_pdf_normalization`

[Source](../tests/test_spatial.py#L15) · **Retain** · Role: numerical.

**Purpose.** Normalize a spatial subhalo distribution with the intended radial measure.
**Pass.** Trapezoidal integration of the normalized toy radial density equals one.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A wrong spatial normalization changes expected subhalo counts in apertures.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 123. `test_shell_probabilities_normalize`

[Source](../tests/test_spatial.py#L33) · **Consolidate without losing conditions** · Role: numerical.

**Purpose.** Keep the total probability of all radial shells equal to unity.
**Pass.** The ten-shell toy measure sums to one.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate with radial-PDF normalization if one parametrized test still covers both point-density and finite-shell representations.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 124. `test_kepler_orbit_turning_points_period_and_measure`

[Source](../tests/test_spatial.py#L52) · **Retain** · Role: numerical.

**Purpose.** Recover an analytic bound orbit and its residence-time distribution.
**Pass.** For GM=a=1 and eccentricity 0.5, turning points are 0.5 and 1.5, period is 2*pi, shell weights are nonnegative and normalized, and missing endpoints are rejected.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This independently validates orbital geometry and time weighting rather than a SASHIMI fit.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 125. `test_turning_points_reject_unbound_scan`

[Source](../tests/test_spatial.py#L146) · **Retain** · Role: numerical.

**Purpose.** Avoid assigning a bound radial period to an orbit without two turning points.
**Pass.** The positive-energy radial Kepler case raises ValueError on the supplied scan interval.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain the bound-orbit domain check; this does not guarantee arbitrary scans find all roots.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_state.py

#### 126. `test_host_and_profile_state_validate_shapes`

[Source](../tests/test_state.py#L7) · **Retain** · Role: guard.

**Purpose.** Keep mass, density and profile fields attached to the same population nodes.
**Pass.** Aligned host arrays and an NFW subhalo state are accepted; a mismatched host-array length raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Misaligned physical state can produce finite but unrelated masses and profiles.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 127. `test_orbital_state_requires_representation`

[Source](../tests/test_state.py#L32) · **Retain** · Role: guard.

**Purpose.** Require an explicit meaning for a stored orbital radius.
**Pass.** A radius-only state without a representation is rejected; an instantaneous-radius state is accepted.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Instantaneous radius and an orbit-averaged measure are not interchangeable observables.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_units.py

#### 128. `test_native_units_reject_nonfinite`

[Source](../tests/test_units.py#L7) · **Consolidate without losing conditions** · Role: numerical.

**Purpose.** Keep canonical units explicit when accepting or returning physical values.
**Pass.** NaN mass, unknown physical dimensions, and an unsupported output-unit string are rejected.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Merge the nonfinite case with public unit-boundary coverage while retaining unknown-unit and dimension cases.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 129. `test_astropy_units_roundtrip`

[Source](../tests/test_units.py#L29) · **Retain** · Role: numerical.

**Purpose.** Convert physical lengths without losing the factor between kpc and Mpc.
**Pass.** 2000 kpc becomes 2 Mpc and converts back; dimensionless arrays are unchanged; seconds cannot be converted to mass.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain a real unit-scale and dimensional-consistency test for the optional backend.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_variance.py

#### 130. `test_callable_variance_model_satisfies_protocol`

[Source](../tests/test_variance.py#L15) · **Retain** · Role: numerical.

**Purpose.** Expose a model-supplied density variance and its mass derivative consistently.
**Pass.** For sigma=M**(-1/4)/(1+z), variance equals sigma squared and its analytic derivative is returned correctly.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain analytic value checks; protocol identity and label spelling alone would be insufficient scientific coverage.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 131. `test_integrated_sharp_k_variance_matches_analytic_power_law`

[Source](../tests/test_variance.py#L49) · **Retain** · Role: numerical.

**Purpose.** Recover the sharp-k integral and growth scaling for a constant spectrum.
**Pass.** At M=1, variance equals (1-1e-9)/(6*pi**2), derivative is -1/(6*pi**2), and both scale as D(z)**2.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** The integral and boundary derivative have exact references; strict algebraic tolerances do not claim equivalent accuracy for an arbitrary physical spectrum.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 132. `test_sharp_k_variance_and_derivative_are_continuous_between_grid_nodes`

[Source](../tests/test_variance.py#L94) · **Retain** · Role: numerical.

**Purpose.** Avoid artificial steps or zero derivatives as the smoothing cutoff moves.
**Pass.** For P(k)=k**(-2), off-grid variance and derivative follow their analytic expressions and the derivative remains negative.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This protects the EPS abundance from artifacts of a discretized wavenumber cutoff.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 133. `test_sharp_k_truncated_domain_has_zero_boundary_derivative`

[Source](../tests/test_variance.py#L138) · **Retain** · Role: numerical.

**Purpose.** Respect the explicitly finite power-spectrum integration domain.
**Pass.** Cutoffs above k_max give a positive constant variance with zero derivative; cutoffs below k_min give zero variance and derivative.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain because a finite input spectrum has different boundary behavior from an extrapolated cosmological spectrum.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 134. `test_integrator_preserves_strict_tabulated_endpoints`

[Source](../tests/test_variance.py#L169) · **Consolidate without losing conditions** · Role: numerical.

**Purpose.** Keep quadrature endpoints inside a non-extrapolating input spectrum.
**Pass.** Integrating a strictly bounded constant-power table returns a finite variance.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Consolidate this smoke case with the analytic sharp-k integral test while preserving strict endpoint enforcement.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 135. `test_variance_cache_is_content_addressed_and_round_trips`

[Source](../tests/test_variance.py#L196) · **Retain** · Role: numerical.

**Purpose.** Tie saved variance data to the filter settings used to compute them.
**Pass.** Changing filter-scale settings changes the key; saved masses and variances reload exactly; a mismatched key is rejected.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Lossless storage and parameter identity are necessary for reproducible halo abundance, independently of quadrature accuracy.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

### tests/test_variance_input_contract.py

#### 136. `test_invalid_coordinates_fail`

[Source](../tests/test_variance_input_contract.py#L114) · **Retain** · Role: guard.

**Purpose.** Reject nonphysical mass and undefined redshift at every variance entry point.
**Pass.** Callable and integrated sigma, variance and derivative reject complex coordinates, nonpositive mass, infinite mass and NaN redshift.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Public adapters must not bypass a shared physical-domain guard.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 137. `test_callable_configuration_is_validated`

[Source](../tests/test_variance_input_contract.py#L147) · **Retain** · Role: guard.

**Purpose.** Require identifiable and callable prescriptions for model-supplied variance.
**Pass.** Blank identity or noncallable sigma/derivative providers raises an error.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** The library cannot interpret an unidentifiable or unevaluable variance model.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 138. `test_callback_results_cannot_corrupt_or_expand_output`

[Source](../tests/test_variance_input_contract.py#L173) · **Retain** · Role: guard.

**Purpose.** Keep returned variance values finite and aligned with the requested halo masses.
**Pass.** NaN, infinity, complex values and an incompatible result length raise ValueError for sigma and derivative callbacks.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** A wrong broadcast axis can alter an entire population without an obvious numerical failure.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 139. `test_negative_sigma_and_square_overflow_fail`

[Source](../tests/test_variance_input_contract.py#L202) · **Retain** · Role: guard.

**Purpose.** Keep the rms density contrast nonnegative and its variance finite.
**Pass.** A negative sigma or sigma whose square overflows raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These are physical and representability constraints, not arbitrary stylistic restrictions.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 140. `test_scalar_callback_can_broadcast_but_invalid_input_never_reaches_callback`

[Source](../tests/test_variance_input_contract.py#L221) · **Retain** · Role: guard.

**Purpose.** Allow a constant physical prescription while rejecting invalid masses early.
**Pass.** A scalar sigma broadcasts over a 2-by-3 mass/redshift grid; negative mass fails before calling the supplied model.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** This positive control ensures validation does not exclude legitimate constant approximations.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 141. `test_complex_power_is_rejected_before_cast`

[Source](../tests/test_variance_input_contract.py#L265) · **Retain** · Role: guard.

**Purpose.** Prevent complex power from being silently used in real density-variance quadrature.
**Pass.** Sharp-k and top-hat paths reject complex spectra for both variance and its derivative.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Retain even if tabulated spectra validate input: custom power callbacks are another entry point.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 142. `test_invalid_growth_and_overflow_fail`

[Source](../tests/test_variance_input_contract.py#L329) · **Retain** · Role: guard.

**Purpose.** Reject growth factors that cannot define real finite variance at each epoch.
**Pass.** Complex, wrong-shape or overflowing growth factors raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** An invalid growth callback otherwise changes redshift-dependent abundance silently.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 143. `test_invalid_window_shape_is_not_silently_broadcast`

[Source](../tests/test_variance_input_contract.py#L353) · **Retain** · Role: guard.

**Purpose.** Keep smoothing-window values on the intended mass-by-wavenumber grid.
**Pass.** A window returning an incompatible shape raises ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** Accidental broadcasting can mix filter radii and produce plausible but wrong integrals.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 144. `test_complex_integration_knots_and_configuration_fail`

[Source](../tests/test_variance_input_contract.py#L406) · **Retain** · Role: guard.

**Purpose.** Require real quadrature breakpoints and mean density for mass smoothing.
**Pass.** Complex integration knots or complex rho_mean raise ValueError.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** These define the integration domain and mass-radius conversion and must not be truncated to real parts.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

#### 145. `test_valid_broadcast_empty_and_derivative_signs_are_preserved`

[Source](../tests/test_variance_input_contract.py#L427) · **Retain** · Role: guard.

**Purpose.** Preserve legitimate variance broadcasting and user-supplied derivative signs.
**Pass.** Analytic sigma and variance broadcast correctly, empty queries stay empty, and a positive supplied derivative is retained.
**Fail.** Any stated comparison/identity mismatch, missing expected exception/warning, or unexpected calculation error fails the case.
**Necessity / proposal.** The generic library must not impose monotonic CDM-like behavior on every supplied physical model.
**Limit.** This controlled calculation does not calibrate a subhalo population against simulations or observations.

## Endpoint-boundary regression follow-up — 2026-10-22

This follow-up reviews two new definitions (ten parametrized cases) against
review head `b7a9fde9d6d8b380c94512337c806851ac5c78f7`. The 145-definition
inventory and validation counts above remain the historical October 14 audit.
No existing test body, assertion, tolerance or historical numerical artifact is
changed or weakened by this follow-up.

### `test_endpoint_respects_bounded_callback_extent`

[Source](../tests/test_picard.py#L371) · **Retain** · Role: numerical-domain guard.

**Purpose.** Permit valid in-domain tidal calculations when all coefficient
callbacks reject extrapolation. The six cases distinguish absent versus explicit
validation domains, a domain smaller than the solver extent, linear versus cubic
tables, and interior/exact-upper/preceding-float accretion queries.
**Pass.** Endpoint, history and direct results reproduce the independent
constant-rate analytic solution at `rtol=5e-13`, without extending the sampled
table or taking a fallback. The six cases fail on the prior unconditional
upper-bound padding.
**Necessity / proposal.** Existing analytic tests have unbounded callbacks and
cannot detect accidental sampling beyond a model's supported redshift range.
Retain these range-rejecting callbacks alongside the analytic controls.
**Limit.** This detects unsupported callback evaluation; it does not establish
quadrature accuracy or physical calibration for varying model coefficients.

### `test_endpoint_integrates_interval_without_distinct_table_nodes`

[Source](../tests/test_picard.py#L430) · **Retain** · Role: representability guard.

**Purpose.** Handle adjacent floating-point epochs when the requested table grid
cannot contain enough distinct redshifts. Four cases cover different floating
spacings and bounded versus extrapolating coefficient callbacks.
**Pass.** A counted and warned direct-ODE fallback preserves the requested
duration, both input masses and agreement with the analytic/history/direct
solutions at `rtol=5e-14`; no degenerate table is cached.
**Necessity / proposal.** Upper-domain padding previously hid this degenerate
grid for permissive callbacks. Removing the padding must preserve their valid
near-zero-duration queries without evaluating bounded callbacks out of range
or silently replacing a nonzero interval with zero evolution.
**Limit.** This protects safe numerical routing and diagnostics, not a
finite-time accuracy or physical-model calibration claim.
