"""Failure boundaries and deterministic transport, independent of SASHIMI physics."""

from dataclasses import replace

import numpy as np
import pytest

from itamae.execution import PopulationPipeline, concatenate_accretion_batches
from itamae.measure import build_accretion_batch
from itamae.types import CatalogMetadata


def batch(size=4, **metadata):
    """Create an accretion batch with explicit independent host and concentration factors.

    Parameters
    ----------
    size : int
        Number of grid or population nodes requested.
    **metadata : CatalogMetadata or Mapping
        Physical model, backend, unit and source identities needed to interpret or combine
        the nodes.

    Returns
    -------
    AccretionBatch
        size aligned nodes with caller-specified metadata.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    return build_accretion_batch(
        np.arange(size, dtype=float) + 1,
        0.5,
        8.0,
        0.25,
        0.5,
        mvir_acc=np.arange(size, dtype=float) + 2,
        weight_host_history=0.8,
        metadata=metadata,
    )


def pipeline(**overrides):
    """Compose deterministic initialization, half-mass evolution and a survival threshold.

    Parameters
    ----------
    **overrides : dict
        Deliberate changes to the otherwise valid controlled model, stage or numerical
        settings.

    Returns
    -------
    PopulationPipeline
        Toy pipeline with optional stage replacements for invalid-contract checks.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    stages = dict(
        initialize=lambda batch, context: {"mass": batch.mvir_acc},
        evolve=lambda batch, initial, context: {"mass": initial["mass"] * 0.5},
        survival=lambda batch, initial, evolved, context: evolved["mass"] > 1.5,
        columns=lambda batch, initial, evolved, masks, context: {
            "m_acc": batch.mvir_acc,
            "m_bound": evolved["mass"],
        },
    )
    return PopulationPipeline(**(stages | overrides))


@pytest.mark.parametrize(
    "key,left,right",
    [
        ("backend", "native", "colossus"),
        ("units", {"mass": "Msun"}, {"mass": "Msun/h"}),
        ("source_revision", "a" * 40, "b" * 40),
        ("schema", 1, 2),
        ("physics_mode", "legacy", "consistent"),
    ],
)
def test_mixed_metadata_rejected_before_callbacks(key, left, right):
    """Prevent merging populations with different cosmologies, units, revisions, or
    prescriptions.

    Parameters
    ----------
    key : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.
    left : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.
    right : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Each incompatible pair raises ValueError before any initialization
    callback runs.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    calls = []
    model = pipeline(initialize=lambda *args: calls.append(True) or {})
    with pytest.raises(ValueError, match=f"batch 1.*metadata.*{key}"):
        model.execute([batch(**{key: left}), batch(**{key: right})])
    assert calls == []


def test_concatenation_checks_nested_metadata_independent_of_mapping_order():
    """Distinguish physical metadata content from dictionary insertion order.

    Notes
    -----
    Pass criteria: Equivalent nested cosmology metadata concatenate to eight nodes; changed
    or absent cosmology data are rejected.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    first = batch(cosmology={"grid": np.array([0.0, 1.0]), "h": 0.7})
    second = batch(cosmology={"h": 0.7, "grid": np.array([0.0, 1.0])})
    result = concatenate_accretion_batches([first, second])
    assert result.m200_acc.size == 8
    with pytest.raises(ValueError, match="metadata.*cosmology.*grid"):
        concatenate_accretion_batches(
            [first, batch(cosmology={"h": 0.7, "grid": np.array([0.0, 2.0])})]
        )
    with pytest.raises(ValueError, match="metadata"):
        concatenate_accretion_batches([first, batch()])


@pytest.mark.parametrize("mask", [[0.0, 0.5, 1.0, 1.0], [-1, 0, 1, 0], ["yes"] * 4])
def test_survival_rejects_truthy_non_masks(mask):
    """Keep a disruption decision distinct from a fractional population weight.

    Parameters
    ----------
    mask : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Fractional, negative, or text-valued survival masks raise an error.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises((ValueError, TypeError), match="survival"):
        pipeline(survival=lambda *args: mask).execute([batch()])


@pytest.mark.parametrize("mask", [[False, True, False, True], [0, 1, 0, 1]])
def test_boolean_and_binary_masks_retain_selection(mask):
    """Preserve the same disrupted and surviving nodes for supported mask encodings.

    Parameters
    ----------
    mask : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Boolean and exact 0/1 masks are transported without changing selected
    nodes.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    execution = pipeline(survival=lambda *args: mask).execute([batch()])
    np.testing.assert_array_equal(execution.survival["default"], mask)


@pytest.mark.parametrize(
    "stage,callback",
    [
        ("initial", lambda *args: {"bad": np.ones(3)}),
        ("evolved", lambda *args: {"bad": np.full(4, np.nan)}),
        ("survival", lambda *args: {}),
        ("columns", lambda *args: {}),
        ("columns", lambda *args: {"bad": np.full(4, "x")}),
    ],
)
def test_errors_identify_batch_stage_and_original_cause(stage, callback):
    """Make a failed structure or survival stage traceable to its physical slice.

    Parameters
    ----------
    stage : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.
    callback : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Malformed stage outputs raise ValueError naming the batch and stage while
    retaining the original exception.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    field = {"initial": "initialize", "evolved": "evolve"}.get(stage, stage)
    with pytest.raises(ValueError, match=f"batch 0.*{stage}") as exc:
        pipeline(**{field: callback}).execute([batch()])
    assert exc.value.__cause__ is not None


def test_callback_exception_reports_component_and_retains_cause():
    """Preserve a model-domain failure rather than returning a partial population.

    Notes
    -----
    Pass criteria: An ArithmeticError from the evolver becomes a contextual error and
    remains the chained cause.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    def broken_evolver(*args):
        """Raise a controlled physical-domain error to verify that execution retains its cause.

        Parameters
        ----------
        *args : tuple
            Additional positional parameters forwarded unchanged to the supplied numerical or
            physical callback.

        Notes
        -----
        Preserve a model-domain failure rather than returning a partial population. Pass
        criteria: An ArithmeticError from the evolver becomes a contextual error and remains the
        chained cause. Failure occurs when the stated comparison or expected rejection is not
        met. This controlled calculation does not calibrate a subhalo population against
        simulations or observations.
        """
        raise ArithmeticError("model domain")

    with pytest.raises(ValueError, match="batch 0.*evolved.*broken_evolver") as exc:
        pipeline(evolve=broken_evolver).execute([batch()])
    assert isinstance(exc.value.__cause__, ArithmeticError)


def test_column_order_is_not_a_schema_difference():
    """Join physically identical catalogues despite different field insertion order.

    Notes
    -----
    Pass criteria: Two toy slices with reversed dictionary order retain the expected
    accreted-mass sequence.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """

    def columns(batch, initial, evolved, masks, context):
        """Return identical physical fields in different mapping orders.

        Parameters
        ----------
        batch : AccretionBatch
            Aligned accretion nodes with separate population, concentration and optional
            host/orbit factors.
        initial : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.
        evolved : bool
            Use present bound masses/profiles when true or accretion masses/profiles when false,
            under the observable population selection.
        masks : Mapping
            Named survival masks from the controlled population evolution.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        dict
            Accretion and bound-mass arrays; ordering varies without changing the scientific
            schema.

        Notes
        -----
        Join physically identical catalogues despite different field insertion order. Pass
        criteria: Two toy slices with reversed dictionary order retain the expected
        accreted-mass sequence. Failure occurs when the stated comparison or expected rejection
        is not met. This controlled calculation does not calibrate a subhalo population against
        simulations or observations.
        """
        fields = {"m_acc": batch.mvir_acc, "m_bound": evolved["mass"]}
        return fields if context == 0 else dict(reversed(list(fields.items())))

    result = pipeline(columns=columns).execute([batch(), batch()], contexts=[0, 1])
    np.testing.assert_array_equal(result.columns["m_acc"], np.tile(np.arange(4) + 2, 2))


def test_optional_weight_mismatch_precedes_callbacks():
    """Reject a population mix where only some slices include host-history weighting.

    Notes
    -----
    Pass criteria: Presence mismatch of weight_host_history raises ValueError before
    callbacks run.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    calls = []
    with pytest.raises(ValueError, match="weight_host_history"):
        pipeline(initialize=lambda *args: calls.append(True) or {}).execute(
            [batch(), replace(batch(), weight_host_history=None)]
        )
    assert calls == []


def test_scalar_batch_has_an_explicit_population_axis_error():
    """Reject a batch that lacks a node axis needed for population alignment.

    Notes
    -----
    Pass criteria: Replacing every one-node array by a scalar raises a population-axis
    error.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    scalar = replace(
        batch(1),
        **{
            name: np.asarray(getattr(batch(1), name)[0])
            for name in (
                "m200_acc",
                "mvir_acc",
                "z_acc",
                "concentration_acc",
                "weight_base",
                "weight_concentration",
                "weight_host_history",
            )
        },
    )
    with pytest.raises(ValueError, match="population axis"):
        pipeline().execute([scalar])


def test_empty_batches_zero_survivors_and_deterministic_partition():
    """Keep population results independent of how accretion nodes are partitioned.

    Notes
    -----
    Pass criteria: Empty slices and zero survivors are handled, no-batch input fails, and
    partitioned versus whole populations agree exactly including a fixed-seed realization.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    with pytest.raises(ValueError, match="At least one"):
        pipeline().execute([])
    empty = pipeline().execute([batch(0)])
    assert empty.columns["m_acc"].size == 0
    source = batch(7, specification="toy:1")
    chunks = []
    for section in (slice(0, 2), slice(2, 2), slice(2, 5), slice(5, 7)):
        chunks.append(
            replace(
                source,
                **{
                    name: getattr(source, name)[section]
                    for name in (
                        "m200_acc",
                        "mvir_acc",
                        "z_acc",
                        "concentration_acc",
                        "weight_base",
                        "weight_concentration",
                        "weight_host_history",
                    )
                },
            )
        )
    complete = pipeline().execute([source])
    partitioned = pipeline().execute(chunks)
    for name in complete.columns:
        np.testing.assert_array_equal(partitioned.columns[name], complete.columns[name])
    metadata = CatalogMetadata(model_identifier="toy", backend_identifier="numpy")
    left, right = complete.to_catalog(metadata), partitioned.to_catalog(metadata)
    np.testing.assert_array_equal(left.weight_final, right.weight_final)
    np.testing.assert_array_equal(
        left.poisson_realization(np.random.default_rng(17))["m_bound"],
        right.poisson_realization(np.random.default_rng(17))["m_bound"],
    )
    zero = pipeline(survival=lambda batch, *args: np.zeros(batch.m200_acc.shape, bool))
    assert np.all(zero.execute([source]).to_catalog(metadata).weight_final == 0)
