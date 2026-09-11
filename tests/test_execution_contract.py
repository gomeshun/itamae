"""Failure boundaries and deterministic transport, independent of SASHIMI physics."""

from dataclasses import replace

import numpy as np
import pytest

from itamae.execution import PopulationPipeline, concatenate_accretion_batches
from itamae.measure import build_accretion_batch
from itamae.types import CatalogMetadata


def batch(size=4, **metadata):
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
    calls = []
    model = pipeline(initialize=lambda *args: calls.append(True) or {})
    with pytest.raises(ValueError, match=f"batch 1.*metadata.*{key}"):
        model.execute([batch(**{key: left}), batch(**{key: right})])
    assert calls == []


def test_concatenation_checks_nested_metadata_independent_of_mapping_order():
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
    with pytest.raises((ValueError, TypeError), match="survival"):
        pipeline(survival=lambda *args: mask).execute([batch()])


@pytest.mark.parametrize("mask", [[False, True, False, True], [0, 1, 0, 1]])
def test_boolean_and_binary_masks_retain_selection(mask):
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
    field = {"initial": "initialize", "evolved": "evolve"}.get(stage, stage)
    with pytest.raises(ValueError, match=f"batch 0.*{stage}") as exc:
        pipeline(**{field: callback}).execute([batch()])
    assert exc.value.__cause__ is not None


def test_callback_exception_reports_component_and_retains_cause():
    def broken_evolver(*args):
        raise ArithmeticError("model domain")

    with pytest.raises(ValueError, match="batch 0.*evolved.*broken_evolver") as exc:
        pipeline(evolve=broken_evolver).execute([batch()])
    assert isinstance(exc.value.__cause__, ArithmeticError)


def test_column_order_is_not_a_schema_difference():
    def columns(batch, initial, evolved, masks, context):
        fields = {"m_acc": batch.mvir_acc, "m_bound": evolved["mass"]}
        return fields if context == 0 else dict(reversed(list(fields.items())))

    result = pipeline(columns=columns).execute([batch(), batch()], contexts=[0, 1])
    np.testing.assert_array_equal(result.columns["m_acc"], np.tile(np.arange(4) + 2, 2))


def test_optional_weight_mismatch_precedes_callbacks():
    calls = []
    with pytest.raises(ValueError, match="weight_host_history"):
        pipeline(initialize=lambda *args: calls.append(True) or {}).execute(
            [batch(), replace(batch(), weight_host_history=None)]
        )
    assert calls == []


def test_scalar_batch_has_an_explicit_population_axis_error():
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
