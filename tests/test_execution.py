import numpy as np

from itamae.execution import PopulationPipeline
from itamae.measure import build_accretion_batch
from itamae.protocols import VarianceModel
from itamae.protocols.execution import VarianceModel as ExecutionVarianceModel
from itamae.types import CatalogMetadata


def _batch(offset: float):
    return build_accretion_batch(
        np.array([1.0, 2.0]) + offset,
        0.5 + offset,
        np.array([4.0, 5.0]),
        np.array([0.2, 0.3]),
        np.array([0.7, 0.8]),
        mvir_acc=np.array([1.1, 2.2]) + offset,
        metadata={"batch": offset},
    )


def test_execution_reuses_canonical_variance_protocol():
    assert VarianceModel is ExecutionVarianceModel


def test_population_pipeline_transports_stages_and_factorized_weights():
    pipeline = PopulationPipeline(
        initialize=lambda batch, context: {"initial": batch.mvir_acc},
        evolve=lambda batch, initial, context: {"m_bound": initial["initial"] * context},
        survival=lambda batch, initial, evolved, context: {
            "cdm": evolved["m_bound"] > 1.0,
            "sidm": evolved["m_bound"] > 2.0,
        },
        columns=lambda batch, initial, evolved, survival, context: {
            "m200_acc": batch.m200_acc,
            "m_bound": evolved["m_bound"],
        },
    )
    execution = pipeline.execute([_batch(0.0), _batch(1.0)], contexts=[1.0, 2.0])

    np.testing.assert_allclose(execution.columns["m_bound"], [1.1, 2.2, 4.2, 6.4])
    np.testing.assert_allclose(execution.weight_factors["weight_base"], [0.2, 0.3, 0.2, 0.3])
    cdm = execution.to_catalog(
        CatalogMetadata(model_identifier="toy", backend_identifier="numpy"),
        view="cdm",
    )
    np.testing.assert_allclose(cdm.weight_final, [0.14, 0.24, 0.14, 0.24])


def test_population_pipeline_rejects_misaligned_contexts():
    pipeline = PopulationPipeline(
        initialize=lambda batch, context: {},
        evolve=lambda batch, initial, context: {},
        survival=lambda batch, initial, evolved, context: np.ones(batch.m200_acc.shape, dtype=bool),
        columns=lambda batch, initial, evolved, survival, context: {
            "m200_acc": batch.m200_acc,
        },
    )
    try:
        pipeline.execute([_batch(0.0)], contexts=[])
    except ValueError as exc:
        assert "one entry per" in str(exc)
    else:
        raise AssertionError("Expected a context-length validation error")
