import numpy as np

from itamae.execution import PopulationPipeline
from itamae.measure import build_accretion_batch
from itamae.protocols import VarianceModel
from itamae.protocols.execution import VarianceModel as ExecutionVarianceModel
from itamae.types import CatalogMetadata


def _batch(offset: float):
    """Create two distinguishable accretion nodes with known concentration and abundance
    factors.

    Parameters
    ----------
    offset : float
        Controlled shift separating the toy population batches while preserving their known
        weight factors.

    Returns
    -------
    AccretionBatch
        Two aligned toy nodes with distinct M200 and virial masses.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    return build_accretion_batch(
        np.array([1.0, 2.0]) + offset,
        0.5 + offset,
        np.array([4.0, 5.0]),
        np.array([0.2, 0.3]),
        np.array([0.7, 0.8]),
        mvir_acc=np.array([1.1, 2.2]) + offset,
        metadata={"model": "toy"},
    )


def test_execution_reuses_canonical_variance_protocol():
    """Track the compatibility alias for the variance interface used by population execution.

    Notes
    -----
    Pass criteria: The re-exported VarianceModel is the same Python object as the canonical
    protocol.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This is an API or output-compatibility check; it does not test a physical
    equation or numerical accuracy.
    """
    assert VarianceModel is ExecutionVarianceModel


def test_population_pipeline_transports_stages_and_factorized_weights():
    """Keep each accretion node aligned with its evolved mass and abundance factors.

    Notes
    -----
    Pass criteria: The two toy slices yield bound masses [1.1,2.2,4.2,6.4], unchanged base
    factors, and final CDM weights [0.14,0.24,0.14,0.24].

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
    """Prevent redshift or host context from being assigned to the wrong accretion slice.

    Notes
    -----
    Pass criteria: A context list with a different length from the batch list raises
    ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
