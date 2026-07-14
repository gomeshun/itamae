import numpy as np

from itamae.measure import build_accretion_batch


def test_build_accretion_batch_broadcasts_and_preserves_weights():
    batch = build_accretion_batch(
        np.array([1.0, 2.0]),
        0.5,
        10.0,
        np.array([0.2, 0.3]),
        1.0,
        mvir_acc=np.array([1.1, 2.2]),
        weight_host_history=0.5,
        metadata={"model_identifier": "test"},
    )
    assert batch.m200_acc.shape == (2,)
    np.testing.assert_allclose(batch.mvir_acc, [1.1, 2.2])
    np.testing.assert_allclose(batch.weight_concentration, 1.0)
    np.testing.assert_allclose(batch.weight_host_history, 0.5)
