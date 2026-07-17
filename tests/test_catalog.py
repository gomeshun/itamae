import numpy as np
import pytest

from itamae.types import CatalogMetadata, WeightedSubhaloCatalog


def metadata():
    return CatalogMetadata(
        model_identifier="test-model:v1",
        backend_identifier="array=numpy;cosmology=test;units=test",
        source_identifier="fixture:test-catalog",
    )


def catalog():
    return WeightedSubhaloCatalog(
        columns={"mass": np.array([1.0, 2.0, 3.0])},
        weights={
            "weight_base": np.array([1.0, 2.0, 3.0]),
            "weight_survival": np.array([1.0, 0.5, 0.0]),
        },
        metadata=metadata(),
    )


def test_catalog_factorized_weights_and_selection():
    weighted = catalog()
    assert np.allclose(weighted.weight_final, [1.0, 1.0, 0.0])
    assert weighted.weighted_sum(weighted.columns["mass"]) == pytest.approx(3.0)
    assert weighted.select(np.array([True, False, True])).shape == (2,)


def test_catalog_rejects_invalid_schema_and_weights():
    with pytest.raises(ValueError, match="weight_base"):
        WeightedSubhaloCatalog(
            columns={"x": np.ones(2)},
            weights={"weight_survival": np.ones(2)},
            metadata=metadata(),
        )
    with pytest.raises(ValueError, match="same shape"):
        WeightedSubhaloCatalog(
            columns={"x": np.ones(2)},
            weights={"weight_base": np.ones(3)},
            metadata=metadata(),
        )
    with pytest.raises(ValueError, match="negative"):
        WeightedSubhaloCatalog(
            columns={"x": np.ones(2)},
            weights={"weight_base": np.array([1.0, -1.0])},
            metadata=metadata(),
        )
    with pytest.raises(ValueError, match="missing required fields"):
        WeightedSubhaloCatalog(
            columns={"x": np.ones(2)},
            weights={"weight_base": np.ones(2)},
            metadata={},
        )
    with pytest.raises(ValueError, match="reserved fields"):
        CatalogMetadata(
            model_identifier="test",
            backend_identifier="test",
            extra={"schema_version": "override"},
        )


def test_catalog_histogram_realization_and_concatenation():
    weighted = catalog()
    counts, edges = weighted.weighted_histogram("mass", bins=[0.0, 2.0, 4.0])
    np.testing.assert_allclose(counts, [1.0, 1.0])
    np.testing.assert_allclose(edges, [0.0, 2.0, 4.0])

    realized = weighted.poisson_realization(np.random.default_rng(12))
    assert set(realized) == {"mass"}
    assert realized["mass"].ndim == 1

    combined = WeightedSubhaloCatalog.concatenate([weighted, weighted])
    assert combined.shape == (6,)
    np.testing.assert_allclose(combined.weight_final, np.tile(weighted.weight_final, 2))


def test_poisson_realization_requires_explicit_generator():
    with pytest.raises(TypeError, match="Generator"):
        catalog().poisson_realization(None)


def test_catalog_npz_round_trip_without_pickle(tmp_path):
    original = catalog()
    path = tmp_path / "catalog.npz"
    original.to_npz(path)
    restored = WeightedSubhaloCatalog.from_npz(path)

    assert dict(restored.metadata) == dict(original.metadata)
    assert set(restored.columns) == set(original.columns)
    assert set(restored.weights) == set(original.weights)
    for name in original.columns:
        np.testing.assert_array_equal(restored.columns[name], original.columns[name])
    for name in original.weights:
        np.testing.assert_array_equal(restored.weights[name], original.weights[name])


def test_catalog_npz_rejects_wrong_archive_schema(tmp_path):
    path = tmp_path / "invalid.npz"
    np.savez(
        path,
        manifest_json=np.asarray(
            '{"archive_schema":"unknown","catalog_metadata":{},"columns":[],"weights":[]}'
        ),
    )
    with pytest.raises(ValueError, match="schema"):
        WeightedSubhaloCatalog.from_npz(path)
