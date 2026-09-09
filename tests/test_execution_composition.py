from dataclasses import dataclass

import numpy as np

from itamae.execution import PopulationComponents, PopulationPipeline
from itamae.measure import build_accretion_batch
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


@dataclass(frozen=True)
class ToyInitializer:
    def initialize(self, batch, context):
        return {"initial": batch.mvir_acc + context["offset"]}


@dataclass(frozen=True)
class ToyEvolver:
    def evolve(self, batch, initial, context):
        return {"m_bound": initial["initial"] * context["scale"]}


@dataclass(frozen=True)
class ToySurvival:
    threshold: float

    def select(self, batch, initial, evolved, context):
        return {
            "default": evolved["m_bound"] > self.threshold,
            "strict": evolved["m_bound"] > self.threshold + 1.0,
        }


@dataclass(frozen=True)
class ToyColumns:
    def build(self, batch, initial, evolved, survival, context):
        return {
            "m200_acc": batch.m200_acc,
            "initial": initial["initial"],
            "m_bound": evolved["m_bound"],
        }


def _callback_pipeline():
    return PopulationPipeline(
        initialize=lambda batch, context: {
            "initial": batch.mvir_acc + context["offset"]
        },
        evolve=lambda batch, initial, context: {
            "m_bound": initial["initial"] * context["scale"]
        },
        survival=lambda batch, initial, evolved, context: {
            "default": evolved["m_bound"] > 2.0,
            "strict": evolved["m_bound"] > 3.0,
        },
        columns=lambda batch, initial, evolved, survival, context: {
            "m200_acc": batch.m200_acc,
            "initial": initial["initial"],
            "m_bound": evolved["m_bound"],
        },
    )


def test_population_components_delegate_to_canonical_pipeline():
    components = PopulationComponents(
        initializer=ToyInitializer(),
        evolver=ToyEvolver(),
        survival=ToySurvival(threshold=2.0),
        columns=ToyColumns(),
    )

    assert isinstance(components.to_pipeline(), PopulationPipeline)


def test_population_components_match_callback_execution_and_catalogs():
    batches = [_batch(0.0), _batch(1.0)]
    contexts = [
        {"offset": 0.1, "scale": 1.0},
        {"offset": 0.2, "scale": 2.0},
    ]
    components = PopulationComponents(
        initializer=ToyInitializer(),
        evolver=ToyEvolver(),
        survival=ToySurvival(threshold=2.0),
        columns=ToyColumns(),
    )

    callback_execution = _callback_pipeline().execute(batches, contexts=contexts)
    component_execution = components.execute(batches, contexts=contexts)

    assert tuple(component_execution.columns) == tuple(callback_execution.columns)
    assert tuple(component_execution.survival) == tuple(callback_execution.survival)
    for name in callback_execution.columns:
        np.testing.assert_array_equal(
            component_execution.columns[name], callback_execution.columns[name]
        )
    for name in callback_execution.survival:
        np.testing.assert_array_equal(
            component_execution.survival[name], callback_execution.survival[name]
        )
    for name in callback_execution.weight_factors:
        np.testing.assert_array_equal(
            component_execution.weight_factors[name], callback_execution.weight_factors[name]
        )

    metadata = CatalogMetadata(model_identifier="toy", backend_identifier="numpy")
    callback_catalog = callback_execution.to_catalog(metadata, view="strict")
    component_catalog = component_execution.to_catalog(metadata, view="strict")
    np.testing.assert_array_equal(
        component_catalog.weight_final, callback_catalog.weight_final
    )
    for name in callback_catalog.columns:
        np.testing.assert_array_equal(
            component_catalog.columns[name], callback_catalog.columns[name]
        )
