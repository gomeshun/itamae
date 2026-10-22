from dataclasses import dataclass

import numpy as np

from itamae.execution import PopulationComponents, PopulationPipeline
from itamae.measure import build_accretion_batch
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


@dataclass(frozen=True)
class ToyInitializer:
    """Assign a controlled initial virial-mass offset for stage-order checks.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """

    def initialize(self, batch, context):
        """Add the context offset to each virial accretion mass.

        Parameters
        ----------
        batch : AccretionBatch
            Aligned accretion nodes with separate population, concentration and optional
            host/orbit factors.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        dict
            Initial toy masses, aligned with the batch.

        Notes
        -----
        This fixture/helper supports the controlled checks in this test module; its synthetic
        values are not a calibrated halo model.
        """
        return {"initial": batch.mvir_acc + context["offset"]}


@dataclass(frozen=True)
class ToyEvolver:
    """Apply a known multiplicative mass change for population-transport checks.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """

    def evolve(self, batch, initial, context):
        """Scale the initialized toy masses using the batch context.

        Parameters
        ----------
        batch : AccretionBatch
            Aligned accretion nodes with separate population, concentration and optional
            host/orbit factors.
        initial : Mapping or array_like
            Per-node physical state for the indicated evolution stage; its fields and units are
            owned by the supplied physical model.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        dict
            Known bound-mass array without modifying the statistical weights.

        Notes
        -----
        This fixture/helper supports the controlled checks in this test module; its synthetic
        values are not a calibrated halo model.
        """
        return {"m_bound": initial["initial"] * context["scale"]}


@dataclass(frozen=True)
class ToySurvival:
    """Define two threshold survival views on one transported toy population.

    Parameters
    ----------
    threshold : float
        Toy survival threshold in the fixture mass convention; used to distinguish two masks
        on the same nodes.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """

    threshold: float

    def select(self, batch, initial, evolved, context):
        """Select ordinary and stricter survival masks from the same bound masses.

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
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        dict
            Boolean masks at threshold and threshold+1.

        Notes
        -----
        This fixture/helper supports the controlled checks in this test module; its synthetic
        values are not a calibrated halo model.
        """
        return {
            "default": evolved["m_bound"] > self.threshold,
            "strict": evolved["m_bound"] > self.threshold + 1.0,
        }


@dataclass(frozen=True)
class ToyColumns:
    """Expose initial and evolved toy masses for independent stage comparisons.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """

    def build(self, batch, initial, evolved, survival, context):
        """Retain accretion, initialized and bound masses as separate physical columns.

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
        survival : callable or Mapping
            Model survival selector, or the named boolean masks it produced; survival enters the
            final weight once.
        context : object
            Variant-owned state for this batch, for example host history, epoch or numerical
            settings; it is not interpreted by the generic executor.

        Returns
        -------
        dict
            Three aligned arrays used to verify stage transport.

        Notes
        -----
        This fixture/helper supports the controlled checks in this test module; its synthetic
        values are not a calibrated halo model.
        """
        return {
            "m200_acc": batch.m200_acc,
            "initial": initial["initial"],
            "m_bound": evolved["m_bound"],
        }


def _callback_pipeline():
    """Construct the same toy physical stages as the component-based execution fixture.

    Returns
    -------
    PopulationPipeline
        Callback executor for a controlled equivalence comparison.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    return PopulationPipeline(
        initialize=lambda batch, context: {"initial": batch.mvir_acc + context["offset"]},
        evolve=lambda batch, initial, context: {"m_bound": initial["initial"] * context["scale"]},
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
    """Check the concrete pipeline class returned by component composition.

    Notes
    -----
    Pass criteria: The composed object is an instance of PopulationPipeline.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This is an API or output-compatibility check; it does not test a physical
    equation or numerical accuracy.
    """
    components = PopulationComponents(
        initializer=ToyInitializer(),
        evolver=ToyEvolver(),
        survival=ToySurvival(threshold=2.0),
        columns=ToyColumns(),
    )

    assert isinstance(components.to_pipeline(), PopulationPipeline)


def test_population_components_match_callback_execution_and_catalogs():
    """Preserve the physical population when the same stages are supplied as objects or
    callbacks.

    Notes
    -----
    Pass criteria: Both routes return exactly the same toy columns, survival views, weight
    factors, and strict-view catalogue.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
            component_execution.weight_factors[name],
            callback_execution.weight_factors[name],
        )

    metadata = CatalogMetadata(model_identifier="toy", backend_identifier="numpy")
    callback_catalog = callback_execution.to_catalog(metadata, view="strict")
    component_catalog = component_execution.to_catalog(metadata, view="strict")
    np.testing.assert_array_equal(component_catalog.weight_final, callback_catalog.weight_final)
    for name in callback_catalog.columns:
        np.testing.assert_array_equal(
            component_catalog.columns[name], callback_catalog.columns[name]
        )
