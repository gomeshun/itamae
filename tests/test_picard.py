"""Independent analytic, domain, history and cache contracts for shared Picard."""

import json

import numpy as np
import pytest

from itamae.evolution import PicardDomain, PicardFallbackWarning, PicardMassLossSolver


def make_solver(rate=0.3, exponent=0.04, domain=None):
    state = {"rate": rate, "host": 1e12}
    solver = PicardMassLossSolver(
        lambda z: np.full_like(np.asarray(z, dtype=float), state["rate"]),
        lambda z: np.full_like(np.asarray(z, dtype=float), exponent),
        lambda z: np.full_like(np.asarray(z, dtype=float), state["host"]),
        z_min=0.0,
        z_max=7.0,
        state_key=lambda: tuple(state.items()),
        domain=domain,
    )
    return solver, state


def exact(mass, za, z, rate=0.3, exponent=0.04, host=1e12):
    if exponent == 0:
        return mass * np.exp(rate * (np.asarray(z) - za))
    return host * ((mass / host) ** -exponent - exponent * rate * (np.asarray(z) - za)) ** (
        -1 / exponent
    )


@pytest.mark.parametrize("exponent", [0.0, 0.04, -0.01])
def test_endpoint_and_history_match_independent_analytic_solution(exponent):
    solver, _ = make_solver(exponent=exponent)
    masses = np.geomspace(1e2, 1e10, 17)
    z = np.array([5.0, 4.3, 2.0, 2.0, 0.0])
    expected = exact(masses[None, :], 5.0, z[:, None], exponent=exponent)
    np.testing.assert_allclose(solver.history(masses, 5.0, z), expected, rtol=5e-4)
    np.testing.assert_allclose(solver.endpoint(masses, 5.0, 0.0), expected[-1], rtol=5e-4)
    np.testing.assert_allclose(solver.direct(masses, 5.0, z), expected, rtol=5e-10)
    # Integrate from accretion even when its time is omitted from the output.
    np.testing.assert_allclose(solver.direct(masses, 5.0, z[1:]), expected[1:], rtol=5e-10)


def test_grid_refinement_reduces_time_dependent_coefficient_error():
    solver, _ = make_solver(exponent=0.0)
    solver._rate = lambda z: 0.1 * np.asarray(z) ** 2
    mass, za = 1e6, 5.0
    expected = mass * np.exp(-0.1 * za**3 / 3)
    coarse = solver.history(mass, za, 0.0, n_integration=33)
    fine = solver.history(mass, za, 0.0, n_integration=129)
    assert abs(fine / expected - 1) < abs(coarse / expected - 1) / 10
    assert abs(fine / expected - 1) < 2e-4


def test_endpoint_broadcasting_and_history_shape():
    solver, _ = make_solver()
    mass = np.array([[1e6], [1e8]])
    za = np.array([0.5, 2.0, 5.0])
    assert solver.endpoint(mass, za, 0.0).shape == (2, 3)
    assert solver.history(mass, 5.0, [5.0, 0.0]).shape == (2, 2, 1)
    assert solver.history(mass, 5.0, 0.0).shape == (2, 1)
    np.testing.assert_array_equal(solver.history(mass, 2.0, [2.0, 2.0]), np.stack([mass, mass]))
    assert solver.history(mass, 2.0, []).shape == (0, 2, 1)
    assert solver.endpoint([], 2.0, 0.0).shape == (0,)


def test_cache_reuse_and_stale_table_rejection_after_provider_change():
    solver, state = make_solver()
    old = solver.endpoint(1e6, 3.0, 0.0)
    table = next(iter(solver._picard_tables.values()))
    solver.endpoint(1e7, 3.0, 0.0)
    assert next(iter(solver._picard_tables.values())) is table
    state["rate"] = 0.6
    with pytest.raises(ValueError, match="stale"):
        table.mass(1e6, 3.0)
    new = solver.endpoint(1e6, 3.0, 0.0)
    assert new < old
    assert next(iter(solver._picard_tables.values())) is not table


@pytest.mark.parametrize("which", ["endpoint", "history"])
def test_domain_fallback_is_explicit_and_recorded(which):
    solver, _ = make_solver(domain=PicardDomain((0.0, 7.0), (-12.0, 0.0)))
    with pytest.warns(PicardFallbackWarning, match="outside"):
        actual = getattr(solver, which)([1e6, 1e7], 8.0, 0.0)
    np.testing.assert_allclose(actual, exact(np.array([1e6, 1e7]), 8.0, 0.0), rtol=5e-10)
    report = solver.diagnostics()
    assert report["fallback_calls"] == 1
    assert report["fallback_mass_count"] == 2
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("which", ["endpoint", "history"])
def test_unconverged_iteration_uses_direct_reference(which):
    solver, _ = make_solver(rate=2.0, exponent=0.1)
    with pytest.warns(PicardFallbackWarning, match="convergence"):
        actual = getattr(solver, which)(1e10, 7.0, 0.0, n_iterations=1)
    np.testing.assert_allclose(actual, exact(1e10, 7.0, 0.0, rate=2.0, exponent=0.1), rtol=5e-10)
    assert solver.diagnostics()["fallback_mass_count"] == 1


@pytest.mark.parametrize(
    "option",
    [
        {"n_iterations": 0},
        {"n_z_acc": 2.5},
        {"n_log_ratio": True},
        {"convergence_tolerance": 0},
        {"interpolation": "unknown"},
    ],
)
def test_invalid_endpoint_options_fail(option):
    solver, _ = make_solver()
    with pytest.raises(ValueError):
        solver.endpoint(1e6, 2.0, 0.0, **option)


def test_bad_inputs_fail_without_clipping_or_extrapolation():
    solver, _ = make_solver()
    for mass, za, z in [(0.0, 2.0, 0.0), (1.0, 2.0, 3.0), (1.0, 2.0, np.nan)]:
        with pytest.raises(ValueError):
            solver.endpoint(mass, za, z)
    with pytest.raises(ValueError, match="non-increasing"):
        solver.history(1e6, 2.0, [0.0, 1.0])
    table = solver.table(z_final=0.0, z_acc_max=2.0, log10_ratio_min=-12.0, log10_ratio_max=0.0)
    with pytest.raises(ValueError, match="range"):
        table.mass(1e6, 3.0)
    with pytest.raises(ValueError, match="range"):
        table.mass(1e14, 1.0)


@pytest.mark.parametrize("which", ["endpoint", "history", "direct"])
def test_invalid_physical_callback_fails_in_every_path(which):
    solver, state = make_solver()
    state["rate"] = -1.0
    with pytest.raises(ValueError, match="nonnegative"):
        getattr(solver, which)(1e6, 2.0, 0.0)
    state["rate"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        getattr(solver, which)(1e6, 2.0, 0.0)


def test_mass_unit_rescaling_preserves_fractional_evolution():
    reference, _ = make_solver()
    changed, state = make_solver()
    state["host"] *= 1e-30
    mass = np.geomspace(1e3, 1e10, 8)
    np.testing.assert_allclose(
        changed.endpoint(mass * 1e-30, 3.0, 0.0) / 1e-30,
        reference.endpoint(mass, 3.0, 0.0),
        rtol=5e-14,
    )
