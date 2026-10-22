"""Independent analytic, domain, history and cache contracts for shared Picard."""

import json

import numpy as np
import pytest
from scipy.interpolate import interp1d

from itamae.evolution import PicardDomain, PicardFallbackWarning, PicardMassLossSolver


def make_solver(rate=0.3, exponent=0.04, domain=None):
    """Bind a constant-coefficient tidal law to mutable state for cache and analytic-solution
    checks.

    Parameters
    ----------
    rate : float
        Constant nonnegative tidal coefficient Phi in the analytic test law.
    exponent : float
        Constant dimensionless mass-ratio exponent zeta in the analytic test law.
    domain : PicardDomain, optional
        Caller-validated redshift and initial mass-ratio bounds. None declares no such
        restriction; it does not establish validation.

    Returns
    -------
    tuple
        Picard solver and the mutable rate/host dictionary used by its physical-state key.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
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
    """Evaluate the closed-form mass-loss solution for constant host, rate and exponent.

    Parameters
    ----------
    mass : float or array_like
        Positive halo mass. Use solar masses with the canonical backend; a supplied variance
        prescription must use the same mass convention and factors of h.
    za : float or array_like
        Dimensionless accretion redshift. History integration requires one scalar accretion
        epoch.
    z : float or array_like
        Dimensionless redshift. The cosmological background requires z > -1; population
        nodes use nonnegative redshifts.
    rate : float
        Constant nonnegative tidal coefficient Phi in the analytic test law.
    exponent : float
        Constant dimensionless mass-ratio exponent zeta in the analytic test law.
    host : float
        Constant positive host mass in the same unit as mass.

    Returns
    -------
    float or numpy.ndarray
        Analytic evolved mass, using the exponential limit when zeta=0.

    Notes
    -----
    This fixture/helper supports the controlled checks in this test module; its synthetic
    values are not a calibrated halo model.
    """
    if exponent == 0:
        return mass * np.exp(rate * (np.asarray(z) - za))
    return host * ((mass / host) ** -exponent - exponent * rate * (np.asarray(z) - za)) ** (
        -1 / exponent
    )


@pytest.mark.parametrize("exponent", [0.0, 0.04, -0.01])
def test_endpoint_and_history_match_independent_analytic_solution(exponent):
    """Recover the analytic tidal mass history for constant host and loss coefficients.

    Parameters
    ----------
    exponent : float
        Dimensionless exponent of the subhalo-to-host mass ratio in the analytic tidal-loss
        law.

    Notes
    -----
    Pass criteria: For exponents 0,0.04,-0.01, Picard masses agree within rtol=5e-4 and
    DOP853 within 5e-10, including repeated outputs and an omitted accretion-time output.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
    """Separate redshift-quadrature error from nonlinear Picard convergence.

    Notes
    -----
    Pass criteria: For rate=0.1*z**2 and zero exponent, refinement from 33 to 129 nodes
    reduces analytic endpoint error by more than tenfold and below 2e-4.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    solver, _ = make_solver(exponent=0.0)
    solver._rate = lambda z: 0.1 * np.asarray(z) ** 2
    mass, za = 1e6, 5.0
    expected = mass * np.exp(-0.1 * za**3 / 3)
    coarse = solver.history(mass, za, 0.0, n_integration=33)
    fine = solver.history(mass, za, 0.0, n_integration=129)
    assert abs(fine / expected - 1) < abs(coarse / expected - 1) / 10
    assert abs(fine / expected - 1) < 2e-4


def test_endpoint_broadcasting_and_history_shape():
    """Keep accretion masses and requested epochs aligned in vectorized stripping queries.

    Notes
    -----
    Pass criteria: Endpoint broadcasting and time-leading history shapes are correct;
    zero-duration histories preserve mass; empty inputs preserve empty axes.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
    """Recompute stripping when the physical mass-loss rate changes.

    Notes
    -----
    Pass criteria: Unchanged state reuses a table; doubling the rate rejects the old table,
    builds a replacement and lowers final mass.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
    """Use a direct physical integration outside the caller-validated Picard domain.

    Parameters
    ----------
    which : str
        Endpoint, full-history or direct tidal integration path under test.

    Notes
    -----
    Pass criteria: A z_acc=8 query outside the declared domain warns, matches the analytic
    solution to 5e-10, and records one fallback covering two masses in JSON-safe
    diagnostics.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
    """Avoid accepting an unconverged nonlinear mass-loss approximation.

    Parameters
    ----------
    which : str
        Endpoint, full-history or direct tidal integration path under test.

    Notes
    -----
    Pass criteria: A deliberately difficult one-iteration case warns, returns the analytic
    reference to 5e-10, and counts the fallback mass.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
    """Reject resolution or interpolation settings that cannot define a valid stripping table.

    Parameters
    ----------
    option : object
        Controlled input supplied by the parametrization or fixture for the specific
        physical/reference case described below.

    Notes
    -----
    Pass criteria: Zero iterations, fractional node counts, Boolean counts, zero tolerance
    and unknown interpolation each raise ValueError.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    solver, _ = make_solver()
    with pytest.raises(ValueError):
        solver.endpoint(1e6, 2.0, 0.0, **option)


def test_bad_inputs_fail_without_clipping_or_extrapolation():
    """Keep invalid masses and reversed cosmic evolution out of the stripping result.

    Notes
    -----
    Pass criteria: Nonpositive mass, final redshift above accretion, NaN epochs, increasing
    history grids, and direct table queries beyond bounds are rejected.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
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
    """Apply the same physical coefficient checks to Picard and direct integration.

    Parameters
    ----------
    which : str
        Endpoint, full-history or direct tidal integration path under test.

    Notes
    -----
    Pass criteria: Negative or NaN stripping rates raise ValueError for endpoint, history
    and direct paths.

    Failure criteria: a stated numerical/identity mismatch, a missing required exception or
    warning, or an unexpected calculation failure fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    solver, state = make_solver()
    state["rate"] = -1.0
    with pytest.raises(ValueError, match="nonnegative"):
        getattr(solver, which)(1e6, 2.0, 0.0)
    state["rate"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        getattr(solver, which)(1e6, 2.0, 0.0)


def test_mass_unit_rescaling_preserves_fractional_evolution():
    """Make stripping depend on the subhalo-to-host mass ratio rather than a chosen mass unit.

    Notes
    -----
    Pass criteria: Scaling host and subhalo masses together by 1e-30 preserves final masses
    after rescaling to rtol=5e-14.

    Failure criteria: any mismatch with the stated value, shape or identity requirements
    fails the test.

    Scope: This controlled calculation does not calibrate a subhalo population against
    simulations or observations.
    """
    reference, _ = make_solver()
    changed, state = make_solver()
    state["host"] *= 1e-30
    mass = np.geomspace(1e3, 1e10, 8)
    np.testing.assert_allclose(
        changed.endpoint(mass * 1e-30, 3.0, 0.0) / 1e-30,
        reference.endpoint(mass, 3.0, 0.0),
        rtol=5e-14,
    )


@pytest.mark.parametrize("declared_upper", [None, 7.0, 4.25])
@pytest.mark.parametrize("n_z_acc", [2, 48])
def test_endpoint_respects_bounded_callback_extent(declared_upper, n_z_acc):
    """Integrate valid tidal queries without sampling beyond the model's redshift extent.

    Parameters
    ----------
    declared_upper : float or None
        Upper redshift of an explicit validation domain; None uses the solver's
        z_max=7 declaration alone. The smaller 4.25 domain tests that callbacks
        need not be defined up to z_max outside their validated range.
    n_z_acc : int
        Two nodes select linear interpolation; 48 nodes exercise the default
        cubic endpoint table and its upper boundary.

    Notes
    -----
    Setup: all three callbacks are bounded interpolants of constant rate 0.3,
    zero exponent and host mass 1e12 in the same arbitrary mass unit as the
    subhalo masses. Any out-of-range callback evaluation raises immediately.
    Queries include an interior accretion epoch, the upper endpoint, and its
    immediately preceding representable float.

    Pass criteria: endpoint, history and direct masses agree with the independent
    solution ma*exp[0.3*(z-z_acc)] to rtol=5e-13. The cached table includes the
    exact upper epoch without extending its callback samples beyond that bound;
    no fallback is needed. Any range error, mass mismatch or fallback fails.
    The tolerance checks numerical integration and interpolation of this simple
    law, not calibration or accuracy for a time-dependent physical model.
    """
    upper = 7.0 if declared_upper is None else declared_upper
    domain = (
        None
        if declared_upper is None
        else PicardDomain((0.0, declared_upper), (-12.0, 0.0))
    )
    solver = PicardMassLossSolver(
        interp1d([0.0, upper], [0.3, 0.3], bounds_error=True),
        interp1d([0.0, upper], [0.0, 0.0], bounds_error=True),
        interp1d([0.0, upper], [1e12, 1e12], bounds_error=True),
        z_min=0.0,
        z_max=7.0,
        state_key=lambda: (0.3, 0.0, 1e12, upper),
        domain=domain,
    )
    masses = np.array([[1e6], [1e8]])
    accretion = np.array([upper / 2, np.nextafter(upper, -np.inf), upper])
    expected = masses * np.exp(-0.3 * accretion)
    np.testing.assert_allclose(
        solver.endpoint(masses, accretion, 0.0, n_z_acc=n_z_acc),
        expected,
        rtol=5e-13,
    )
    for start in accretion:
        redshifts = np.array([start, start / 2, 0.0])
        analytic = masses[None, ...] * np.exp(0.3 * (redshifts - start))[:, None, None]
        np.testing.assert_allclose(solver.history(masses, start, redshifts), analytic, rtol=5e-13)
        np.testing.assert_allclose(solver.direct(masses, start, redshifts), analytic, rtol=5e-13)
    diagnostics = solver.diagnostics()
    assert diagnostics["tables"][0]["z_acc_max"] == upper
    assert diagnostics["fallback_calls"] == 0


@pytest.mark.parametrize("upper", [0.1, 7.0])
@pytest.mark.parametrize("bounded", [False, True])
def test_endpoint_integrates_interval_without_distinct_table_nodes(upper, bounded):
    """Retain physical evolution when adjacent floating-point epochs cannot form a table.

    Parameters
    ----------
    upper : float
        Accretion redshift and callback-domain upper bound; 0.1 and 7 exercise
        different floating-point spacings.
    bounded : bool
        Whether the coefficient interpolants reject evaluation outside [0,upper]
        or permit extrapolation. Both callback contracts must remain usable.

    Notes
    -----
    Setup: accretion is at upper and the final epoch is its immediately preceding
    representable float. The default 48-node table cannot have distinct coordinates
    over that interval. Constant rate 0.3, zero exponent and host mass 1e12 give
    the independent solution ma*exp[0.3*(z-z_acc)].

    Pass criteria: the endpoint emits one explicit representability fallback,
    records both input masses, builds no degenerate table, and agrees with the
    analytic solution, history and direct integration to rtol=5e-14. Invalid
    callback evaluation or an interpolation exception fails. This roundoff-scale
    test establishes safe numerical routing, not scientific convergence for
    finite-time mass loss.
    """
    fill_value = np.nan if bounded else "extrapolate"
    solver = PicardMassLossSolver(
        interp1d([0.0, upper], [0.3, 0.3], bounds_error=bounded, fill_value=fill_value),
        interp1d([0.0, upper], [0.0, 0.0], bounds_error=bounded, fill_value=fill_value),
        interp1d([0.0, upper], [1e12, 1e12], bounds_error=bounded, fill_value=fill_value),
        z_min=0.0,
        z_max=upper,
        state_key=lambda: (0.3, 0.0, 1e12, upper),
    )
    masses = np.array([1e8, 1e10])
    final = np.nextafter(upper, -np.inf)
    expected = masses * np.exp(0.3 * (final - upper))
    with pytest.warns(PicardFallbackWarning, match="representable redshift nodes"):
        endpoint = solver.endpoint(masses, upper, final)
    np.testing.assert_allclose(endpoint, expected, rtol=5e-14)
    np.testing.assert_allclose(solver.history(masses, upper, final), expected, rtol=5e-14)
    np.testing.assert_allclose(solver.direct(masses, upper, final), expected, rtol=5e-14)
    diagnostics = solver.diagnostics()
    assert diagnostics["tables"] == []
    assert diagnostics["fallback_calls"] == 1
    assert diagnostics["fallback_mass_count"] == 2
    assert diagnostics["fallback_events"][0]["reason"] == (
        "insufficient representable redshift nodes for Picard table"
    )
