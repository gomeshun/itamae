"""Independent public-boundary cases that bypass population execution."""

import json
import warnings

import numpy as np
import pytest

from itamae.cosmology import NativeFlatLCDM
from itamae.halo import invert_nfw_mass_function
from itamae.types import CatalogMetadata, WeightedSubhaloCatalog
from itamae.units import NativeUnits


def make_catalog(column=(1.0, 2.0), weights=None):
    return WeightedSubhaloCatalog(
        columns={"physical": np.asarray(column)},
        weights={"weight_base": np.ones(2)} if weights is None else weights,
        metadata=CatalogMetadata(model_identifier="analytic-test", backend_identifier="numpy"),
    )


@pytest.mark.parametrize("column", ([1.0, np.nan], [1.0, np.inf], [1j, 2j], ["a", "b"]))
def test_direct_catalog_rejects_nonreal_or_nonfinite_columns(column):
    with pytest.raises(ValueError, match="column.*real|column.*finite"):
        make_catalog(column)


def test_direct_catalog_rejects_complex_weights_before_cast():
    with pytest.raises(ValueError, match="Weight.*real"):
        make_catalog(weights={"weight_base": np.array([1 + 1j, 2 + 0j])})


def test_factor_product_overflow_is_an_error():
    with pytest.raises(ValueError, match="product.*finite|product.*overflow"):
        make_catalog(
            weights={"weight_base": np.full(2, 1e200), "weight_concentration": np.full(2, 1e200)}
        )


@pytest.mark.parametrize("values", ([1.0, np.nan], [1.0, np.inf], [1 + 1j, 2 + 0j]))
def test_weighted_sum_rejects_invalid_values(values):
    with pytest.raises(ValueError, match="Values.*real|Values.*finite"):
        make_catalog().weighted_sum(values)


def test_weighted_sum_overflow_is_an_error():
    with pytest.raises(ValueError, match="sum.*finite|sum.*overflow"):
        make_catalog().weighted_sum([1e308, 1e308])


def test_empty_density_histogram_has_no_silent_nan():
    empty = make_catalog().select(np.array([False, False]))
    counts, _ = empty.weighted_histogram("physical", bins=[0, 1, 2])
    np.testing.assert_array_equal(counts, [0, 0])
    with pytest.raises(ValueError, match="density|histogram"):
        empty.weighted_histogram("physical", bins=[0, 1, 2], density=True)


def test_npz_rejects_duplicate_names_in_manifest(tmp_path):
    manifest = {
        "archive_schema": "itamae-weighted-catalog-npz:1.0",
        "catalog_metadata": dict(make_catalog().metadata),
        "columns": ["physical", "physical"],
        "weights": ["weight_base"],
    }
    path = tmp_path / "duplicate.npz"
    np.savez(
        path,
        manifest_json=json.dumps(manifest),
        column_0=np.ones(2),
        column_1=np.zeros(2),
        weight_0=np.ones(2),
    )
    with pytest.raises(ValueError, match="duplicate|unique"):
        WeightedSubhaloCatalog.from_npz(path)


def test_npz_rejects_nonmapping_manifest(tmp_path):
    path = tmp_path / "list.npz"
    np.savez(path, manifest_json="[]")
    with pytest.raises(ValueError, match="manifest"):
        WeightedSubhaloCatalog.from_npz(path)


def test_npz_checks_physical_arrays_on_load(tmp_path):
    catalog = make_catalog()
    path = tmp_path / "corrupt.npz"
    catalog.to_npz(path)
    with np.load(path, allow_pickle=False) as archive:
        payload = dict(archive)
    payload["column_0"] = np.array([1.0, np.nan])
    np.savez(path, **payload)
    with pytest.raises(ValueError, match="column.*finite"):
        WeightedSubhaloCatalog.from_npz(path)


@pytest.mark.parametrize("value", ([1 + 1j], [np.inf], [np.nan]))
def test_native_units_reject_unrepresentable_inputs_and_outputs(value):
    with pytest.raises(ValueError, match="real|finite"):
        NativeUnits().to_internal(value, "mass")
    with pytest.raises(ValueError, match="real|finite"):
        NativeUnits().from_internal(value)


def test_astropy_units_reject_complex_before_cast():
    u = pytest.importorskip("astropy.units")
    from itamae.units.astropy import AstropyUnits

    with pytest.raises(ValueError, match="real"):
        AstropyUnits().to_internal(np.array([1 + 1j]) * u.Msun, "mass")
    with pytest.raises(ValueError, match="real"):
        AstropyUnits().from_internal(np.array([1 + 1j]), u.kpc)


@pytest.mark.parametrize(
    "parameter,value",
    [
        ("omega_m0", 0.0),
        ("omega_m0", -0.1),
        ("omega_m0", 1.1),
        ("omega_m0", np.nan),
        ("h", 0.0),
        ("h", -0.1),
        ("h", np.inf),
    ],
)
def test_native_cosmology_rejects_invalid_parameters(parameter, value):
    with pytest.raises(ValueError, match=parameter):
        NativeFlatLCDM(**{parameter: value})


@pytest.mark.parametrize("z", (-1.0, -2.0, np.nan, np.inf, 1 + 1j))
def test_native_cosmology_rejects_invalid_redshift(z):
    model = NativeFlatLCDM()
    for method in (
        model.e2,
        model.H,
        model.rho_m,
        model.rho_crit,
        model.omega_m,
        model.growth_factor,
        model.collapse_threshold,
        model.cosmic_time,
        model.lookback_time,
    ):
        with pytest.raises(ValueError, match="[Rr]edshift"):
            method(z)


def test_native_cosmology_identifies_distinct_float_parameters():
    assert NativeFlatLCDM(h=0.674).identifier != NativeFlatLCDM(h=0.67400000001).identifier


def test_background_density_overflow_is_explicit():
    for method in (NativeFlatLCDM().rho_m, NativeFlatLCDM().rho_crit):
        with pytest.raises(ValueError, match="finite numerical background"):
            method(1e100)


def test_subnormal_nfw_target_has_a_finite_warning_free_bracket():
    value = 1e-310
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        radius = invert_nfw_mass_function(value)
    # At this radius every higher-order series term is below float precision.
    assert radius == pytest.approx(np.sqrt(2 * value), rel=1e-10, abs=0.0)


def test_einstein_de_sitter_limits_and_future_redshift():
    model = NativeFlatLCDM(omega_m0=1.0)
    z = np.array([-0.5, 0.0, 1.0, 3.0])
    np.testing.assert_allclose(model.growth_factor(z), 1 / (1 + z), rtol=2e-15)
    np.testing.assert_allclose(model.omega_m(z), np.ones(4), rtol=2e-15)
    hubble_in_gyr = 100 * model.h * 1.0227121650537077e-3
    np.testing.assert_allclose(
        model.cosmic_time(z), 2 / (3 * hubble_in_gyr) * (1 + z) ** -1.5, rtol=2e-10
    )
