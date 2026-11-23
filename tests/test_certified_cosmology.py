import json
from fractions import Fraction
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import patch

import pytest

from itamae.cosmology.certified import evaluate_mass, evaluate_scalar, evaluate_table


def reply(data, code=0):
    return CompletedProcess([], code, json.dumps(data), "")


def mock_project(tmp_path: Path):
    (tmp_path / "scripts").mkdir()
    (tmp_path / "lean-toolchain").write_text("leanprover/lean4:v4.36.0")
    for name in ("spectrum_table.py", "sashimi.py", "mass_variance.py"):
        (tmp_path / "scripts" / name).touch()
    return tmp_path


def test_scalar_wheel_with_explicit_project_and_replay(tmp_path):
    project = mock_project(tmp_path)
    with patch(
        "itamae.cosmology.certified.subprocess.run",
        return_value=reply(
            {
                "status": "enclosed",
                "unit": "km/s/Mpc",
                "tolerance": "1/10",
                "width": "1/20",
                "kernel_replayed": True,
                "interval": {"lower": "1/3", "upper": "23/60"},
            }
        ),
    ) as run:
        s = evaluate_scalar("hubble", {"z": "2.5"}, project=project, tolerance="0.10", replay=True)
    assert s.lower == Fraction(1, 3)
    assert s.upper == Fraction(23, 60)
    assert s.width == Fraction(1, 20)
    assert s.kernel_replayed
    argv = run.call_args.args[0]
    assert "z=5/2" in argv and "--replay" in argv
    assert run.call_args.kwargs["cwd"] == project


def test_table_returns_exact_bounds_and_source_metadata():
    data = {
        "status": "enclosed",
        "unit": "Mpc^3",
        "absolute_width_budget": "1/10",
        "lower": "2",
        "upper": "41/20",
        "width": "1/20",
        "table": {"sha256": "source"},
    }
    with patch("itamae.cosmology.certified.subprocess.run", return_value=reply(data)) as run:
        s = evaluate_table(
            "example.dat",
            layout="camb",
            h="0.7",
            k="0.25",
            tolerance="0.1",
            interpolation="linear",
            transfer="cdm",
            placement="continuous",
        )
    assert s.width == Fraction(1, 20) and not s.kernel_replayed
    assert s.metadata["table"] == {"sha256": "source"}
    argv = run.call_args.args[0]
    assert str(Path("example.dat").resolve()) in argv and "7/10" in argv
    assert run.call_args.kwargs["cwd"].name == "physicslean"


def test_unresolved_has_no_numeric_value():
    data = {"status": "unresolved", "absolute_width_budget": "1/1000", "unit": "1"}
    with patch("itamae.cosmology.certified.subprocess.run", return_value=reply(data, 2)):
        s = evaluate_table(
            "p.dat",
            layout="canonical",
            k=1,
            k_max=2,
            tolerance="1e-3",
            interpolation="log-log",
            placement="at-knots",
            transfer="cdm",
        )
    assert s.status == "unresolved" and s.tolerance == Fraction(1, 1000)
    assert s.upper is None and s.lower is None and s.width is None


def test_table_can_use_wdm_particle_parameters():
    data = {"status": "unresolved", "absolute_width_budget": "1/1000", "unit": "Mpc^3"}
    with patch("itamae.cosmology.certified.subprocess.run", return_value=reply(data, 2)) as run:
        s = evaluate_table(
            "p.dat",
            layout="canonical",
            k=1,
            tolerance=".001",
            interpolation="linear",
            placement="at-knots",
            transfer="wdm-q10",
            wdm_mass_kev=3,
            omega_wdm=".25",
            thermal_h=".7",
        )
    argv = run.call_args.args[0]
    assert argv[argv.index("--wdm-mass-kev") :] == [
        "--wdm-mass-kev",
        "3",
        "--omega-wdm",
        "1/4",
        "--thermal-h",
        "7/10",
    ]
    assert s.status == "unresolved"


def test_normalized_table_replay_carries_anchor_metadata():
    args = dict(
        layout="canonical",
        k="1.5",
        tolerance="0.1",
        interpolation="linear",
        transfer="cdm",
        placement="continuous",
        sigma8="0.8",
        anchor_radius="0.25",
        anchor_k_min="1",
        anchor_k_max="2",
        anchor_depth=4,
        replay=True,
    )
    normalization = {"target_sigma": "4/5", "radius_Mpc": "1/4"}
    data = {
        "status": "enclosed",
        "unit": "Mpc^3",
        "absolute_width_budget": "1/10",
        "lower": "2",
        "upper": "41/20",
        "width": "1/20",
        "normalization": normalization,
        "kernel_replayed": True,
    }
    with patch("itamae.cosmology.certified.subprocess.run", return_value=reply(data)) as run:
        s = evaluate_table("input.dat", **args)
    argv = run.call_args.args[0]
    assert argv[argv.index("--sigma8") :] == [
        "--sigma8",
        "4/5",
        "--anchor-radius",
        "1/4",
        "--anchor-k-min",
        "1",
        "--anchor-k-max",
        "2",
        "--anchor-depth",
        "4",
        "--replay",
    ]
    assert s.kernel_replayed and s.metadata["normalization"] == normalization
    data["kernel_replayed"] = False
    with (
        patch("itamae.cosmology.certified.subprocess.run", return_value=reply(data)),
        pytest.raises(RuntimeError, match="did not confirm"),
    ):
        evaluate_table("input.dat", **args)


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"anchor_depth": 2}, "require --sigma8"),
        ({"sigma8": "4/5"}, "needs --anchor-radius"),
        ({"sigma8": "4/5", "anchor_radius": 0, "anchor_k_min": 2}, "anchor band"),
        ({"sigma8": "4/5", "anchor_radius": 0, "anchor_depth": 0.5}, "nonnegative integer"),
        ({"tolerance": 0.3}, "decimal/fraction string"),
    ],
)
def test_invalid_table_options_before_lean(tmp_path, options, message):
    table = tmp_path / "power.dat"
    table.write_text("1 2\n2 2\n")
    args = {
        "layout": "canonical",
        "k": 1,
        "tolerance": "0.1",
        "interpolation": "linear",
        "transfer": "cdm",
        "placement": "continuous",
        **options,
    }
    with pytest.raises(ValueError, match=message):
        evaluate_table(table, **args)


def test_malformed_driver_output_is_not_a_result():
    with (
        patch(
            "itamae.cosmology.certified.subprocess.run",
            return_value=reply(
                {
                    "status": "enclosed",
                    "tolerance": "0",
                    "unit": "1",
                    "width": "20",
                    "interval": {"lower": "5", "upper": "25"},
                }
            ),
        ),
        pytest.raises(RuntimeError, match="violates the requested width"),
    ):
        evaluate_scalar("growth", {"z": 0}, tolerance=0)


def test_mass_variance_derivative_carries_growth_convention(tmp_path):
    project = mock_project(tmp_path)
    options = dict(
        project=project,
        layout="canonical",
        quantity="sharp-derivative",
        interpolation="linear",
        placement="continuous",
        transfer="cdm",
        density="3e10",
        mass="1e14",
        k_min=".6",
        k_max="1.1",
        tolerance="1e-12",
        derivative_growth="legacy-f-linear",
        omega_m=".315",
        redshift=2,
    )
    data = {
        "status": "enclosed",
        "unit": "1/Msun",
        "tolerance": "1e-12",
        "lower": "-2e-15",
        "upper": "-1e-15",
        "width": "1e-15",
        "growth_power": 1,
        "growth_model": "CPT",
        "redshift": "2",
    }
    with patch("itamae.cosmology.certified.subprocess.run", return_value=reply(data)) as run:
        result = evaluate_mass("p.dat", **options)
    argv = run.call_args.args[0]
    assert "--growth" not in argv
    for flag, expected in [("--omega-m", "63/200"), ("--derivative-growth", "legacy-f-linear")]:
        assert argv[argv.index(flag) + 1] == expected
    assert result.upper == Fraction(-1, 10**15)
    assert result.unit == "1/Msun" and result.metadata["growth_power"] == 1
    with pytest.raises(ValueError, match="either growth or redshift"):
        evaluate_mass("p.dat", **options, growth=".5")


def test_mass_preserves_w_calibration_and_particle_units(tmp_path):
    project = mock_project(tmp_path)
    options = dict(
        project=project,
        layout="canonical",
        quantity="sharp",
        interpolation="linear",
        placement="continuous",
        transfer="wdm-q5",
        density="3e10",
        mass="1e14",
        k_min=".6",
        k_max="1.1",
        tolerance=".01",
        growth=1,
        omega_wdm=".2231",
        thermal_h=".7",
    )
    data = {"status": "unresolved", "unit": "1", "tolerance": "1/100"}
    with patch("itamae.cosmology.certified.subprocess.run", return_value=reply(data, 2)) as run:
        for calibration in ("w-family-legacy", "w-main", "sashimi-w", 1):
            evaluate_mass("p.dat", **options, calibration=calibration, wdm_mass_kev="1.5")
            argv = run.call_args.args[0]
            assert argv[argv.index("--calibration") + 1] == str(calibration)
            assert argv[argv.index("--wdm-mass-kev") + 1] == "3/2" and "--parameter" not in argv
        with pytest.raises(ValueError, match="wdm_mass_kev"):
            evaluate_mass("p.dat", **options, parameter=".045")
        assert run.call_count == 4


def test_explicit_missing_project_and_inexact_inputs(tmp_path):
    with pytest.raises(FileNotFoundError, match="project="):
        evaluate_scalar("growth", {"z": 0}, project=tmp_path)
    with pytest.raises(ValueError):
        evaluate_scalar("growth", {"z": 0.3}, project=tmp_path)
    with pytest.raises(ValueError):
        evaluate_scalar("growth", {"z": "0.3"}, bits=True)
    with pytest.raises(ValueError, match="strictly positive"):
        evaluate_scalar("growth", {"z": 0, "h": 0})
