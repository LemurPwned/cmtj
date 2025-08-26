import pytest
import numpy as np
from cmtj.models import LayerSB, VectorObj, Solver
from cmtj.utils import mu0

import json
import os


def angular_difference(angle1, angle2, period=2 * np.pi):
    """
    Calculate the minimum angular difference between two angles, accounting for periodicity.

    Args:
        angle1, angle2: angles in radians
        period: period of the angle (2*pi for phi, pi for some theta cases)

    Returns:
        Minimum angular difference
    """
    diff = np.abs(angle1 - angle2)
    # Check if the difference is smaller when considering the periodic boundary
    diff_periodic = period - diff
    return np.minimum(diff, diff_periodic)


def assert_angles_close(
    actual, expected, rtol=1e-1, atol=1e-1, period=2 * np.pi, err_msg=""
):
    """
    Assert that two angles are close, accounting for periodicity.

    Args:
        actual, expected: angles to compare
        rtol, atol: relative and absolute tolerance
        period: period of the angle (2*pi for phi, no period handling for theta in [0,pi])
        err_msg: error message
    """
    if period is None:
        # No periodicity handling (for theta in [0, π])
        np.testing.assert_allclose(
            actual, expected, rtol=rtol, atol=atol, err_msg=err_msg
        )
    else:
        # Handle periodicity
        diff = angular_difference(actual, expected, period)
        tolerance = atol + rtol * np.abs(expected)
        assert (
            diff <= tolerance
        ), f"{err_msg}. Angular difference: {diff:.6f}, tolerance: {tolerance:.6f}, actual: {actual:.6f}, expected: {expected:.6f}"


@pytest.fixture
def tutorial_two_layer_system():
    """Create the two-layer system from the SBModel tutorial"""
    Ms1 = 1.0 / mu0  # saturation magnetisation in A/m
    Ms2 = 1.2 / mu0

    layerA = LayerSB(
        _id=0,
        thickness=1e-9,
        Kv=VectorObj(
            np.deg2rad(0.0), np.deg2rad(0), 1e1
        ),  # only phi angle counts for Kv
        Ks=3e4,
        Ms=Ms1,
    )

    layerB = LayerSB(
        _id=1,
        thickness=1.3e-9,
        Kv=VectorObj(np.deg2rad(0.0), np.deg2rad(0), 1e4),
        Ks=1e1,
        Ms=Ms2,
    )

    return layerA, layerB


# Expected results for different solver configurations
# These are reference values from a specific field point for consistency testing
EXPECTED_FREQUENCIES = {
    "analytical_no_lu": {"freq_count": 4, "min_freq": 0.0, "max_freq": 50.0},
    "analytical_lu": {"freq_count": 4, "min_freq": 0.0, "max_freq": 50.0},
    "numerical_no_lu": {"freq_count": 4, "min_freq": 0.0, "max_freq": 50.0},
    "numerical_lu": {"freq_count": 4, "min_freq": 0.0, "max_freq": 50.0},
}

EXPECTED_COORDINATES = {
    "theta_1_range": (0.0, np.pi),
    "phi_1_range": (0.0, 2 * np.pi),
    "theta_2_range": (0.0, np.pi),
    "phi_2_range": (0.0, 2 * np.pi),
}


# Add pytest marks for better control
@pytest.mark.fast
@pytest.mark.parametrized
@pytest.mark.parametrize(
    "prefer_numerical_roots,use_LU_decomposition,config_name",
    [
        (False, False, "analytical_no_lu"),
        (False, True, "analytical_lu"),
        (True, False, "numerical_no_lu"),
        (True, True, "numerical_lu"),
    ],
)
def test_sb_solver_configurations(
    tutorial_two_layer_system, prefer_numerical_roots, use_LU_decomposition, config_name
):
    """Test SB solver with different configurations for consistency of frequencies and coordinates"""
    layerA, layerB = tutorial_two_layer_system

    # Test at a specific field point for consistency
    test_field = 100e3  # 100 kA/m
    current_position = [
        np.deg2rad(89),
        np.deg2rad(0.1),
        np.deg2rad(180),
        np.deg2rad(0.1),
    ]

    solver = Solver(
        layers=[layerA, layerB],
        J1=[1e-4],
        J2=[0.0],
        H=VectorObj(np.deg2rad(89), np.deg2rad(0.1), test_field),
        prefer_numerical_roots=prefer_numerical_roots,
        use_LU_decomposition=use_LU_decomposition,
    )

    # Solve the system
    eq, frequencies = solver.solve(init_position=current_position, perturbation=0)
    t1, p1, t2, p2 = eq

    # Test that coordinates are within valid ranges
    assert 0 <= t1 <= np.pi, f"theta_1 out of range for {config_name}: {t1}"
    assert 0 <= p1 <= 2 * np.pi, f"phi_1 out of range for {config_name}: {p1}"
    assert 0 <= t2 <= np.pi, f"theta_2 out of range for {config_name}: {t2}"
    assert 0 <= p2 <= 2 * np.pi, f"phi_2 out of range for {config_name}: {p2}"

    # Test that frequencies are reasonable (positive, finite, within expected range)
    assert len(frequencies) > 0, f"No frequencies returned for {config_name}"
    for freq in frequencies:
        assert np.isfinite(freq), f"Non-finite frequency for {config_name}: {freq}"
        assert freq >= 0, f"Negative frequency for {config_name}: {freq}"
        assert freq <= 100, f"Unreasonably high frequency for {config_name}: {freq} GHz"

    # Test that frequencies are in GHz (as documented)
    max_freq = max(frequencies) if np.all(frequencies) else 0
    assert (
        max_freq < 100
    ), f"Frequencies seem too high for {config_name}, might not be in GHz: {max_freq}"


@pytest.mark.fast
def test_solver_reproducibility(tutorial_two_layer_system):
    """Test that solver produces reproducible results with same configuration"""
    layerA, layerB = tutorial_two_layer_system

    test_configs = [
        (False, False),
        (False, True),
        (True, False),
        (True, True),
    ]

    for prefer_numerical, use_LU in test_configs:
        # Run the same solver twice
        results1 = []
        results2 = []

        for run in range(2):
            solver = Solver(
                layers=[layerA, layerB],
                J1=[1e-4],
                J2=[0.0],
                H=VectorObj(np.deg2rad(89), np.deg2rad(0.1), 100e3),
                prefer_numerical_roots=prefer_numerical,
                use_LU_decomposition=use_LU,
            )

            init_pos = [
                np.deg2rad(89),
                np.deg2rad(0.1),
                np.deg2rad(180),
                np.deg2rad(0.1),
            ]
            result = solver.solve(init_position=init_pos)

            if run == 0:
                results1 = result
            else:
                results2 = result

        # Compare results
        coords1, freqs1 = results1
        coords2, freqs2 = results2

        # Coordinates should be identical - using angular comparison for phi angles
        t1_1, p1_1, t2_1, p2_1 = coords1
        t1_2, p1_2, t2_2, p2_2 = coords2

        # Compare theta angles (no periodicity needed as they're in [0, π])
        assert_angles_close(
            t1_1,
            t1_2,
            rtol=1e-2,
            atol=1e-2,
            period=None,
            err_msg=f"theta_1 not reproducible for config: numerical={prefer_numerical}, LU={use_LU}",
        )
        assert_angles_close(
            t2_1,
            t2_2,
            rtol=1e-2,
            atol=1e-2,
            period=None,
            err_msg=f"theta_2 not reproducible for config: numerical={prefer_numerical}, LU={use_LU}",
        )

        # Compare phi angles (with 2π periodicity)
        assert_angles_close(
            p1_1,
            p1_2,
            rtol=1e-2,
            atol=1e-2,
            period=2 * np.pi,
            err_msg=f"phi_1 not reproducible for config: numerical={prefer_numerical}, LU={use_LU}",
        )
        assert_angles_close(
            p2_1,
            p2_2,
            rtol=1e-2,
            atol=1e-2,
            period=2 * np.pi,
            err_msg=f"phi_2 not reproducible for config: numerical={prefer_numerical}, LU={use_LU}",
        )

        # Frequencies should be identical (sorted to handle order differences)
        freqs1_sorted = sorted(freqs1)
        freqs2_sorted = sorted(freqs2)

        np.testing.assert_allclose(
            freqs1_sorted,
            freqs2_sorted,
            rtol=1e-2,
            atol=1e-2,
            err_msg=f"Frequencies not reproducible for config: numerical={prefer_numerical}, LU={use_LU}",
        )


@pytest.fixture
def reference_data():
    """Load reference data once for all tests"""
    assets_dir = os.path.join(os.path.dirname(__file__), "assets")
    ref_file = os.path.join(assets_dir, "bilayer_sb_fmr_8_point_sweep.json")

    with open(ref_file, "r") as f:
        return json.load(f)


# Mark the reference data test as potentially slower
@pytest.mark.slow
@pytest.mark.parametrized
@pytest.mark.parametrize(
    "prefer_numerical,use_LU,config_name",
    [
        (True, False, "analytical_no_lu"),
        (True, True, "numerical_lu"),
    ],
)
def test_sb_solver_against_reference_data(
    tutorial_two_layer_system, reference_data, prefer_numerical, use_LU, config_name
):
    """Test SB solver configurations against stashed reference data for exact validation"""
    layerA, layerB = tutorial_two_layer_system

    # Extract reference field values and expected results
    ref_Hmag = reference_data["Hmag"]
    ref_frequencies = reference_data["frequency"]
    ref_theta_1 = reference_data["theta_1"]
    ref_phi_1 = reference_data["phi_1"]
    ref_theta_2 = reference_data["theta_2"]
    ref_phi_2 = reference_data["phi_2"]

    # Create field points from reference data (unique values)
    unique_fields = sorted(list(set(ref_Hmag)))
    # Initialize starting position
    current_position = [
        np.deg2rad(89),
        np.deg2rad(0.1),
        np.deg2rad(180),
        np.deg2rad(0.1),
    ]

    # Store results for this configuration
    config_results = {
        "Hmag": [],
        "theta_1": [],
        "phi_1": [],
        "theta_2": [],
        "phi_2": [],
        "frequency": [],
    }

    # Run field sweep
    for Hmag in unique_fields:
        solver = Solver(
            layers=[layerA, layerB],
            J1=[1e-4],
            J2=[0.0],
            H=VectorObj(np.deg2rad(89), np.deg2rad(0.1), Hmag),
            prefer_numerical_roots=prefer_numerical,
            use_LU_decomposition=use_LU,
        )

        (t1, p1, t2, p2), frequencies = solver.solve(
            init_position=current_position, perturbation=0
        )

        # Store results
        config_results["theta_1"].append(t1)
        config_results["phi_1"].append(p1)
        config_results["theta_2"].append(t2)
        config_results["phi_2"].append(p2)
        for frequency in frequencies:
            config_results["frequency"].append(frequency)

        # Use previous solution as initial guess for next iteration
        current_position = [t1, p1, t2, p2]
    for k in config_results:
        config_results[k] = np.asarray(config_results[k]).tolist()

    assert len(config_results["frequency"]) == len(ref_frequencies)
    assert len(config_results["theta_1"]) == len(ref_theta_1)
    assert len(config_results["phi_1"]) == len(ref_phi_1)
    assert len(config_results["theta_2"]) == len(ref_theta_2)
    assert len(config_results["phi_2"]) == len(ref_phi_2)

    for k in ("frequency", "theta_1", "phi_1", "theta_2", "phi_2"):
        assert np.allclose(
            config_results[k], reference_data[k], rtol=1e-1, atol=1e-1
        ), f"Mismatch in {k} for {config_name}"
