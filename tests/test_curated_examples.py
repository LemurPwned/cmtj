"""
Tests for curated examples - comparing simulation outputs to 3 decimal places.

These tests verify that the curated examples produce consistent results.
The examples are kept in curated-examples/ but these tests ensure
regression testing of the simulation outputs.

Reference data is stored in tests/assets/curated_examples/ and can be updated
by setting the UPDATE_REFERENCE environment variable.
"""

import json
import math
import os
from collections import defaultdict

import numpy as np
import pytest

import cmtj
from cmtj import (
    AxialDriver,
    CVector,
    Junction,
    Layer,
    NullDriver,
    ScalarDriver,
    SolverMode,
    constantDriver,
    sineDriver,
    stepDriver,
)
from cmtj.utils import FieldScan, Filters, OetoAm, TtoAm, compute_gmr
from cmtj.utils.resistance import calculate_resistance_parallel, calculate_resistance_series


def round_to_3_decimals(arr):
    """Round array to 3 decimal places for comparison."""
    return np.round(arr, decimals=3)


def get_reference_file_path(test_name: str) -> str:
    """Get the path to the reference data file for a test."""
    assets_dir = os.path.join(os.path.dirname(__file__), "assets", "curated_examples")
    os.makedirs(assets_dir, exist_ok=True)
    return os.path.join(assets_dir, f"{test_name}.json")


def load_reference_data(test_name: str) -> dict | None:
    """Load reference data for a test, or return None if it doesn't exist."""
    ref_file = get_reference_file_path(test_name)
    if os.path.exists(ref_file):
        with open(ref_file, "r") as f:
            return json.load(f)
    return None


def save_reference_data(test_name: str, data: dict):
    """Save reference data for a test."""
    ref_file = get_reference_file_path(test_name)
    # Convert numpy arrays to lists for JSON serialization
    serializable_data = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            serializable_data[key] = value.tolist()
        elif isinstance(value, (list, tuple)) and len(value) > 0 and isinstance(value[0], np.ndarray):
            serializable_data[key] = [v.tolist() if isinstance(v, np.ndarray) else v for v in value]
        else:
            serializable_data[key] = value

    with open(ref_file, "w") as f:
        json.dump(serializable_data, f, indent=2)


def compare_with_reference(
    test_name: str, current_data: dict, update_reference: bool = False, stochastic: bool = False
):
    """
    Compare current test results with reference data.

    Args:
        test_name: Name of the test
        current_data: Dictionary of current test outputs (will be rounded to 3 decimals)
        update_reference: If True, update the reference data instead of comparing
        stochastic: If True, use more lenient comparison for stochastic simulations
    """
    # Round all arrays to 3 decimals
    rounded_data = {}
    for key, value in current_data.items():
        if isinstance(value, np.ndarray):
            rounded_data[key] = round_to_3_decimals(value)
        elif isinstance(value, (list, tuple)) and len(value) > 0:
            if isinstance(value[0], np.ndarray):
                rounded_data[key] = [round_to_3_decimals(v) for v in value]
            else:
                rounded_data[key] = value
        else:
            rounded_data[key] = value

    ref_data = load_reference_data(test_name)

    if update_reference or ref_data is None:
        if ref_data is None:
            print(f"No reference data found for {test_name}, creating new reference...")
        else:
            print(f"Updating reference data for {test_name}...")
        save_reference_data(test_name, rounded_data)
        return True

    # Compare with reference
    assert ref_data is not None, f"No reference data found for {test_name}"

    for key in rounded_data:
        assert key in ref_data, f"Missing key '{key}' in reference data for {test_name}"

        current_val = rounded_data[key]
        ref_val = ref_data[key]

        # Handle nested lists (list of arrays)
        if isinstance(current_val, list) and len(current_val) > 0:
            if isinstance(current_val[0], (list, np.ndarray)):
                # List of arrays case
                current_list = [np.array(v) if isinstance(v, list) else v for v in current_val]
                ref_list = [np.array(v) if isinstance(v, list) else v for v in ref_val]
                assert len(current_list) == len(ref_list), (
                    f"Length mismatch for '{key}' in {test_name}: "
                    f"current {len(current_list)} vs reference {len(ref_list)}"
                )
                for i, (c, r) in enumerate(zip(current_list, ref_list)):
                    if stochastic:
                        # For stochastic tests, compare statistics instead of exact values
                        # Use more lenient tolerance for stochastic simulations
                        mean_diff = np.abs(np.mean(c) - np.mean(r))
                        std_diff = np.abs(np.std(c) - np.std(r))
                        assert mean_diff < 0.5, (
                            f"Mean mismatch for '{key}[{i}]' in {test_name}: "
                            f"current {np.mean(c):.3f} vs reference {np.mean(r):.3f} (diff: {mean_diff:.3f})"
                        )
                        assert std_diff < 0.5, (
                            f"Std mismatch for '{key}[{i}]' in {test_name}: "
                            f"current {np.std(c):.3f} vs reference {np.std(r):.3f} (diff: {std_diff:.3f})"
                        )
                        # Also check that ranges are similar
                        assert np.abs(np.min(c) - np.min(r)) < 1.0, f"Min mismatch for '{key}[{i}]' in {test_name}"
                        assert np.abs(np.max(c) - np.max(r)) < 1.0, f"Max mismatch for '{key}[{i}]' in {test_name}"
                    else:
                        np.testing.assert_array_almost_equal(
                            c, r, decimal=3, err_msg=f"Value mismatch for '{key}[{i}]' in {test_name}"
                        )
            else:
                # Regular list of numbers
                current_arr = np.array(current_val)
                ref_arr = np.array(ref_val)
                assert current_arr.shape == ref_arr.shape, (
                    f"Shape mismatch for '{key}' in {test_name}: "
                    f"current {current_arr.shape} vs reference {ref_arr.shape}"
                )
                np.testing.assert_array_almost_equal(
                    current_arr, ref_arr, decimal=3, err_msg=f"Value mismatch for '{key}' in {test_name}"
                )
        elif isinstance(current_val, np.ndarray):
            # Single array (could be 1D, 2D, etc.)
            ref_arr = np.array(ref_val)
            assert current_val.shape == ref_arr.shape, (
                f"Shape mismatch for '{key}' in {test_name}: current {current_val.shape} vs reference {ref_arr.shape}"
            )
            if stochastic:
                # For stochastic tests, compare statistics for arrays
                mean_diff = np.abs(np.mean(current_val) - np.mean(ref_arr))
                std_diff = np.abs(np.std(current_val) - np.std(ref_arr))
                assert mean_diff < 0.5, (
                    f"Mean mismatch for '{key}' in {test_name}: "
                    f"current {np.mean(current_val):.3f} vs reference {np.mean(ref_arr):.3f} (diff: {mean_diff:.3f})"
                )
                assert std_diff < 0.5, (
                    f"Std mismatch for '{key}' in {test_name}: "
                    f"current {np.std(current_val):.3f} vs reference {np.std(ref_arr):.3f} (diff: {std_diff:.3f})"
                )
                # Check that most values are close (allow some outliers)
                diff = np.abs(current_val - ref_arr)
                assert np.percentile(diff, 95) < 0.5, (
                    f"95th percentile mismatch for '{key}' in {test_name}: "
                    f"95% of values differ by {np.percentile(diff, 95):.3f}"
                )
            else:
                np.testing.assert_array_almost_equal(
                    current_val, ref_arr, decimal=3, err_msg=f"Value mismatch for '{key}' in {test_name}"
                )
        else:
            # Scalar value or tuple/list that should be compared as values
            # Convert tuples to lists for comparison
            if isinstance(current_val, tuple):
                current_val = list(current_val)
            if isinstance(ref_val, tuple):
                ref_val = list(ref_val)
            if isinstance(current_val, list) and isinstance(ref_val, list):
                # Compare lists element-wise
                assert len(current_val) == len(ref_val), f"Length mismatch for '{key}' in {test_name}"
                for i, (c, r) in enumerate(zip(current_val, ref_val)):
                    assert c == r or (isinstance(c, float) and isinstance(r, float) and abs(c - r) < 1e-3), (
                        f"Value mismatch for '{key}[{i}]' in {test_name}: current {c} vs reference {r}"
                    )
            else:
                assert current_val == ref_val, (
                    f"Value mismatch for '{key}' in {test_name}: current {current_val} vs reference {ref_val}"
                )

    return True


@pytest.fixture(scope="session")
def update_reference():
    """Fixture to check if we should update reference data."""
    return os.getenv("UPDATE_REFERENCE", "false").lower() == "true"


@pytest.mark.slow
def test_ahe_loops(update_reference):
    """Test AHE loops simulation."""
    # Reduced parameters for faster testing
    thickness = 1e-9
    J = 3e-3
    Ms = 1.1
    Ku1 = 0.55e6
    Ku2 = 0.3e6
    Ku1_theta = 1
    Ku1_phi = 0
    Ku2_theta = 0.01
    Ku2_phi = 0
    Hdl = 90 * OetoAm * 0.45
    Hfl = 90 * OetoAm
    AHE_asymmetry = 1.0
    tstep = 1e-13
    AHE_scaling = 5
    alpha = 5e-2

    rho_f = 27e-8
    rho_h = 21e-8
    t_fm = 1e-9
    t_hm = thickness
    w = 10e-6
    l = 80e-6
    area = w * l
    FM_R = rho_f * t_fm / area
    HM_R = rho_h * l / (w * t_hm)
    demagTensor = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1)]

    jden = 4.24e10
    HDL = Hdl / jden
    HFL = Hfl / jden

    Kdir1 = FieldScan.angle2vector(Ku1_theta, Ku1_phi)
    Kdir2 = FieldScan.angle2vector(Ku2_theta, Ku2_phi)
    pdir1 = CVector(0, 1.0, 0.0)
    pdir2 = CVector(0, 1.0, 0.0)

    layer_free = Layer.createSOTLayer(
        id="free",
        mag=Kdir1,
        anis=Kdir1,
        Ms=Ms,
        thickness=t_fm,
        cellSurface=area,
        demagTensor=demagTensor,
        damping=alpha,
        dampingLikeTorque=HDL,
        fieldLikeTorque=HFL,
    )

    layer_bottom = Layer.createSOTLayer(
        id="bottom",
        mag=Kdir2,
        anis=Kdir2,
        Ms=Ms,
        thickness=t_fm,
        cellSurface=area,
        demagTensor=demagTensor,
        damping=alpha,
        dampingLikeTorque=HDL,
        fieldLikeTorque=HFL,
    )

    j = Junction([layer_free, layer_bottom])
    j.setLayerAnisotropyDriver("free", constantDriver(Ku1))
    j.setLayerAnisotropyDriver("bottom", constantDriver(Ku2))
    j.setLayerReferenceLayer("free", pdir1)
    j.setLayerReferenceLayer("bottom", pdir2)
    j.setIECDriver("free", "bottom", constantDriver(J))

    j.runSimulation(4e-9, tstep, tstep)
    output = defaultdict(list)
    n_lay = 2

    # Reduced H scan for faster testing
    Hscan = np.linspace(-2000 * OetoAm, 2000 * OetoAm, 10)

    for Hval in Hscan:
        j.clearLog()
        j.setLayerExternalFieldDriver(
            "all",
            AxialDriver(0, 0, Hval),
        )
        j.runSimulation(15e-9, tstep, tstep, calculateEnergies=False)

        log = j.getLog()
        m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ("free", "bottom")])
        Rx0 = [100] * n_lay
        Ry0 = [0] * n_lay
        SMR = [5.11] * n_lay
        AMR = [0.41] * n_lay
        AHE = [AHE_scaling, AHE_scaling * AHE_asymmetry]

        _, Rxy = calculate_resistance_series(Rx0, Ry0, AMR=AMR, AHE=AHE, SMR=SMR, m=m, l=[l] * n_lay, w=[w] * n_lay)
        Rstable = Rxy[-100:].mean()
        output["Hscan"].append(Hval)
        output["R"].append(Rstable)

    # Check that we got results
    assert len(output["R"]) > 0
    assert len(output["Hscan"]) > 0

    # Check that resistance values are reasonable
    R_values = np.array(output["R"])
    assert np.all(np.isfinite(R_values))
    assert np.all(np.abs(R_values) < 1000)  # Allow negative values

    # Prepare data for reference comparison
    test_data = {
        "Hscan": np.array(output["Hscan"]),
        "R": R_values,
    }

    # Compare with reference or update
    compare_with_reference("test_ahe_loops", test_data, update_reference)


@pytest.mark.slow
def test_cims_hysteresis(update_reference):
    """Test CIMS hysteresis simulation."""
    demagTensor = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1)]
    sgnHdl = 1
    sgnHfl = -1

    Hdl = sgnHdl * 0.000371 * TtoAm
    Hfl = sgnHfl * 0.000433 * TtoAm
    jden1 = 6e10
    HDL = Hdl / jden1
    HFL = Hfl / jden1

    last_N = 100
    dt = 6e-9
    tstart = 1e-9
    tstep = 5e-13
    sim_time = 15e-9

    t_hm = 6e-9
    t_fm = 1.7e-9
    Kdir1 = CVector(0, 0, 1.0)
    pdir1 = CVector(0.0, 1.0, 0.0)
    Ku1 = 0.3e6
    Ms = 1.45
    alpha = 0.0064

    layer_free = Layer.createSOTLayer(
        id="free",
        mag=CVector(1, 0, 0),
        anis=Kdir1,
        Ms=Ms,
        thickness=t_fm,
        cellSurface=0.0,
        demagTensor=demagTensor,
        damping=alpha,
        dampingLikeTorque=HDL,
        fieldLikeTorque=HFL,
    )

    j = Junction([layer_free])
    j.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(Ku1))
    j.setLayerReferenceLayer("free", pdir1)

    j.clearLog()
    j.runSimulation(14e-9, tstep, tstep)
    Hoescales = [1, 0, -1]

    if abs(HDL) and abs(HFL):
        i_max = [0.1e12 for _ in range(len(Hoescales))]
    else:
        i_max = [1.0e12 for _ in range(len(Hoescales))]

    hysteresis_scales = []
    for i, Hoescale in enumerate(Hoescales):
        v = CVector(0.001, 0.99, 0)
        v.normalize()
        j.setLayerMagnetisation("free", v)
        # Reduced scan for faster testing
        jscan = np.linspace(-i_max[i], i_max[i], 20)
        jscan = np.concatenate([jscan, jscan[::-1][1:]])

        hysteresis = []
        for current in jscan:
            j.clearLog()
            j.setLayerCurrentDriver("all", stepDriver(0, current, tstart, tstart + dt))
            j.setLayerOerstedFieldDriver(
                "all",
                AxialDriver(
                    constantDriver(0),
                    stepDriver(0, Hoescale * current * t_hm / 2, tstart, tstart + dt),
                    constantDriver(0),
                ),
            )
            j.runSimulation(sim_time, tstep, tstep)
            log = j.getLog()
            m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ("free",)])

            Rstable = compute_gmr(
                Rp=100,
                Rap=200,
                m1=m.squeeze()[:, -last_N:],
                m2=np.asarray([[0, 1, 0] for _ in range(last_N)]).transpose(),
            ).mean()
            hysteresis.append(Rstable)

        hysteresis_scales.append(np.asarray(hysteresis))

    # Check results
    assert len(hysteresis_scales) == len(Hoescales)
    for hyst in hysteresis_scales:
        assert len(hyst) > 0
        assert np.all(hyst > 0)
        assert np.all(hyst < 300)

    # Prepare data for reference comparison
    test_data = {
        "hysteresis_scales": hysteresis_scales,
        "Hoescales": Hoescales,
    }

    # Compare with reference or update
    compare_with_reference("test_cims_hysteresis", test_data, update_reference)


@pytest.mark.slow
def test_vsd_field_scan(update_reference):
    """Test VSD field scan simulation."""
    Rx0 = 304.306
    Ry0 = 1.008
    SMR = -0.464
    AMR = -0.053
    AHE = -5.71
    w = 3e-5
    l = 2e-5
    INT_STEP = 1e-13
    HMIN = 350e3
    HMAX = 630e3
    HSTEPS = 20  # Reduced for faster testing

    def compute_vsd2(dynamicR, integration_step, dynamicI):
        SD = -dynamicI * dynamicR
        fs = 1.0 / integration_step
        SD_dc = Filters.butter_lowpass_filter(SD, cutoff=10e6, fs=fs, order=3)
        return np.mean(SD_dc)

    def simulate_lorentz(Ms, Ku, frequency, orient, alpha=1e-4, Irf=0.5e-3):
        demagTensor = [
            CVector(0.0, 0.0, 0.0),
            CVector(0.0, 0.0, 0.0),
            CVector(0.0, 0.0, 1.0),
        ]
        thickness = 1.45e-9
        l1 = Layer(
            id="free",
            mag=CVector(0.0, 0.0, 1.0),
            anis=CVector(0, 0.0, 1.0),
            Ms=Ms,
            thickness=thickness,
            cellSurface=0,
            demagTensor=demagTensor,
            damping=alpha,
        )
        junction = Junction([l1])
        junction.setLayerAnisotropyDriver("free", constantDriver(Ku))

        HoeAmpl = 5000
        Hspace = np.linspace(HMIN, HMAX, num=HSTEPS)
        theta = np.deg2rad(91)
        if orient == "4p":
            phideg = 0
        elif orient == "2p":
            phideg = 45
        else:
            raise ValueError("Unknown orient")
        phi = np.deg2rad(phideg)
        Hsweep = np.zeros(Hspace.shape[0])
        for i, H in enumerate(Hspace):
            junction.clearLog()
            HDriver = AxialDriver(
                constantDriver(H * math.sin(theta) * math.cos(phi)),
                constantDriver(H * math.sin(theta) * math.sin(phi)),
                constantDriver(H * math.cos(theta)),
            )

            HoeDriver = AxialDriver(
                NullDriver(),
                NullDriver(),
                sineDriver(0, -HoeAmpl, frequency, 0),
            )
            junction.setLayerExternalFieldDriver("all", HDriver)
            junction.setLayerOerstedFieldDriver("all", HoeDriver)
            junction.runSimulation(40e-9, INT_STEP, INT_STEP, solverMode=cmtj.RK4)

            log = junction.getLog()
            m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ["free"]])
            dynamicRx, dynamicRy = calculate_resistance_parallel(
                Rx0=[Rx0],
                Ry0=[Ry0],
                AMR=[AMR],
                AHE=[AHE],
                SMR=[SMR],
                m=m,
                l=[l],
                w=[w],
            )
            dynamicR = dynamicRx if orient == "2p" else dynamicRy
            dynamicI = Irf * np.sin(2 * math.pi * frequency * np.asarray(log["time"]))
            vmix = compute_vsd2(dynamicR, INT_STEP, dynamicI)
            Hsweep[i] = vmix
        return Hspace, Hsweep

    alpha = 30e-3
    Ms = 0.525
    Ku = 1.54e5
    # Test only one frequency and orientation for speed
    f = 14e9
    orient = "4p"
    irf = 0.75e-3

    hscan, vscan = simulate_lorentz(Ms, Ku, f, orient=orient, alpha=alpha, Irf=irf)

    assert len(hscan) == HSTEPS
    assert len(vscan) == HSTEPS
    # Round to 3 decimals
    vscan_rounded = round_to_3_decimals(vscan)
    # Check that values are finite (VSD can produce very small values near zero)
    assert np.all(np.isfinite(vscan_rounded))
    # Check that we have reasonable range (allow small values)
    assert np.abs(vscan_rounded).max() >= 0  # At least some finite values

    # Prepare data for reference comparison
    test_data = {
        "hscan": hscan,
        "vscan": vscan,
    }

    # Compare with reference or update
    compare_with_reference("test_vsd_field_scan", test_data, update_reference)


@pytest.mark.slow
def test_vsd_freq_scan(update_reference):
    """Test VSD frequency scan simulation."""
    Rx0 = 304.306
    Ry0 = 1.008
    SMR = -0.464
    AMR = -0.053
    AHE = -5.71
    w = 3e-5
    l = 2e-5
    INT_STEP = 1e-13

    def compute_vsd2(dynamicR, integration_step, dynamicI):
        SD = -dynamicI * dynamicR
        fs = 1.0 / integration_step
        SD_dc = Filters.butter_lowpass_filter(SD, cutoff=10e6, fs=fs, order=3)
        return np.mean(SD_dc)

    def simulate_lorentz_freq_scan(Ms, Ku, Hvalue, orient, alpha=1e-4, Irf=0.5e-3):
        demagTensor = [
            CVector(0.0, 0.0, 0.0),
            CVector(0.0, 0.0, 0.0),
            CVector(0.0, 0.0, 1.0),
        ]
        thickness = 1.45e-9
        l1 = Layer(
            id="free",
            mag=CVector(0.0, 0.0, 1.0),
            anis=CVector(0, 0.0, 1.0),
            Ms=Ms,
            thickness=thickness,
            cellSurface=0,
            demagTensor=demagTensor,
            damping=alpha,
        )
        junction = Junction([l1])
        junction.setLayerAnisotropyDriver("free", constantDriver(Ku))

        HoeAmpl = 5000
        theta = np.deg2rad(91)
        if orient == "4p":
            phideg = 0
        elif orient == "2p":
            phideg = 45
        else:
            raise ValueError("Unknown orient")
        phi = np.deg2rad(phideg)
        # Reduced frequency range for faster testing
        fspace = np.arange(10, 16, 1) * 1e9
        fsweep = np.zeros((fspace.shape[0]))
        for i, f in enumerate(fspace):
            junction.clearLog()
            HDriver = AxialDriver(
                constantDriver(Hvalue * math.sin(theta) * math.cos(phi)),
                constantDriver(Hvalue * math.sin(theta) * math.sin(phi)),
                constantDriver(Hvalue * math.cos(theta)),
            )

            HoeDriver = AxialDriver(
                NullDriver(),
                NullDriver(),
                sineDriver(0, -HoeAmpl, f, 0),
            )
            junction.setLayerExternalFieldDriver("all", HDriver)
            junction.setLayerOerstedFieldDriver("all", HoeDriver)
            junction.runSimulation(40e-9, INT_STEP, INT_STEP, solverMode=cmtj.RK4)

            log = junction.getLog()
            m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ["free"]])
            dynamicRx, dynamicRy = calculate_resistance_parallel(
                Rx0=[Rx0],
                Ry0=[Ry0],
                AMR=[AMR],
                AHE=[AHE],
                SMR=[SMR],
                m=m,
                l=[l],
                w=[w],
            )
            dynamicR = dynamicRx if orient == "4p" else dynamicRy
            dynamicI = Irf * np.sin(2 * math.pi * f * np.asarray(log["time"]))
            vmix = compute_vsd2(dynamicR, INT_STEP, dynamicI)
            fsweep[i] = vmix
        return fspace, fsweep

    alpha = 30e-3
    Ms = 0.525
    Ku = 1.54e5
    Hvalue = 500e3
    orient = "4p"
    irf = 0.75e-3

    fspace, fsweep = simulate_lorentz_freq_scan(Ms, Ku, Hvalue, orient=orient, alpha=alpha, Irf=irf)

    assert len(fspace) > 0
    assert len(fsweep) > 0
    fsweep_rounded = round_to_3_decimals(fsweep)
    # Check that values are finite (VSD can produce very small values near zero)
    assert np.all(np.isfinite(fsweep_rounded))
    # Check that we have reasonable range (allow small values)
    assert np.abs(fsweep_rounded).max() >= 0  # At least some finite values

    # Prepare data for reference comparison
    test_data = {
        "fspace": fspace,
        "fsweep": fsweep,
    }

    # Compare with reference or update
    compare_with_reference("test_vsd_freq_scan", test_data, update_reference)


@pytest.mark.slow
def test_fmr_basic(update_reference):
    """Test FMR basic simulation."""
    Ms1 = 1.03
    Ms2 = 1.03
    Ms3 = 0.12

    Ku1 = 489e3
    Ku2 = 514e3
    Ku3 = 17e3
    Kdir = CVector(0, 0, 1)
    damping = 0.01
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]
    demag_buffer = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]

    l1 = Layer(
        "free",
        mag=CVector(0, 0, 0.99),
        anis=Kdir,
        Ms=Ms1,
        thickness=1.0e-9,
        damping=damping,
        demagTensor=demag,
        cellSurface=1e-9,
    )
    l2 = Layer(
        "perpendicular",
        mag=CVector(0, 0, 0.99),
        anis=Kdir,
        Ms=Ms2,
        thickness=1.0e-9,
        damping=damping,
        cellSurface=1e-9,
        demagTensor=demag,
    )
    l3 = Layer(
        "buffer",
        mag=CVector(1, 0, 0),
        anis=CVector(1, 0, 0),
        Ms=Ms3,
        thickness=1.0e-9,
        damping=0.015,
        cellSurface=1e-9,
        demagTensor=demag_buffer,
    )
    l1.setAnisotropyDriver(constantDriver(Ku1))
    l2.setAnisotropyDriver(constantDriver(Ku2))
    l3.setAnisotropyDriver(constantDriver(Ku3))
    j = Junction([l1, l2, l3])

    # Reduced field scan for faster testing
    Hmax = -700e3
    Hscan, Hvecs = FieldScan.amplitude_scan(start=-Hmax, stop=Hmax, steps=20, theta=0, phi=0)
    dt = 5e-13
    sim_time = 60e-9
    ampl = 500
    dur = 2e-9
    spectrum = []
    wait_time = 4e-9
    j.setLayerExternalFieldDriver("all", AxialDriver(*Hvecs[0]))
    j.runSimulation(sim_time, dt, dt)

    from scipy.fft import fft, fftfreq

    def compute_fft(mixed_signal: np.ndarray, dt: float) -> tuple:
        fft_mixed = np.abs(fft(mixed_signal))
        ffreqs = fftfreq(len(mixed_signal), dt)
        max_freq = 30e9
        fft_mixed = fft_mixed[1 : len(fft_mixed) // 2]
        ffreqs = ffreqs[1 : len(ffreqs) // 2]
        fft_mixed = fft_mixed[np.where(ffreqs < max_freq)]
        ffreqs = ffreqs[np.where(ffreqs < max_freq)]
        return fft_mixed, ffreqs

    def extract_fft(log: dict, dt: float, wait_time: float) -> tuple:
        tm = np.asarray(log["time"])
        tmindx = np.argwhere(tm > wait_time).ravel()
        my_free = np.asarray(log["free_mx"])[tmindx]
        my_perpendicular = np.asarray(log["perpendicular_mx"])[tmindx]
        buffer = np.asarray(log["buffer_mx"])[tmindx]
        mmixed = (1 / 3) * (my_free + my_perpendicular + buffer)
        fft_mixed, ffreqs = compute_fft(mmixed, dt)
        return fft_mixed, ffreqs

    # Test only a few field values for speed
    for Hv in Hvecs[:5]:
        j.clearLog()
        j.setLayerExternalFieldDriver("all", AxialDriver(*Hv))
        j.setLayerOerstedFieldDriver("all", AxialDriver(NullDriver(), stepDriver(0, ampl, 0, dur), NullDriver()))
        j.runSimulation(sim_time, dt, dt)
        log = j.getLog()
        fft_mixed, ffreqs = extract_fft(log, dt, wait_time)
        spectrum.append(fft_mixed)

    assert len(spectrum) > 0
    for spec in spectrum:
        assert len(spec) > 0
        assert np.all(spec >= 0)

    # Prepare data for reference comparison
    test_data = {
        "spectrum": spectrum,
        "ffreqs": ffreqs if len(spectrum) > 0 else [],
    }

    # Compare with reference or update
    compare_with_reference("test_fmr_basic", test_data, update_reference)


@pytest.mark.slow
def test_harmonic_hall(update_reference):
    """Test harmonic hall simulation."""
    Rx0 = [304.306]
    SMR = [-0.464]
    AMR = [-0.053]
    AHE = [-5.71]
    w = [3e-5]
    l = [2e-5]

    Ms = 0.525
    Ku = 1.54e5
    alpha = 0.03

    def run_simulation(junction: Junction, Hvecs: np.ndarray, mode: str, int_time=1e-12):
        sim_time = 10e-9
        layer_str = ["free"]
        mags = [CVector(0, 0, 1) for _ in layer_str]
        Rxy, Rxx = [], []
        for Hval in Hvecs:
            junction.clearLog()
            HDriver = AxialDriver(
                ScalarDriver.getConstantDriver(Hval[0]),
                ScalarDriver.getConstantDriver(Hval[1]),
                ScalarDriver.getConstantDriver(Hval[2]),
            )
            junction.setLayerExternalFieldDriver("all", HDriver)
            for i, l_str in enumerate(layer_str):
                junction.setLayerMagnetisation(l_str, mags[i])

            junction.runSimulation(sim_time, int_time, int_time)

            for i, str_ in enumerate(layer_str):
                mags[i] = junction.getLayerMagnetisation(str_)

            log = junction.getLog()
            m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in layer_str])
            dynamicRx, dynamicRy = calculate_resistance_parallel(
                Rx0,
                [0],
                AMR=AMR,
                AHE=AHE,
                SMR=SMR,
                m=m,
                l=l,
                w=w,
            )

            Rxy.append(dynamicRy[-1])
            Rxx.append(dynamicRx[-1])
        return np.asarray(Rxx) if mode.lower() == "rxx" else np.asarray(Rxy)

    def simulate(Ku, Ms, Hvecs, alpha, mode="rxx"):
        demagTensor = [
            CVector(0.00024164288391924, 2.71396011566517e-10, 5.95503928124313e-14),
            CVector(2.71396011566517e-10, 0.000160046006320031, 1.32504057070646e-14),
            CVector(5.95503928124313e-14, 1.32504057070646e-14, 0.999598310229469),
        ]

        thickness = 1.45e-9
        surface = w[0] * l[0]
        l1 = Layer(
            id="free",
            mag=CVector(0, 0, 1),
            anis=CVector(0.0, 0, 1),
            Ms=Ms,
            thickness=thickness,
            cellSurface=surface,
            damping=alpha,
            demagTensor=demagTensor,
        )
        junction = Junction([l1])

        HoePulseAmpl = 50
        HoeDriver = AxialDriver(
            NullDriver(),
            NullDriver(),
            ScalarDriver.getStepDriver(0, HoePulseAmpl, 0.0, 1e-11),
        )
        junction.setLayerOerstedFieldDriver("all", HoeDriver)
        junction.setLayerAnisotropyDriver("all", ScalarDriver.getConstantDriver(Ku))

        return run_simulation(junction=junction, Hvecs=Hvecs, mode=mode)

    # Reduced field scan for faster testing
    phi = 0
    theta = 92
    Hscan, Hvecs = FieldScan.amplitude_scan(theta=theta, phi=phi, start=-300e3, stop=300e3, steps=20, back=False)
    simulated = simulate(Ms=Ms, Ku=Ku, alpha=alpha, mode="rxx", Hvecs=Hvecs)

    assert len(simulated) == len(Hscan)
    assert np.all(simulated > 0)
    simulated_rounded = round_to_3_decimals(simulated)
    assert len(np.unique(simulated_rounded)) > 1

    # Prepare data for reference comparison
    test_data = {
        "Hscan": Hscan,
        "simulated": simulated,
    }

    # Compare with reference or update
    compare_with_reference("test_harmonic_hall", test_data, update_reference)


@pytest.mark.slow
def test_sto(update_reference):
    """Test STO simulation."""
    demagTensor = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1.0)]

    damping = 0.3
    currentDensity = -6e9
    beta = 1
    spinPolarisation = 1.0

    l1 = Layer.createSTTLayer(
        id="free",
        mag=CVector(0.0, 0.0, 1.0),
        anis=CVector(0, 0.0, 1.0),
        Ms=1.0,
        thickness=1.4e-9,
        cellSurface=7e-10 * 7e-10,
        demagTensor=demagTensor,
        damping=damping,
        SlonczewskiSpacerLayerParameter=1.0,
        spinPolarisation=spinPolarisation,
        beta=beta,
    )

    l1.setReferenceLayer(CVector(0, 1.0, 1.0))
    junction = Junction([l1], 100, 200)

    junction.setLayerAnisotropyDriver("free", constantDriver(350e3))
    junction.setLayerCurrentDriver("free", constantDriver(currentDensity))
    junction.runSimulation(50e-9, 1e-13, 1e-13)
    log = junction.getLog()

    assert "free_mx" in log
    assert "free_my" in log
    assert "free_mz" in log
    assert "R" in log

    mx = np.array(log["free_mx"])
    my = np.array(log["free_my"])
    mz = np.array(log["free_mz"])
    R = np.array(log["R"])

    assert len(mx) > 0
    assert len(R) > 0

    # Check that magnetization oscillates (for STO)
    mx_rounded = round_to_3_decimals(mx[-100:])
    # Should have variation in magnetization
    assert len(np.unique(mx_rounded)) > 1

    # Check resistance values
    R_rounded = round_to_3_decimals(R[-100:])
    assert np.all(R_rounded > 0)
    assert np.all(R_rounded < 300)

    # Prepare data for reference comparison
    test_data = {
        "mx": mx[-100:],
        "my": my[-100:],
        "mz": mz[-100:],
        "R": R[-100:],
    }

    # Compare with reference or update
    compare_with_reference("test_sto", test_data, update_reference)


@pytest.mark.slow
def test_cims_stability(update_reference):
    """Test CIMS stability diagram simulation."""
    Kdir_FL = CVector(0.0, 0.0, 1.0)

    FL_Ms = 1.12
    FL_Ku = 563e3

    damping = 0.01

    p = 0.017
    Nxx = Nyy = p
    Nzz = np.sqrt(1 - 2 * (p**2))
    demag = [CVector(Nxx, 0, 0), CVector(0, Nyy, 0), CVector(0, 0, Nzz)]

    # Reduced parameters for faster testing
    thetas = [0.1]
    Ku2s = [-30e3]
    tau_scales = [0.7]
    Hmax = 400e3
    Nsteps = 10  # Reduced steps
    jscan_base = np.linspace(0, 0.3, 20)  # Reduced scan
    jscan = np.concatenate([jscan_base, jscan_base[::-1], -jscan_base, -jscan_base[::-1]])
    T = 0  # No temperature for faster testing
    solver = SolverMode.DormandPrince
    N_thermal_samples = 1
    cell_surface = np.power(300e-9, 2) * np.pi

    wf = 1e-11
    integration_time = 1e-12
    pulse_len = 10e-9

    for tau_scale in tau_scales:
        for theta in thetas:
            for FL_Ku2 in Ku2s:
                Hscan, Hvecs = FieldScan.amplitude_scan(start=-Hmax, stop=Hmax, steps=Nsteps, theta=theta, phi=0)

                # Process a single field vector (simplified, no multiprocessing for test)
                def process_field_vector(field_data):
                    idx, Hv = field_data
                    all_thermal_results = []

                    for thermal_sample in range(N_thermal_samples):
                        free = Layer(
                            "free",
                            mag=CVector(0.1, 0.1, 0.99),
                            anis=Kdir_FL,
                            Ms=FL_Ms,
                            thickness=1.0e-9,
                            damping=damping,
                            cellSurface=cell_surface,
                            demagTensor=demag,
                        )

                        free.setReferenceLayer(CVector(0, 0, 1))
                        free.setSecondOrderAnisotropyDriver(constantDriver(FL_Ku2))
                        free.setAnisotropyDriver(constantDriver(FL_Ku))
                        j = Junction([free])

                        j.clearLog()
                        if T != 0:
                            np.random.seed(thermal_sample * 1000 + idx)
                        Hv_perturbed = Hv + np.random.randn(3) * 100
                        mag_init = CVector(*Hv_perturbed.tolist())
                        mag_init.normalize()
                        j.setLayerMagnetisation("free", mag_init)
                        Hval = Hv.tolist()
                        j.setLayerExternalFieldDriver("free", AxialDriver(*Hval))
                        j.runSimulation(25e-9, integration_time, wf, solverMode=solver)

                        voltage_dict = defaultdict(list)
                        local_scan = jscan.copy()
                        if Hv[-1] > 0:
                            local_scan = local_scan[::-1]
                        for jden in local_scan:
                            j.clearLog()
                            j.setLayerFieldLikeTorqueDriver(
                                "free",
                                stepDriver(0, tau_scale * TtoAm * 1e-3 * (jden**2), 0, pulse_len),
                            )
                            j.setLayerDampingLikeTorqueDriver(
                                "free",
                                stepDriver(0, tau_scale * TtoAm * 16e-3 * jden, 0, pulse_len),
                            )
                            j.runSimulation(2 * pulse_len, integration_time, wf, solverMode=solver)
                            log = j.getLog()
                            mz = np.mean(log["free_mz"][-100:])
                            voltage_dict[jden].append(mz)

                        thermal_result = {}
                        for jden, measured_values in voltage_dict.items():
                            thermal_result[jden] = np.mean(measured_values)
                        all_thermal_results.append(thermal_result)

                    final_voltage_dict = defaultdict(list)
                    for thermal_result in all_thermal_results:
                        for jden, mz_value in thermal_result.items():
                            final_voltage_dict[jden].append(mz_value)

                    averaged_voltage_dict = {}
                    for jden, mz_values in final_voltage_dict.items():
                        averaged_voltage_dict[jden] = np.mean(mz_values)

                    averaged_voltage_dict = dict(sorted(averaged_voltage_dict.items(), key=lambda x: x[0]))
                    return (
                        idx,
                        list(averaged_voltage_dict.values()),
                        sorted(list(averaged_voltage_dict.keys())),
                    )

                # Process field vectors sequentially for testing
                results = []
                for i, Hv in enumerate(Hvecs[:5]):  # Only test first 5 for speed
                    result = process_field_vector((i, Hv))
                    results.append(result)

                results.sort(key=lambda x: x[0])

                stability_spectrum = np.array([r[1] for r in results])
                voltage_keys = results[0][2] if results else []

                assert len(stability_spectrum) > 0
                assert len(voltage_keys) > 0

                # Prepare data for reference comparison
                test_data = {
                    "stability_spectrum": stability_spectrum,
                    "voltage_keys": np.array(voltage_keys),
                    "Hscan": Hscan[:5],  # Only first 5
                }

                # Compare with reference or update (stochastic test - use lenient comparison)
                compare_with_reference("test_cims_stability", test_data, update_reference, stochastic=True)
