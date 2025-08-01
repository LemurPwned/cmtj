from cmtj.stack import ParallelStack, SeriesStack
from cmtj.utils.procedures import ResistanceParameters
from cmtj import Layer, CVector, ScalarDriver
from typing import Tuple
from cmtj import Junction
import pytest
import numpy as np
from math import pi


def test_invalid_stack_indices(
    single_layer_mtj_fictious: Tuple[Junction, ResistanceParameters],
):
    junction, _ = single_layer_mtj_fictious
    with pytest.raises(RuntimeError, match="Asking for id of a non-existing junction!"):
        ParallelStack([junction, junction]).getLog(2)

    with pytest.raises(RuntimeError, match="Asking for id of a non-existing junction!"):
        SeriesStack([junction, junction]).getLog(10)

    with pytest.raises(TypeError):
        ParallelStack([junction, junction]).getLog(-1)

    with pytest.raises(TypeError):
        SeriesStack([junction, junction]).getLog(-10)


def test_invalid_stack(
    single_layer_mtj_fictious: Tuple[Junction, ResistanceParameters],
):
    junction, _ = single_layer_mtj_fictious
    with pytest.raises(RuntimeError, match="must have at least 2 junctions"):
        ParallelStack([junction])
    with pytest.raises(RuntimeError, match="must have at least 2 junctions"):
        SeriesStack([junction])
    with pytest.raises(RuntimeError, match="must have at least 2 junctions"):
        ParallelStack([])


@pytest.mark.parametrize(
    "arg", ["single_layer_mtj_fictious", "two_layer_mtj"], indirect=True
)
def test_basic_parallel_stack(arg: Tuple[Junction, ResistanceParameters]):
    junction, _ = arg
    stack = SeriesStack([junction, junction])
    for c in (0, 0.1, 0.5, 1):
        stack.setCouplingStrength(c)
        stack.runSimulation(5e-9, 1e-12, 1e-12)
        log = stack.getLog()
        assert "Resistance" in log.keys()


@pytest.mark.parametrize(
    "arg", ["single_layer_mtj_fictious", "two_layer_mtj"], indirect=True
)
def test_basic_series_stack(arg: Tuple[Junction, ResistanceParameters]):
    junction, _ = arg
    stack = ParallelStack([junction, junction])
    for c in (0, 0.1, 0.5, 1):
        stack.setCouplingStrength(c)
        stack.runSimulation(5e-9, 1e-12, 1e-12)
        log = stack.getLog()
        assert "Resistance" in log.keys()


def test_stack_simulation_parameters():
    """Test stack behavior with invalid simulation parameters"""
    junction = Junction(
        [
            Layer(
                "free",
                CVector(0, 0, 1),
                CVector(0, 0, 1),
                Ms=1e6,
                thickness=2e-9,
                cellSurface=np.pi * (20e-9) ** 2,
                demagTensor=[CVector(0, 0, 0)] * 3,
            )
        ],
    )
    with pytest.raises(RuntimeError, match="must have at least 2 junctions"):
        ParallelStack([junction])
    with pytest.raises(RuntimeError, match="must have at least 2 junctions"):
        SeriesStack([junction])


@pytest.mark.parametrize(
    "arg", ["single_layer_mtj_fictious", "two_layer_mtj"], indirect=True
)
def test_stack_log_consistency(arg: Tuple[Junction, ResistanceParameters]):
    """Test that stack logs maintain consistency across simulations"""
    junction, _ = arg
    stack = ParallelStack([junction, junction])

    # Run first simulation
    stack.setCouplingStrength(0.5)
    stack.runSimulation(5e-9, 1e-12, 1e-12)
    log1 = stack.getLog()

    # Clear and run second simulation
    stack.clearLogs()
    stack.runSimulation(5e-9, 1e-12, 1e-12)
    log2 = stack.getLog()

    # Verify log structure remains consistent
    assert set(log1.keys()) == set(log2.keys())
    assert len(log1["time"]) == len(log2["time"])
    assert "Resistance" in log1 and "Resistance" in log2


@pytest.mark.parametrize(
    "arg", ["single_layer_mtj_fictious", "two_layer_mtj"], indirect=True
)
def test_stack_resistance_behavior(arg: Tuple[Junction, ResistanceParameters]):
    """Test that resistance values follow expected patterns"""
    junction, params = arg
    parallel_stack = ParallelStack([junction, junction])
    series_stack = SeriesStack([junction, junction])

    # Test with no coupling
    parallel_stack.setCouplingStrength(0)
    series_stack.setCouplingStrength(0)

    parallel_stack.runSimulation(5e-9, 1e-12, 1e-12)
    series_stack.runSimulation(5e-9, 1e-12, 1e-12)

    p_log = parallel_stack.getLog()
    s_log = series_stack.getLog()

    # Basic sanity checks for resistance values
    assert all(r > 0 for r in p_log["Resistance"])  # Resistance should be positive
    assert all(r > 0 for r in s_log["Resistance"])

    # Series resistance should be larger than parallel
    assert np.mean(s_log["Resistance"]) > np.mean(p_log["Resistance"])


@pytest.mark.parametrize("stack_type", [ParallelStack, SeriesStack])
def test_kcl_vs_non_kcl_current_behavior(stack_type):
    """Test that KCL and non-KCL modes behave differently for current distribution"""
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    alpha = 0.005
    surface = 150e-9 * 150e-9 * pi
    Kdir = CVector(1, 0, 0)
    l1 = Layer(
        "free",
        mag=CVector(0, 0, 1),
        anis=Kdir,
        Ms=1.65,
        thickness=3e-9,
        cellSurface=surface,
        demagTensor=demag,
        damping=alpha,
    )
    K1 = 1.05e3
    l1.setReferenceLayer(CVector(1, 0, 0))
    junction = Junction([l1], 100, 200)
    junction.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K1))

    Kdir = CVector(1, 0.1, 0.3)
    l1 = Layer(
        "free",
        mag=CVector(0.3, 0.2, 0.9),
        anis=Kdir,
        Ms=1.35,
        thickness=3e-9,
        cellSurface=surface,
        demagTensor=demag,
        damping=alpha,
    )
    l1.setReferenceLayer(CVector(1, 0, 0))
    junction2 = Junction([l1], 100, 200)
    junction2.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K1 * 0.1))

    test_current = 2e12
    current_driver = ScalarDriver.getConstantDriver(test_current)

    stack_non_kcl = stack_type([junction, junction2], useKCL=False)
    stack_kcl = stack_type([junction, junction2], useKCL=True)

    stack_non_kcl.setCouplingStrength(0.3)
    stack_kcl.setCouplingStrength(0.3)

    # Set the current driver for both stacks
    stack_non_kcl.setCoupledCurrentDriver(current_driver)
    stack_kcl.setCoupledCurrentDriver(current_driver)

    # Run short simulations
    simulation_time = 0.1e-9  # 1 ns
    time_step = 1e-12  # 1 ps
    write_freq = 1e-12  # 1 ps

    stack_non_kcl.runSimulation(simulation_time, time_step, write_freq)
    stack_kcl.runSimulation(simulation_time, time_step, write_freq)

    # Get the logs for both stacks
    log_non_kcl = stack_non_kcl.getLog()
    log_kcl = stack_kcl.getLog()

    # For non-KCL: current in first junction should be the same as set current
    current_first_junction_non_kcl = np.array(log_non_kcl["I_0"]) / 1e12

    # For KCL: current in first junction should be different from set current
    current_first_junction_kcl = np.array(log_kcl["I_0"]) / 1e12

    mean_I0_non_kcl = np.mean(current_first_junction_non_kcl)
    mean_I0_kcl = np.mean(current_first_junction_kcl)
    test_current /= 1e12

    # Verify that non-KCL preserves the input current in the first junction
    # (allowing for small numerical errors)
    assert np.allclose(
        current_first_junction_non_kcl, test_current, rtol=1e-6
    ), f"Non-KCL mode: Expected current {test_current}, got {current_first_junction_non_kcl}"

    # Verify that KCL mode changes the current distribution
    # The current should NOT be the same as the input current
    assert not np.allclose(
        mean_I0_kcl, test_current, rtol=1e-6
    ), f"KCL mode: Current should be different from input {test_current}, but got {mean_I0_kcl}"

    # Verify that the two modes produce different current distributions
    assert not np.allclose(
        mean_I0_non_kcl, mean_I0_kcl, rtol=1e-6
    ), "KCL and non-KCL modes should produce different current distributions"

    # Additional verification: Both stacks should have the same number of logged data points
    assert len(log_non_kcl["time"]) == len(log_kcl["time"])
    assert (
        "I_1" in log_non_kcl and "I_1" in log_kcl
    )  # Both should log current for second junction
