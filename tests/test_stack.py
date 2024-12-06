from cmtj.stack import ParallelStack, SeriesStack
from cmtj.utils.procedures import ResistanceParameters
from cmtj import Layer, CVector
from typing import Tuple
from cmtj import Junction
import pytest
import numpy as np


def test_invalid_stack_indices(
    single_layer_mtj_fictious: Tuple[Junction, ResistanceParameters]
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
    single_layer_mtj_fictious: Tuple[Junction, ResistanceParameters]
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
