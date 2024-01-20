from cmtj.stack import ParallelStack, SeriesStack
from cmtj.utils.procedures import ResistanceParameters
from typing import Tuple
from cmtj import Junction
import pytest


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
