import pytest
from cmtj import *
from cmtj.utils.procedures import ResistanceParameters
from typing import Tuple


@pytest.mark.parametrize(
    "arg", ["single_layer_mtj", "two_layer_mtj", "tri_layer_mtj"], indirect=True
)
def test_base_solver(arg: Tuple[Junction, ResistanceParameters]):
    junction, _ = arg
    junction.runSimulation(10e-9, 1e-12, 1e-12)
    junction.runSimulation(10e-9, 1e-12, 1e-12, solverMode=SolverMode.EulerHeun)
    junction.runSimulation(10e-9, 1e-12, 1e-12, solverMode=SolverMode.Heun)


@pytest.mark.parametrize(
    "arg", ["single_layer_mtj", "two_layer_mtj", "tri_layer_mtj"], indirect=True
)
def test_base_solver(arg: Tuple[Junction, ResistanceParameters]):
    junction, _ = arg
    junction.setLayerTemperatureDriver("all", ScalarDriver.getConstantDriver(300))
    junction.runSimulation(10e-9, 1e-12, 1e-12)
    junction.runSimulation(10e-9, 1e-12, 1e-12, solverMode=SolverMode.EulerHeun)
    junction.runSimulation(10e-9, 1e-12, 1e-12, solverMode=SolverMode.Heun)
