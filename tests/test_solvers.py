import numpy as np
import pytest
from cmtj import *
from cmtj import constants
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


def _precessing_junction(n_layers: int, H: float):
    demag = [CVector(0, 0, 0)] * 3
    layers = [
        Layer(f"L{i}", CVector(1, 0, 0), CVector(0, 0, 1), 1.0, 1e-9, 1e-16, demag, 0.0)
        for i in range(n_layers)
    ]
    junction = Junction(layers)
    junction.setLayerExternalFieldDriver("all", AxialDriver(0, 0, H))
    return junction


def test_dormand_prince_time_axis():
    """Logged time must match the precession phase of a free spin."""
    H = 100e3
    total_time = 1e-9
    junction = _precessing_junction(1, H)
    junction.runSimulation(total_time, 1e-13, 1e-12, solverMode=SolverMode.DormandPrince)
    log = junction.getLog()
    t = np.asarray(log["time"])
    phase = np.unwrap(np.arctan2(log["L0_my"], log["L0_mx"]))
    omega = np.abs(phase[-1] - phase[0]) / (t[-1] - t[0])
    expected = constants.PhysicalConstants.gyromagnetic_ratio() * H
    assert omega == pytest.approx(expected, rel=1e-3)
    assert t[-1] <= total_time * (1 + 1e-9)


def test_dormand_prince_multilayer():
    total_time = 2e-10
    junction = _precessing_junction(2, 100e3)
    junction.setIECDriver("L0", "L1", ScalarDriver.getConstantDriver(1e-3))
    junction.runSimulation(total_time, 1e-13, 1e-12, solverMode=SolverMode.DormandPrince)
    t = np.asarray(junction.getLog()["time"])
    assert np.all(np.diff(t) > 0)
    assert t[-1] == pytest.approx(total_time, rel=1e-6)


def test_dormand_prince_multilayer_rejection_rollback():
    """A layer rejecting steps must not let its neighbour run ahead of time.

    L1 has a tight tolerance and a pulsed field, so it rejects steps that the
    freely precessing L0 accepts. L0's phase must still follow the logged time.
    """
    H = 100e3
    demag = [CVector(0, 0, 0)] * 3
    free = Layer("L0", CVector(1, 0, 0), CVector(0, 0, 1), 1.0, 1e-9, 1e-16, demag, 0.0)
    pulsed = Layer("L1", CVector(1, 0, 0), CVector(0, 0, 1), 1.0, 1e-9, 1e-16, demag, 0.01)
    params = AdaptiveIntegrationParams()
    params.abs_tol = 1e-9
    params.rel_tol = 1e-9
    pulsed.setAdaptiveParams(params)
    junction = Junction([free, pulsed])
    junction.setLayerExternalFieldDriver("L0", AxialDriver(0, 0, H))
    junction.setLayerExternalFieldDriver(
        "L1",
        AxialDriver(NullDriver(), NullDriver(), ScalarDriver.getPulseDriver(0, 300e3, 2e-11, 0.5)),
    )
    junction.runSimulation(5e-10, 1e-13, 1e-12, solverMode=SolverMode.DormandPrince)
    log = junction.getLog()
    t = np.asarray(log["time"])
    phase = np.abs(np.unwrap(np.arctan2(log["L0_my"], log["L0_mx"])))
    omega = constants.PhysicalConstants.gyromagnetic_ratio() * H
    assert np.abs(phase - omega * t).max() < 0.05
