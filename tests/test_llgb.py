from cmtj.llgb import LLGBLayer, LLGBJunction
from cmtj import CVector, ScalarDriver
import numpy as np
import pytest


def test_basic():
    Ms = 0.27
    Tc = 448
    susceptibility = 0.04
    me = 0.9
    damping = 0.0275
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    l1 = LLGBLayer(
        "free",
        CVector(1, 0, 0),
        CVector(1, 0, 0),
        Ms,
        2e-9,
        100e-9 * 100e-9,
        demag,
        damping,
        Tc,
        susceptibility,
        me,
    )
    junction = LLGBJunction([l1])
    sim_time = 10e-9
    for i in range(18):
        junction.setLayerTemperatureDriver(
            "all", ScalarDriver.getConstantDriver(40 * i)
        )
        junction.runSimulation(sim_time, 1e-13, 1e-13)
    log = junction.getLog()
    assert "free_T" in log


def test_temperature_dependence():
    """Test magnetization response to temperature changes"""
    Ms = 0.27
    Tc = 448  # Curie temperature
    susceptibility = 0.04
    me = 0.9
    damping = 0.0275
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]

    layer = LLGBLayer(
        "free",
        CVector(1, 0, 0),
        CVector(1, 0, 0),
        Ms,
        2e-9,
        100e-9 * 100e-9,
        demag,
        damping,
        Tc,
        susceptibility,
        me,
    )
    junction = LLGBJunction([layer])

    # Test at room temperature (300K)
    junction.setLayerTemperatureDriver("all", ScalarDriver.getConstantDriver(300))
    junction.runSimulation(1e-9, 1e-13, 1e-13)
    log_room = junction.getLog()
    mx_room = log_room["free_mx"][-1]

    # Test near Curie temperature
    junction.setLayerTemperatureDriver("all", ScalarDriver.getConstantDriver(Tc - 1))
    junction.runSimulation(1e-9, 1e-13, 1e-13)
    log_hot = junction.getLog()
    mx_hot = log_hot["free_mx"][-1]

    # Magnetization should be significantly reduced near Tc
    assert abs(mx_hot) < abs(mx_room)


def test_multiple_layers():
    """Test LLGB with multiple layers"""
    Ms = 0.27
    Tc = 448
    susceptibility = 0.04
    me = 0.9
    damping = 0.0275
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]

    layer1 = LLGBLayer(
        "free1",
        CVector(1, 0, 0),
        CVector(1, 0, 0),
        Ms,
        2e-9,
        100e-9 * 100e-9,
        demag,
        damping,
        Tc,
        susceptibility,
        me,
    )

    layer2 = LLGBLayer(
        "free2",
        CVector(-1, 0, 0),  # Opposite initial magnetization
        CVector(1, 0, 0),
        Ms,
        2e-9,
        100e-9 * 100e-9,
        demag,
        damping,
        Tc,
        susceptibility,
        me,
    )

    junction = LLGBJunction([layer1, layer2])
    junction.setLayerTemperatureDriver("all", ScalarDriver.getConstantDriver(300))
    junction.runSimulation(5e-9, 1e-13, 1e-13)

    log = junction.getLog()
    assert "free1_mx" in log
    assert "free2_mx" in log
    # Layers should maintain opposite magnetization
    assert log["free1_mx"][-1] * log["free2_mx"][-1] < 0


def _make_llgb_layer(layer_id="free", seed=None):
    """Helper: build a simple thermally-driven LLGBLayer with an optional seed."""
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    layer = LLGBLayer(
        layer_id,
        CVector(1, 0, 0),
        CVector(1, 0, 0),
        Ms=0.27,
        thickness=2e-9,
        cellSurface=100e-9 * 100e-9,
        demagTensor=demag,
        damping=0.0275,
        Tc=448,
        susceptibility=0.04,
        me=0.9,
    )
    if seed is not None:
        layer.setSeed(seed)
    return layer


def test_llgblayer_set_seed_reproducibility():
    """Same seed must produce identical stochastic LLGBLayer trajectories."""
    sim_time = 1e-9
    dt = 1e-13

    def run_with_seed(seed):
        layer = _make_llgb_layer(seed=seed)
        junction = LLGBJunction([layer])
        junction.setLayerTemperatureDriver("all", ScalarDriver.getConstantDriver(300))
        junction.runSimulation(sim_time, dt, dt)
        return np.asarray(junction.getLog()["free_mx"])

    traj_a = run_with_seed(7)
    traj_b = run_with_seed(7)
    np.testing.assert_array_equal(
        traj_a, traj_b, err_msg="LLGB trajectories with the same seed must be identical"
    )


@pytest.mark.skip(
    reason="LLGBLayer stochastic thermal fields are currently zeroed out in "
           "the C++ implementation (work in progress); enable once stochastic "
           "terms are activated."
)
def test_llgblayer_different_seeds_differ():
    """Different seeds must produce different stochastic LLGBLayer trajectories."""
    sim_time = 1e-9
    dt = 1e-13

    def run_with_seed(seed):
        layer = _make_llgb_layer(seed=seed)
        junction = LLGBJunction([layer])
        junction.setLayerTemperatureDriver("all", ScalarDriver.getConstantDriver(300))
        junction.runSimulation(sim_time, dt, dt)
        return np.asarray(junction.getLog()["free_mx"])

    traj_a = run_with_seed(3)
    traj_b = run_with_seed(4)
    assert not np.array_equal(
        traj_a, traj_b
    ), "LLGB trajectories with different seeds should differ"
