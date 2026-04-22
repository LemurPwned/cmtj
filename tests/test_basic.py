import pytest
import numpy as np
from cmtj import Junction, AxialDriver, CVector, Layer, ScalarDriver


def test_no_p_mtj():
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    alpha = 0.005
    Kdir = CVector(1, 0, 0)
    l1 = Layer(
        "free",
        mag=CVector(0, 0, 1),
        anis=Kdir,
        Ms=1.65,
        thickness=3e-9,
        cellSurface=0,
        demagTensor=demag,
        damping=alpha,
    )
    with pytest.raises(
        ValueError,
        match="must have a pinning",
    ):
        j1 = Junction([l1], 100, 200)


def test_basic_mtj():
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    alpha = 0.005
    Kdir = CVector(1, 0, 0)
    l1 = Layer(
        "free",
        mag=CVector(0, 0, 1),
        anis=Kdir,
        Ms=1.65,
        thickness=3e-9,
        cellSurface=0,
        demagTensor=demag,
        damping=alpha,
    )
    K1 = 1.05e3
    H = CVector(300e3, 10e3, 0)
    l1.setReferenceLayer(CVector(1, 0, 0))
    junction = Junction([l1], 100, 200)

    junction.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K1))
    junction.setLayerExternalFieldDriver("all", AxialDriver(H))
    junction.runSimulation(10e-9, 1e-12, 1e-12)
    log = junction.getLog()
    assert "R" in log.keys()
    R = np.asarray(log["R"])
    assert np.all(R)


def test_second_order_anisotropy():
    # Create a test layer with demagnetization tensor
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]

    alpha = 0.03
    Ms = 0.5
    Kdir = CVector(0, 0, 1)  # Anisotropy axis along z
    l1 = Layer(
        "free",
        mag=CVector(0, 0, 1),  # Initial magnetization along z
        anis=Kdir,
        Ms=Ms,  # Saturation magnetization
        thickness=1e-9,
        cellSurface=70e-9,
        demagTensor=demag,
        damping=alpha,
    )

    # Set up second order anisotropy driver
    K2 = 1e6  # Second order anisotropy constant
    l1.setAnisotropyDriver(ScalarDriver.getConstantDriver(K2))
    l1.setSecondOrderAnisotropyDriver(ScalarDriver.getConstantDriver(K2))
    j = Junction([l1])
    j.runSimulation(1e-9, 1e-12, 1e-12, calculateEnergies=True)
    log = j.getLog()
    assert "free_K2" in log.keys()
    K2_log = np.asarray(log["free_K2"])
    assert np.all(K2_log == K2)

    l1 = Layer(
        "free",
        mag=CVector(0, 0, 1),  # Initial magnetization along z
        anis=Kdir,
        Ms=Ms,  # Saturation magnetization
        thickness=1e-9,
        cellSurface=70e-9,
        demagTensor=demag,
        damping=alpha,
    )

    l1.setAnisotropyDriver(ScalarDriver.getConstantDriver(K2))
    j = Junction([l1])
    j.setLayerSecondOrderAnisotropyDriver("free", ScalarDriver.getConstantDriver(K2))
    j.runSimulation(1e-9, 1e-12, 1e-12, calculateEnergies=True)
    log = j.getLog()
    assert "free_K2" in log.keys()
    K2_log = np.asarray(log["free_K2"])
    assert np.all(K2_log == K2)


def _make_thermal_layer(seed=None):
    """Helper: build a simple thermally-driven Layer with an optional seed."""
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    l = Layer(
        "free",
        mag=CVector(0, 0, 1),
        anis=CVector(0, 0, 1),
        Ms=1.0,
        thickness=2e-9,
        cellSurface=150e-9 * 150e-9,
        demagTensor=demag,
        damping=0.01,
    )
    if seed is not None:
        l.setSeed(seed)
    l.setTemperatureDriver(ScalarDriver.getConstantDriver(300))
    return l


def test_layer_set_seed_reproducibility():
    """Same seed must produce identical stochastic trajectories."""
    K = 1e3
    sim_time = 2e-9
    dt = 1e-12

    def run_with_seed(seed):
        l = _make_thermal_layer(seed=seed)
        j = Junction([l])
        j.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K))
        j.runSimulation(sim_time, dt, dt)
        return np.asarray(j.getLog()["free_mz"])

    traj_a = run_with_seed(42)
    traj_b = run_with_seed(42)
    np.testing.assert_array_equal(
        traj_a, traj_b, err_msg="Trajectories with the same seed must be identical"
    )


def test_layer_different_seeds_differ():
    """Different seeds must produce different stochastic trajectories."""
    K = 1e3
    sim_time = 2e-9
    dt = 1e-12

    def run_with_seed(seed):
        l = _make_thermal_layer(seed=seed)
        j = Junction([l])
        j.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K))
        j.runSimulation(sim_time, dt, dt)
        return np.asarray(j.getLog()["free_mz"])

    traj_a = run_with_seed(1)
    traj_b = run_with_seed(2)
    assert not np.array_equal(
        traj_a, traj_b
    ), "Trajectories with different seeds should differ"
