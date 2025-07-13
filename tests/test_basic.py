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
