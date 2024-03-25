from typing import Tuple

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
