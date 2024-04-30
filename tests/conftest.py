import pytest
from typing import Tuple
from math import pi
from cmtj.utils import mu0
from cmtj import Junction, CVector, Layer, ScalarDriver
from cmtj.utils.procedures import ResistanceParameters
from cmtj.models import LayerDynamic, VectorObj, LayerSB

rp = ResistanceParameters(
    Rxx0=100, Rxy0=1, Rsmr=-0.46, Rahe=-2.7, Ramr=-0.24, l=30, w=20
)
demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
alpha = 0.005
surface = 150e-9 * 150e-9 * pi


@pytest.fixture
def arg(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def two_layer_symbolic_dyn() -> Tuple[LayerDynamic]:
    layerA = LayerDynamic(
        0,
        thickness=2.3e-9,
        Kv=VectorObj(4e3, 0),
        Ks=100,
        Ms=1.8 / mu0,
        alpha=1e-3,
    )
    layerB = LayerDynamic(
        1,
        thickness=3.2e-9,
        Kv=VectorObj(340e3, 0),
        Ks=100,
        Ms=1.2 / mu0,
        alpha=1e-3,
    )
    return (layerA, layerB)


@pytest.fixture
def two_layer_symbolic_classic() -> Tuple[LayerSB]:
    layerA = LayerSB(
        0,
        thickness=2.3e-9,
        Kv=VectorObj(4e3, 0),
        Ks=100,
        Ms=1.8 / mu0,
    )
    layerB = LayerSB(
        1,
        thickness=3.2e-9,
        Kv=VectorObj(340e3, 0),
        Ks=100,
        Ms=1.2 / mu0,
    )
    return (layerA, layerB)


@pytest.fixture
def single_layer_mtj() -> Tuple[Junction, ResistanceParameters]:
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
    junction = Junction([l1])
    junction.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K1))

    return (junction, [rp])


@pytest.fixture
def single_layer_mtj_fictious() -> Tuple[Junction, ResistanceParameters]:
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

    return (junction, None)


@pytest.fixture
def two_layer_mtj() -> Tuple[Junction, ResistanceParameters]:
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

    l2 = Layer(
        "bottom",
        mag=CVector(0, 0, 1),
        anis=Kdir,
        Ms=1.45,
        thickness=4e-9,
        cellSurface=surface,
        demagTensor=demag,
        damping=alpha,
    )
    J1 = -1.78e-3
    J2 = -1.69e-4

    K1 = 1.05e3
    K2 = 50e3

    junction = Junction([l1, l2], 100, 200)

    junction.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K1))

    junction.setLayerAnisotropyDriver("bottom", ScalarDriver.getConstantDriver(K2))
    junction.setIECDriver("free", "bottom", ScalarDriver.getConstantDriver(J1))
    junction.setQuadIECDriver("free", "bottom", ScalarDriver.getConstantDriver(J2))
    return (junction, [rp, rp])


@pytest.fixture
def two_layer_mtj_rz() -> Tuple[Junction, ResistanceParameters]:
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

    l2 = Layer(
        "bottom",
        mag=CVector(0, 0, 1),
        anis=Kdir,
        Ms=1.45,
        thickness=4e-9,
        cellSurface=surface,
        demagTensor=demag,
        damping=alpha,
    )
    J1 = -1.78e-3
    J2 = -1.69e-4

    K1 = 1.05e3
    K2 = 50e3

    junction = Junction([l1, l2], 100, 200)

    junction.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(K1))

    junction.setLayerAnisotropyDriver("bottom", ScalarDriver.getConstantDriver(K2))
    junction.setIECDriver("free", "bottom", ScalarDriver.getConstantDriver(J1))
    junction.setQuadIECDriver("free", "bottom", ScalarDriver.getConstantDriver(J2))
    return (junction, None)


@pytest.fixture
def tri_layer_mtj() -> Tuple[Junction, ResistanceParameters]:
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

    l2 = Layer(
        "middle",
        mag=CVector(0, 0, 1),
        anis=Kdir,
        Ms=1.45,
        thickness=4e-9,
        cellSurface=surface,
        demagTensor=demag,
        damping=alpha,
    )

    l3 = Layer(
        "bottom",
        mag=CVector(0.0, 0, 1.0),
        anis=CVector(0, 0, 1.0),
        Ms=0.7,
        thickness=1e-9,
        cellSurface=surface,
        demagTensor=demag,
        damping=alpha,
    )
    J11 = -1.78e-3
    J12 = -1.69e-4

    J21 = 3e-4
    J22 = 2e-6

    K1 = 1.05e3
    K2 = 50e3
    K3 = 1e4
    junction = Junction([l1, l2, l3])

    for l, anis_val in zip(["free", "middle", "bottom"], [K1, K2, K3]):
        junction.setLayerAnisotropyDriver(l, ScalarDriver.getConstantDriver(anis_val))

    junction.setIECDriver("free", "middle", ScalarDriver.getConstantDriver(J11))
    junction.setQuadIECDriver("free", "middle", ScalarDriver.getConstantDriver(J12))

    junction.setIECDriver("middle", "bottom", ScalarDriver.getConstantDriver(J21))
    junction.setQuadIECDriver("middle", "bottom", ScalarDriver.getConstantDriver(J22))
    return (junction, [rp, rp, rp])
