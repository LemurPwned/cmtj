from cmtj import AxialDriver, CVector, Junction, Layer
from cmtj import constantDriver, sineDriver
import pytest

def test_cvector_operators():
    vec1 = (1.0, 2.0, 3.0)
    vec2 = (1.0, 2.0, 3.0)
    vec3 = (2, 4, 6)

    cvec1 = CVector(*vec1)
    cvec2 = CVector(*vec2)
    cvec3 = CVector(*vec3)

    # test assignment operators
    assert vec1[0] == cvec1[0]
    assert vec1[1] == cvec1[1]
    assert vec1[2] == cvec1[2]

    # test ops
    assert cvec1 == cvec2
    assert cvec1 != cvec3
    assert (cvec1 + cvec2) == cvec3
    assert (cvec3 - cvec1) == cvec1

    # test self-assignment operators
    assert 2 * cvec1 == cvec3
    assert cvec3 * 0.5 == cvec1
    cvec1 += CVector(1.0, 2.0, 3.0)
    assert cvec1 == cvec3
    cvec3 -= cvec2
    assert cvec3 == cvec2


def test_axial_definitions():
    vec = (1.0, 2.0, 3.0)
    cvec = CVector(*vec)
    d1 = AxialDriver(*vec)
    assert d1.getCurrentAxialDrivers(0.0) == cvec
    d1 = AxialDriver(cvec)
    assert d1.getCurrentAxialDrivers(0.0) == cvec


def test_aliases():
    d1 = AxialDriver(constantDriver(1.0), constantDriver(2.0), constantDriver(3.0))
    assert d1.getCurrentAxialDrivers(0.0) == CVector(1.0, 2.0, 3.0)
    assert d1.getCurrentAxialDrivers(1e6) == CVector(1.0, 2.0, 3.0)


def test_driver_ops():
    driver = sineDriver(10, 20, 1, 0)
    assert driver.getCurrentScalarValue(1 / 4) == 30
    driver *= 2
    assert driver.getCurrentScalarValue(1 / 4) == 60

    driver = sineDriver(10, 20, 1, 0)
    driver += 2
    assert driver.getCurrentScalarValue(1 / 4) == 34

    driver = sineDriver(10, 20, 1, 0) * 2
    assert driver.getCurrentScalarValue(1 / 4) == 60

    driver = sineDriver(10, 20, 1, 0) + 2
    assert driver.getCurrentScalarValue(1 / 4) == 34


def test_junction_with_driver():
    Kdir = CVector(1, 0, 0)
    l1 = Layer(
        "free",
        mag=CVector(0, 0, 1),
        anis=Kdir,
        Ms=1.65,
        thickness=3e-9,
        cellSurface=100,
        demagTensor=[CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)],
        damping=0.1,
    )
    K1 = 1.05e3
    l1.setReferenceLayer(CVector(1, 0, 0))
    junction = Junction([l1], 100, 200)
    junction.setLayerAnisotropyDriver("free", constantDriver(K1))
    junction.setLayerExternalFieldDriver(
        "all",
        AxialDriver(constantDriver(0), constantDriver(0), sineDriver(0, 1e3, 7e9, 0)),
    )
    junction.runSimulation(10e-9, 1e-13, 1e-13)
