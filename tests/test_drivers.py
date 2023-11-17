from cmtj import AxialDriver, CVector


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
