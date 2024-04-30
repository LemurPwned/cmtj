from cmtj.llgb import LLGBLayer, LLGBJunction
from cmtj import CVector, ScalarDriver


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
