import math

import pytest

from cmtj import AxialDriver, CVector, Junction, Layer, ScalarDriver


def make_afm_layer(mag1=None, mag2=None, Jafm=-1e-3):
    mag1 = mag1 or CVector(1, 0, 0.01)
    mag2 = mag2 or CVector(-1, 0, -0.01)
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    return Layer.createAFMLayer(
        "afm",
        mag1,
        mag2,
        CVector(1, 0, 0),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
        ScalarDriver.getConstantDriver(Jafm),
    )


def test_afm_layer_is_flagged():
    layer = make_afm_layer()
    assert layer.isAFM
    assert layer.mag2.z == pytest.approx(-0.01, abs=1e-4)


def test_afm_equilibrium_is_static_without_field():
    # antiparallel sublattices with no external field/anisotropy torque
    # should stay at the Neel exchange energy minimum
    layer = make_afm_layer(mag1=CVector(1, 0, 0), mag2=CVector(-1, 0, 0))
    junction = Junction([layer])
    junction.runSimulation(1e-9, 1e-12, 1e-11)
    log = junction.getLog()
    assert log["afm_mx"][-1] == pytest.approx(1.0, abs=1e-6)
    assert log["afm_m2x"][-1] == pytest.approx(-1.0, abs=1e-6)


def test_afm_sublattices_stay_normalised_under_field():
    layer = make_afm_layer()
    layer.setExternalFieldDriver(
        AxialDriver(
            ScalarDriver.getConstantDriver(0),
            ScalarDriver.getConstantDriver(0),
            ScalarDriver.getConstantDriver(1e5),
        )
    )
    junction = Junction([layer])
    junction.runSimulation(1e-9, 1e-12, 1e-11)
    log = junction.getLog()
    for mx, my, mz in zip(log["afm_mx"], log["afm_my"], log["afm_mz"]):
        assert math.sqrt(mx**2 + my**2 + mz**2) == pytest.approx(1.0, abs=1e-6)
    for mx, my, mz in zip(log["afm_m2x"], log["afm_m2y"], log["afm_m2z"]):
        assert math.sqrt(mx**2 + my**2 + mz**2) == pytest.approx(1.0, abs=1e-6)


def test_afm_rejects_temperature_driver():
    layer = make_afm_layer()
    with pytest.raises(RuntimeError, match="stochastic/temperature"):
        layer.setTemperatureDriver(ScalarDriver.getConstantDriver(300))


def test_afm_rejects_temperature_driver_set_before_magnetisation2():
    # order-independence: marking a layer AFM *after* a temperature driver
    # was already set must be rejected too, not just the reverse order.
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    layer = Layer(
        "afm",
        CVector(1, 0, 0),
        CVector(1, 0, 0),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
    )
    layer.setTemperatureDriver(ScalarDriver.getConstantDriver(300))
    with pytest.raises(RuntimeError, match="temperature/noise driver"):
        layer.setMagnetisation2(CVector(-1, 0, 0))


def test_afm_rejects_onef_noise():
    layer = make_afm_layer()
    with pytest.raises(RuntimeError, match="stochastic noise"):
        layer.setOneFNoise(1, 0.0, 0.1)


def test_junction_rejects_stochastic_solver_with_sibling_afm_layer():
    # a temperature driver on a NON-AFM sibling layer must not silently
    # freeze the AFM layer's sublattice 2 under Euler-Heun/Heun.
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    normal = Layer(
        "free",
        CVector(0, 0, 1),
        CVector(0, 0, 1),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
    )
    normal.setTemperatureDriver(ScalarDriver.getConstantDriver(300))
    afm = make_afm_layer()
    junction = Junction([normal, afm])
    with pytest.raises(RuntimeError, match="AFM layer"):
        junction.runSimulation(1e-11, 1e-12, 1e-11)


def test_junction_rejects_dormand_prince_with_afm_layer():
    from cmtj import SolverMode

    layer = make_afm_layer()
    junction = Junction([layer])
    with pytest.raises(RuntimeError, match="RK4"):
        junction.runSimulation(1e-11, 1e-12, 1e-11, solverMode=SolverMode.DormandPrince)


def test_junction_rejects_interlayer_coupling_with_afm_layer():
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    normal = Layer(
        "free",
        CVector(0, 0, 1),
        CVector(0, 0, 1),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
    )
    afm = make_afm_layer()
    junction = Junction([normal, afm])
    with pytest.raises(RuntimeError, match="AFM layer"):
        junction.setIECDriver("free", "afm", ScalarDriver.getConstantDriver(1e-3))


def test_junction_rejects_magnetoresistance_with_afm_layer():
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    normal = Layer(
        "free",
        CVector(0, 0, 1),
        CVector(0, 0, 1),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
    )
    afm = make_afm_layer()
    junction = Junction([normal, afm], 100, 200)
    with pytest.raises(RuntimeError, match="AFM layer"):
        junction.getMagnetoresistance()


def test_junction_all_broadcast_temperature_rejects_without_partial_mutation():
    # a broadcast onto "all" layers must not partially mutate the normal
    # layer before discovering the AFM layer can't take a temperature driver
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    normal = Layer(
        "free",
        CVector(0, 0, 1),
        CVector(0, 0, 1),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
    )
    afm = make_afm_layer()
    junction = Junction([normal, afm])
    with pytest.raises(RuntimeError, match="AFM layer"):
        junction.setLayerTemperatureDriver("all", ScalarDriver.getConstantDriver(300))
    # the normal layer must be untouched -- it should still be usable under RK4
    junction.runSimulation(1e-11, 1e-12, 1e-11)


def test_afm_rejects_sot_torque_drivers():
    layer = make_afm_layer()
    with pytest.raises(RuntimeError, match="SOT/STT"):
        layer.setFieldLikeTorqueDriver(ScalarDriver.getConstantDriver(1e3))


def test_stt_layer_rejects_magnetisation2():
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    layer = Layer.createSTTLayer(
        "stt",
        CVector(0, 0, 1),
        CVector(0, 0, 1),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
        1.0,
        0.0,
        0.5,
    )
    with pytest.raises(RuntimeError, match="STT/SOT"):
        layer.setMagnetisation2(CVector(0, 0, -1))


def test_mr_junction_with_afm_layer_rejects_run():
    # logLayerParams logs R from sublattice 1 only -- running an MR junction
    # containing an AFM layer must fail loudly instead of logging wrong R
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
    normal = Layer(
        "free",
        CVector(0, 0, 1),
        CVector(0, 0, 1),
        0.5,
        1e-9,
        1e-16,
        demag,
        0.01,
    )
    afm = make_afm_layer()
    junction = Junction([normal, afm], 100, 200)
    with pytest.raises(RuntimeError, match="AFM layer"):
        junction.runSimulation(1e-11, 1e-12, 1e-11)


def test_set_layer_magnetisation2_unknown_id_raises():
    layer = make_afm_layer()
    junction = Junction([layer])
    with pytest.raises(RuntimeError, match="Failed to find a layer"):
        junction.setLayerMagnetisation2("nonexistent", CVector(-1, 0, 0))
