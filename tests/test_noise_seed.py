import numpy as np

from cmtj import CVector, Junction, Layer, ScalarDriver, SolverMode


def _build_seeded_junction(seed: int, include_temperature: bool) -> Junction:
    demag = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1.0)]
    layer = Layer(
        "free",
        CVector(1.0, 0.0, 0.0),
        CVector(0.0, 0.0, 1.0),
        1.0,
        1.0e-9,
        50e-9 * 50e-9,
        demag,
        damping=0.02,
    )
    junction = Junction([layer])
    junction.setLayerSeed("free", seed)
    junction.setLayerOneFNoise("free", 32, 0.65, 0.02)
    if include_temperature:
        junction.setLayerTemperatureDriver("free", ScalarDriver.getConstantDriver(300.0))
    return junction


def _run_and_get_mz(seed: int, include_temperature: bool) -> np.ndarray:
    junction = _build_seeded_junction(seed, include_temperature)
    junction.runSimulation(5e-10, 1e-12, 1e-12, solverMode=SolverMode.EulerHeun)
    return np.asarray(junction.getLog()["free_mz"])


def _build_seeded_layer(seed: int) -> Layer:
    demag = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1.0)]
    layer = Layer(
        "free",
        CVector(1.0, 0.0, 0.0),
        CVector(0.0, 0.0, 1.0),
        1.0,
        1.0e-9,
        50e-9 * 50e-9,
        demag,
        damping=0.02,
    )
    layer.setSeed(seed)
    layer.setOneFNoise(32, 0.65, 0.02)
    return layer


def test_seed_reproducibility_with_thermal_and_onef_noise():
    same_seed_a = _run_and_get_mz(1234, include_temperature=True)
    same_seed_b = _run_and_get_mz(1234, include_temperature=True)
    different_seed = _run_and_get_mz(5678, include_temperature=True)

    np.testing.assert_allclose(same_seed_a, same_seed_b)
    assert not np.allclose(same_seed_a, different_seed)


def test_seed_reproducibility_with_onef_noise_only():
    layer_a = _build_seeded_layer(2468)
    layer_b = _build_seeded_layer(2468)
    layer_c = _build_seeded_layer(1357)

    same_seed_a = np.asarray([layer_a.getOneFNoise() for _ in range(8)])
    same_seed_b = np.asarray([layer_b.getOneFNoise() for _ in range(8)])
    different_seed = np.asarray([layer_c.getOneFNoise() for _ in range(8)])

    np.testing.assert_allclose(same_seed_a, same_seed_b)
    assert not np.allclose(same_seed_a, different_seed)
