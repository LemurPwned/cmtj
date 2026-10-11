from itertools import product
from typing import Tuple

import numpy as np
import pytest

from cmtj import Axis, CVector, Junction, Layer
from cmtj.utils import FieldScan
from cmtj.utils.procedures import (PIMM_procedure, ResistanceParameters,
                                   VSD_procedure)
from cmtj.utils.resistance import (calculate_magnetoresistance,
                                   calculate_resistance_parallel,
                                   calculate_resistance_series)


def product_dict(inp):
    return (dict(zip(inp.keys(), values)) for values in product(*inp.values()))


def test_base_pimm(single_layer_mtj: Tuple[Junction, ResistanceParameters]):
    junction, res = single_layer_mtj
    _, Hvecs = FieldScan.amplitude_scan(-500e3,
                                        500e3,
                                        steps=50,
                                        theta=90,
                                        phi=0)
    parameter_variants = {
        "static_only": [True, False],
        "full_output": [True, False],
        "resistance_fn":
        [calculate_resistance_series, calculate_resistance_parallel],
        "Hoe_direction": [Axis.xaxis, Axis.yaxis]
    }
    for param_inp in product_dict(parameter_variants):
        spec, f, out = PIMM_procedure(junction,
                                      Hvecs=Hvecs,
                                      resistance_params=res,
                                      int_step=1e-12,
                                      simulation_duration=10e-9,
                                      disable_tqdm=True,
                                      **param_inp)
        if param_inp["static_only"]:
            assert f is None
            assert len(spec) == 0
        elif param_inp['full_output']:
            assert 'm_traj' in out


@pytest.mark.parametrize(
    'arg', ['single_layer_mtj', 'two_layer_mtj', 'tri_layer_mtj'],
    indirect=True)
def test_pimm(arg: Tuple[Junction, ResistanceParameters]):
    junction, res = arg
    _, Hvecs = FieldScan.amplitude_scan(-500e3,
                                        500e3,
                                        steps=50,
                                        theta=90,
                                        phi=0)
    PIMM_procedure(junction,
                   Hvecs=Hvecs,
                   int_step=1e-12,
                   resistance_params=res,
                   simulation_duration=10e-9,
                   disable_tqdm=True)


@pytest.mark.parametrize(
    'arg', ['single_layer_mtj', 'two_layer_mtj', 'tri_layer_mtj'],
    indirect=True)
def test_vsd(arg: Tuple[Junction, ResistanceParameters]):
    junction, res = arg
    _, Hvecs = FieldScan.amplitude_scan(-500e3,
                                        500e3,
                                        steps=50,
                                        theta=90,
                                        phi=0)
    # TODO: add more tests for Rz
    for rtype in ['Rx']:
        VSD_procedure(junction,
                      Hvecs=Hvecs,
                      resistance_params=res,
                      frequencies=np.arange(1e9, 15e9, 1e9),
                      int_step=1e-12,
                      simulation_duration=10e-9,
                      Rtype=rtype,
                      disable_tqdm=True)


def test_pimm_resistance_uses_unit_magnetisation():
    demag = [CVector(0, 0, 0)] * 3
    layers = [
        Layer(f"L{i}", CVector(1, 0, 0), CVector(1, 0, 0), 1.0, 1e-9, 1e-16, demag, 0.01)
        for i in range(2)
    ]
    junction = Junction(layers)
    rp = ResistanceParameters(Rxx0=100, Ramr=1, w=1, l=1)
    _, _, out = PIMM_procedure(
        junction,
        np.asarray([[500e3, 0, 0]]),
        1e-12,
        [rp, rp],
        simulation_duration=1e-10,
        disturbance=0,
        static_only=True,
        disable_tqdm=True,
    )
    # both layers saturated along x: 101 Ohm each, in series
    assert np.ravel(out["Rx"])[-1] == pytest.approx(202, rel=1e-6)


def test_vsd_rz_single_layer():
    demag = [CVector(0, 0, 0)] * 3
    layer = Layer("free", CVector(1, 0, 0), CVector(1, 0, 0), 1.0, 1e-9, 1e-16, demag, 0.01)
    layer.setReferenceLayer(CVector(1, 0, 0))
    junction = Junction([layer], 100, 200)
    spectrum = VSD_procedure(
        junction,
        np.asarray([[500e3, 0, 0]]),
        [1e9],
        1e-12,
        Rtype="Rz",
        simulation_duration=1e-10,
        disable_tqdm=True,
    )
    assert spectrum.shape == (1, 1)


def test_calculate_magnetoresistance_limits():
    parallel = np.asarray([[[0], [0], [1]], [[0], [0], [1]]], dtype=float)
    antiparallel = np.asarray([[[0], [0], [1]], [[0], [0], [-1]]], dtype=float)
    assert calculate_magnetoresistance(100, 200, parallel)[0] == pytest.approx(100)
    assert calculate_magnetoresistance(100, 200, antiparallel)[0] == pytest.approx(200)
