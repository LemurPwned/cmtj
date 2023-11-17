from itertools import product
from typing import Tuple

import numpy as np
import pytest

from cmtj import Axis, Junction
from cmtj.utils import FieldScan
from cmtj.utils.procedures import (PIMM_procedure, ResistanceParameters,
                                   VSD_procedure)
from cmtj.utils.resistance import (calculate_resistance_parallel,
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
