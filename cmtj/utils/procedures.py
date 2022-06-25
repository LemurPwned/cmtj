import math
from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.fft import fft, fftfreq
from tqdm import tqdm

from cmtj import AxialDriver, Axis, Junction, NullDriver, ScalarDriver

from . import (calculate_magnetoresistance, calculate_resistance_series,
               compute_sd)
from .parallel import distribute


@dataclass
class ResistanceParameters:
    """A data holder for resistance parameters. Not all have to be filled"""
    Rxx0: float = 0
    Rxy0: float = 0
    Rp: float = 100
    Rap: float = 200
    Rahe: float = 0
    Rsmr: float = 0
    Ramr: float = 0
    w: float = 0  # width
    l: float = 0  # length


def PIMM_procedure(
    junction: 'Junction',
    Hvecs: np.ndarray,
    int_step: float,
    resistance_params: List[ResistanceParameters],
    Hoe_direction: Axis = Axis.zaxis,
    Hoe_excitation: float = 50,
    Hoe_duration: int = 3,
    simulation_duration: float = 5e-9,
    max_frequency: float = 80e9,
    output_full_trajectories: bool = False
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Procedure for computing Pulse Induced Microwave Magnetometry
    :param junction: junction to be simulated.
    :param Hvecs: list of cartesian vectors. (use FieldScan.amplitude_scan or alike)
    :param int_step: integration step [s].
    :param resistance_params: list of resistance parameters.
    :param Hoe_direction: direction of oersted field (x, y or z).
    :param simulation_duration: duration of simulation [s].
    :param Hoe_excitation: excitation amplitude of Hoe [A/m].
    :param Hoe_duration: duration of Hoe excitation in multiples of in step
    :param max_frequency: maximum frequency -- larger will be dropped [Hz].
    :param output_full_trajectories: if True, return the full trajectories of the magnets.
    :return: (spectrum, frequencies, other_data)
    other_data is a dictionary with the following keys:
    - 'H': Hext field [A/m]
    - 'Rx': resistance in x direction [Ohm]
    - 'Ry': resistance in y direction [Ohm]
    - 'm_avg': average magnetization [unit]
    - 'm_traj': magnetization trajectories [unit]
    """
    spectrum = []
    extraction_m_component = None
    if Hoe_direction == Axis.zaxis:
        extraction_m_component = 'z'
        oedriver = AxialDriver(
            NullDriver(), NullDriver(),
            ScalarDriver.getStepDriver(0, Hoe_excitation, 0,
                                       int_step * Hoe_duration))
    elif Hoe_direction == Axis.yaxis:
        extraction_m_component = 'y'
        oedriver = AxialDriver(
            NullDriver(),
            ScalarDriver.getStepDriver(0, Hoe_excitation, 0,
                                       int_step * Hoe_duration), NullDriver())
    else:
        extraction_m_component = 'x'
        oedriver = AxialDriver(
            ScalarDriver.getStepDriver(0, Hoe_excitation, 0,
                                       int_step * Hoe_duration), NullDriver(),
            NullDriver())

    # get layer strings
    layer_ids = junction.getLayerIds()
    output = defaultdict(list)
    for H in tqdm(Hvecs, desc="Computing PIMM"):
        junction.clearLog()
        junction.setLayerExternalFieldDriver(
            "all",
            AxialDriver(ScalarDriver.getConstantDriver(H[0]),
                        ScalarDriver.getConstantDriver(H[1]),
                        ScalarDriver.getConstantDriver(H[2])))
        junction.setLayerOerstedFieldDriver("all", oedriver)

        junction.runSimulation(simulation_duration, int_step, int_step)
        log = junction.getLog()

        m_traj = np.asarray([[
            log[f'{layer_ids[i]}_mx'], log[f'{layer_ids[i]}_my'],
            log[f'{layer_ids[i]}_mz']
        ] for i in range(len(layer_ids))])
        m = m_traj[:, :, -100:]  # all layers, all x, y, z, last timestamp
        m_avg = np.mean(m_traj[:, :, -1], 0)
        mixed = [
            np.asarray(log[f"{layer_ids[i]}_m{extraction_m_component}"])
            for i in range(len(layer_ids))
        ]
        mixed = np.mean(np.squeeze(mixed), axis=0)
        yf = np.abs(fft(mixed))
        freqs = fftfreq(len(yf), int_step)
        freqs = freqs[:len(freqs) // 2]
        yf = yf[:len(yf) // 2]

        findx = np.argwhere(freqs <= max_frequency)
        freqs = freqs[findx]
        yf = yf[findx]

        spectrum.append(yf)
        Rx, Ry = calculate_resistance_series(
            [r.Rxx0
             for r in resistance_params], [r.Rxy0 for r in resistance_params],
            [r.Ramr
             for r in resistance_params], [r.Rahe for r in resistance_params],
            [r.Rsmr for r in resistance_params],
            m,
            l=[r.l for r in resistance_params],
            w=[r.w for r in resistance_params])
        # fill the output dict
        output['H'].append(H)
        output['Rx'].append(Rx)
        output['Ry'].append(Ry)
        output['m_avg'].append(m_avg)
        if output_full_trajectories:
            output['m_traj'].append(m_traj)
    spectrum = np.squeeze(np.asarray(spectrum))
    return spectrum, freqs, output


def VSD_procedure(junction: Junction,
                  Hvecs: np.ndarray,
                  frequencies: np.ndarray,
                  int_step: float,
                  resistance_params,
                  Hoe_direction: Axis = Axis.yaxis,
                  Hoe_excitation: float = 50,
                  simulation_duration: float = 30e-9,
                  Rtype: str = 'Rz'):
    """Procedure for computing Voltage-Spin Diode
    :param junction: junction to be simulated.
    :param Hvecs: list of cartesian vectors. (use FieldScan.amplitude_scan or alike)
    :param frequencies: list of frequencies [Hz].
    :param int_step: integration step [s].
    :param resistance_params: list of resistance parameters.
    :param Hoe_direction: direction of oersted field (x, y or z).
    :param Hoe_excitation: excitation amplitude of Hoe [A/m].
    :param simulation_duration: duration of simulation [s].
    :param Rtype: type of resistance to be used. (Rx Ry or Rz)
    """
    layer_ids = junction.getLayerIds()

    def simulate_VSD(H: np.ndarray, frequency: float, resistance_params):

        if Hoe_direction == Axis.zaxis:
            oedriver = AxialDriver(
                NullDriver(), NullDriver(),
                ScalarDriver.getSineDriver(0, Hoe_excitation, frequency, 0))
        elif Hoe_direction == Axis.yaxis:
            oedriver = AxialDriver(
                NullDriver(),
                ScalarDriver.getSineDriver(0, Hoe_excitation, frequency, 0),
                NullDriver())
        else:
            oedriver = AxialDriver(
                ScalarDriver.getSineDriver(0, Hoe_excitation, frequency, 0),
                NullDriver(), NullDriver())

        junction.clearLog()
        junction.setLayerExternalFieldDriver(
            "all",
            AxialDriver(ScalarDriver.getConstantDriver(H[0]),
                        ScalarDriver.getConstantDriver(H[1]),
                        ScalarDriver.getConstantDriver(H[2])))
        junction.setLayerOerstedFieldDriver("all", oedriver)
        junction.runSimulation(simulation_duration, int_step, int_step)
        log = junction.getLog()
        m_traj = np.asarray([[
            log[f'{layer_ids[i]}_mx'], log[f'{layer_ids[i]}_my'],
            log[f'{layer_ids[i]}_mz']
        ] for i in range(len(layer_ids))])
        if Rtype == 'Rz':
            if isinstance(resistance_params,
                          list) and len(resistance_params) > 1:
                resistance_params = resistance_params[0]
            R = log[f'R_{layer_ids[0]}_{layer_ids[1]}']
        else:
            Rx, Ry = calculate_resistance_series(
                [r.Rxx0 for r in resistance_params],
                [r.Rxy0 for r in resistance_params],
                [r.Ramr for r in resistance_params],
                [r.Rahe for r in resistance_params],
                [r.Rsmr for r in resistance_params],
                m_traj,
                l=[r.l for r in resistance_params],
                w=[r.w for r in resistance_params])
            if Rtype == 'Rx':
                R = Rx
            elif Rtype == 'Ry':
                R = Ry
            else:
                raise ValueError("Rtype must be either Rx or Ry or Rz")
        dynamicI = np.sin(2 * math.pi * frequency * np.asarray(log['time']))
        vmix = compute_sd(R, dynamicI, int_step)
        return vmix

    spectrum = np.zeros((len(Hvecs), len(frequencies)))
    for hindx, H in enumerate(tqdm(Hvecs, "Computing VSD")):
        for findx, f in enumerate(frequencies):
            spectrum[hindx, findx] = simulate_VSD(H, f, resistance_params)
    return spectrum
