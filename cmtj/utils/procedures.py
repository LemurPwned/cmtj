import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
from scipy.fft import fft, fftfreq
from tqdm import tqdm

from cmtj import AxialDriver, Axis, CVector, Junction, NullDriver, ScalarDriver

from .resistance import calculate_resistance_series, compute_sd


@dataclass
class ResistanceParameters:
    """A data holder for resistance parameters. Not all have to be filled in."""

    Rxx0: float = 0
    Rxy0: float = 0
    Rahe: float = 0
    Rsmr: float = 0
    Ramr: float = 0
    w: float = 0  # width
    l: float = 0  # length


def compute_spectrum_strip(input_m: np.ndarray, int_step: float, max_frequency: float):
    """Compute the spectrum of a given magnetization trajectory."""
    yf = np.abs(fft(input_m))
    freqs = fftfreq(len(yf), int_step)
    freqs = freqs[: len(freqs) // 2]
    yf = yf[: len(yf) // 2]

    findx = np.argwhere(freqs <= max_frequency)
    freqs = freqs[findx]
    yf = yf[findx]

    return yf, freqs


def PIMM_procedure(
    junction: "Junction",
    Hvecs: np.ndarray,
    int_step: float,
    resistance_params: list[ResistanceParameters],
    Hoe_direction: Axis = Axis.zaxis,
    Hoe_excitation: float = 50,
    Hoe_duration: int = 3,
    simulation_duration: float = 5e-9,
    wait_time: float = 0e-9,
    max_frequency: float = 80e9,
    resistance_fn: Callable = calculate_resistance_series,
    disturbance: float = 1e-3,
    take_last_n: int = 100,
    full_output: bool = False,
    disable_tqdm: bool = False,
    static_only: bool = False,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Procedure for computing Pulse Induced Microwave Magnetometry.
    It computes both PIMM and Resistance (for instance AHE loops).
    Set `static_only` to True to only compute the static resistance.
    :param junction: junction to be simulated.
    :param Hvecs: list of cartesian vectors. (use FieldScan.amplitude_scan or alike)
    :param int_step: integration step [s].
    :param resistance_params: list of resistance parameters.
    :param Hoe_direction: direction of oersted field (x, y or z).
    :param simulation_duration: duration of simulation [s].
    :param wait_time: time to wait before taking vector for the fft [s].
    :param Hoe_duration: duration of Hoe excitation in multiples of in step
    :param max_frequency: maximum frequency -- larger will be dropped [Hz].
    :param resistance_fn: function to be used to compute the resistance
        (either calculate_resistance_series or calculate_resistance_parallel).
    :param disturbance: disturbance to be applied to the magnetization (std of normal distribution).
    :param take_last_n: number of last time steps to be taken for the compuation.
    :param full_output: if True, return the full trajectories and per layer spectra.
    :param disable_tqdm: if True, disable tqdm progress bar.
    :param static_only: if True, only compute the static resistance.
    :return: (spectrum, frequencies, other_data)
    other_data is a dictionary with the following keys:
    - 'H': Hext field [A/m]
    - 'Rx': resistance in x direction [Ohm]
    - 'Ry': resistance in y direction [Ohm]
    - 'm_avg': average magnetization [unit]
    - 'm_traj': magnetization trajectories [unit]
    """
    if wait_time > simulation_duration:
        raise ValueError("wait_time must be smaller than simulation_duration!")
    spectrum = []
    extraction_m_component = None
    if Hoe_direction == Axis.zaxis:
        extraction_m_component = "z"
        oedriver = AxialDriver(
            NullDriver(),
            NullDriver(),
            ScalarDriver.getStepDriver(0, Hoe_excitation, 0, int_step * Hoe_duration),
        )
    elif Hoe_direction == Axis.yaxis:
        extraction_m_component = "y"
        oedriver = AxialDriver(
            NullDriver(),
            ScalarDriver.getStepDriver(0, Hoe_excitation, 0, int_step * Hoe_duration),
            NullDriver(),
        )
    else:
        extraction_m_component = "x"
        oedriver = AxialDriver(
            ScalarDriver.getStepDriver(0, Hoe_excitation, 0, int_step * Hoe_duration),
            NullDriver(),
            NullDriver(),
        )

    # get layer strings
    layer_ids = junction.getLayerIds()
    if len(layer_ids) != len(resistance_params):
        raise ValueError("The number of layers in the junction must match the number of resistance parameters!")
    output = defaultdict(list)
    normalising_factor = np.sum([layer.thickness * layer.Ms for layer in junction.layers])
    freqs = None  # in case of static_only
    for H in tqdm(Hvecs, desc="Computing PIMM", disable=disable_tqdm):
        junction.clearLog()
        junction.setLayerExternalFieldDriver(
            "all",
            AxialDriver(
                ScalarDriver.getConstantDriver(H[0]),
                ScalarDriver.getConstantDriver(H[1]),
                ScalarDriver.getConstantDriver(H[2]),
            ),
        )
        junction.setLayerOerstedFieldDriver("all", oedriver)
        if disturbance:
            for layer_id in layer_ids:
                old_mag = junction.getLayerMagnetisation(layer_id)
                new_mag = CVector(
                    old_mag.x + np.random.normal(0, disturbance),
                    old_mag.y + np.random.normal(0, disturbance),
                    old_mag.z + np.random.normal(0, disturbance),
                )
                new_mag.normalize()
                junction.setLayerMagnetisation(layer_id, new_mag)
        junction.runSimulation(simulation_duration, int_step, int_step)
        log = junction.getLog()
        indx = np.argwhere(np.asarray(log["time"]) >= wait_time).ravel()
        m_traj = np.asarray(
            [
                np.asarray(
                    [
                        log[f"{layer.id}_mx"],
                        log[f"{layer.id}_my"],
                        log[f"{layer.id}_mz"],
                    ]
                )
                * layer.thickness
                * layer.Ms
                / normalising_factor
                for layer in junction.layers
            ]
        )
        m = m_traj[:, :, -take_last_n:]  # all layers, all x, y, z, last 100 steps
        Rx, Ry = resistance_fn(
            [r.Rxx0 for r in resistance_params],
            [r.Rxy0 for r in resistance_params],
            [r.Ramr for r in resistance_params],
            [r.Rahe for r in resistance_params],
            [r.Rsmr for r in resistance_params],
            m,
            l=[r.l for r in resistance_params],
            w=[r.w for r in resistance_params],
        )
        if not static_only:
            mixed = np.asarray(
                [
                    np.asarray(log[f"{layer.id}_m{extraction_m_component}"])[indx]
                    * layer.thickness
                    * layer.Ms
                    / normalising_factor
                    for layer in junction.layers
                ]
            )
            mixed_sum = mixed.sum(axis=0)
            yf, freqs = compute_spectrum_strip(mixed_sum, int_step, max_frequency)

            spectrum.append(yf)

        # fill the output dict
        output["H"].append(H)
        output["Rx"].append(Rx)
        output["Ry"].append(Ry)
        output["m_avg"].append(m_traj[:, :, -1].sum(0))
        if full_output and not static_only:
            output["m_traj"].append(m_traj)
            for li, layer_id in enumerate(layer_ids):
                y, _ = compute_spectrum_strip(mixed[li], int_step, max_frequency)
                output[layer_id].append(y)
    spectrum = np.squeeze(np.asarray(spectrum))
    if full_output:
        for layer_id in layer_ids:
            output[layer_id] = np.asarray(output[layer_id]).squeeze()
    return spectrum, freqs, output


def VSD_procedure(
    junction: Junction,
    Hvecs: np.ndarray,
    frequencies: np.ndarray,
    int_step: float,
    resistance_params: list[ResistanceParameters] = None,
    Hoe_direction: Axis = Axis.yaxis,
    Hoe_excitation: float = 50,
    simulation_duration: float = 30e-9,
    disturbance: float = 1e-3,
    Rtype: str = "Rz",
    resistance_fn: Callable = calculate_resistance_series,
    disable_tqdm: bool = False,
):
    """Procedure for computing Voltage-Spin Diode.
    We use the Oersted field sine exctitation to excite the system.
    :param junction: junction to be simulated.
    :param Hvecs: list of cartesian vectors. (use FieldScan.amplitude_scan or alike)
    :param frequencies: list of frequencies [Hz].
    :param int_step: integration step [s].
    :param resistance_params: list of resistance parameters.
    :param Hoe_direction: direction of oersted field (x, y or z).
    :param Hoe_excitation: excitation amplitude of Hoe [A/m].
    :param simulation_duration: duration of simulation [s].
    :param disturbance: disturbance to be applied to the magnetization (std of normal distribution).
    :param resistance_fn: function to be used to compute the resistance
        (either calculate_resistance_series or calculate_resistance_parallel). Rz forces standard magnetores.
    :param Rtype: type of resistance to be used. (Rx Ry or Rz)
    :param disable_tqdm: if True, disable tqdm progress bar.
    """
    if resistance_params is None:
        resistance_params = []
    layer_ids = junction.getLayerIds()
    if Rtype == "Rz" and len(layer_ids) > 2:
        raise ValueError("Rz can only be used for 2 layer junctions. Use Rx or Ry instead.")
    elif len(resistance_params) != len(layer_ids):
        raise ValueError("The number of layers in the junction must match the number of resistance parameters!")

    def simulate_VSD(H: np.ndarray, frequency: float, resistance_params: ResistanceParameters):
        if Hoe_direction == Axis.zaxis:
            oedriver = AxialDriver(
                NullDriver(),
                NullDriver(),
                ScalarDriver.getSineDriver(0, Hoe_excitation, frequency, 0),
            )
        elif Hoe_direction == Axis.yaxis:
            oedriver = AxialDriver(
                NullDriver(),
                ScalarDriver.getSineDriver(0, Hoe_excitation, frequency, 0),
                NullDriver(),
            )
        else:
            oedriver = AxialDriver(
                ScalarDriver.getSineDriver(0, Hoe_excitation, frequency, 0),
                NullDriver(),
                NullDriver(),
            )

        junction.clearLog()
        junction.setLayerExternalFieldDriver(
            "all",
            AxialDriver(
                ScalarDriver.getConstantDriver(H[0]),
                ScalarDriver.getConstantDriver(H[1]),
                ScalarDriver.getConstantDriver(H[2]),
            ),
        )
        junction.setLayerOerstedFieldDriver("all", oedriver)
        if disturbance:
            for layer_id in layer_ids:
                old_mag = junction.getLayerMagnetisation(layer_id)
                new_mag = CVector(
                    old_mag.x + np.random.normal(0, disturbance),
                    old_mag.y + np.random.normal(0, disturbance),
                    old_mag.z + np.random.normal(0, disturbance),
                )
                new_mag.normalize()
                junction.setLayerMagnetisation(layer_id, new_mag)
        junction.runSimulation(simulation_duration, int_step, int_step)
        log = junction.getLog()
        m_traj = np.asarray(
            [
                [
                    log[f"{layer_ids[i]}_mx"],
                    log[f"{layer_ids[i]}_my"],
                    log[f"{layer_ids[i]}_mz"],
                ]
                for i in range(len(layer_ids))
            ]
        )
        if Rtype == "Rz":
            if len(layer_ids) > 2:
                raise ValueError("Rz can only be used for 2 layer junctions. One layer can be fictisious.")
            elif len(layer_ids) == 2:
                R = log[f"R_{layer_ids[0]}_{layer_ids[1]}"]
            elif len(layer_ids) == 1:
                R = log["Resistance"]
            else:
                raise ValueError(
                    "Resistance definition ambiguous!"
                    "If you want to use Rz, you must provide"
                    "a single resistance parameter set or set Rp Rap"
                    " at junction creation."
                )
        else:
            Rx, Ry = resistance_fn(
                [r.Rxx0 for r in resistance_params],
                [r.Rxy0 for r in resistance_params],
                [r.Ramr for r in resistance_params],
                [r.Rahe for r in resistance_params],
                [r.Rsmr for r in resistance_params],
                m_traj,
                l=[r.l for r in resistance_params],
                w=[r.w for r in resistance_params],
            )
            if Rtype == "Rx":
                R = Rx
            elif Rtype == "Ry":
                R = Ry
            else:
                raise ValueError("Rtype must be either Rx or Ry or Rz")
        dynamicI = np.sin(2 * math.pi * frequency * np.asarray(log["time"]))
        vmix = compute_sd(R, dynamicI, int_step)
        return vmix

    spectrum = np.zeros((len(Hvecs), len(frequencies)))
    for hindx, H in enumerate(tqdm(Hvecs, "Computing VSD", disable=disable_tqdm)):
        for findx, f in enumerate(frequencies):
            spectrum[hindx, findx] = simulate_VSD(H, f, resistance_params)
    return spectrum
