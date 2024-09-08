from typing import Union, list

import numpy as np

from .filters import Filters


def compute_sd(dynamic_r: np.ndarray, dynamic_i: np.ndarray, integration_step: float) -> np.ndarray:
    """Computes the SD voltage.
    :param dynamic_r: magnetoresistance from log
    :param dynamic_i: excitation current
    :param integration_step: integration paramemter from run_simulation
    """
    SD = -dynamic_i * dynamic_r
    fs = 1.0 / integration_step
    SD_dc = Filters.butter_lowpass_filter(SD, cutoff=10e6, fs=fs, order=3)
    return np.mean(SD_dc)


def compute_resistance(
    Rx0: list[float],
    Ry0: list[float],
    AMR: list[float],
    AHE: list[float],
    SMR: list[float],
    m: Union[list[float], np.ndarray],
    l: list[float],
    w: list[float],
):
    """Computes the resistance of the system.
    If you want to compute the resistance for an entire time series, pass m as a 3D array
    with shape [number_of_layers, 3, T], where T is the time component.
    [number_of_layers, 3, T] where T is the time component.
    """
    number_of_layers = len(Rx0)
    if not isinstance(m, np.ndarray):
        m = np.asarray(m)
    if m.ndim == 2:
        SxAll = np.zeros((number_of_layers,))
        SyAll = np.zeros((number_of_layers,))

    elif m.ndim == 3:
        SxAll = np.zeros((number_of_layers, m.shape[2]))
        SyAll = np.zeros((number_of_layers, m.shape[2]))

    for i in range(number_of_layers):
        w_l = w[i] / l[i]
        SxAll[i] = Rx0[i] + (AMR[i] * m[i, 0] ** 2 + SMR[i] * m[i, 1] ** 2)
        SyAll[i] = Ry0[i] + 0.5 * AHE[i] * m[i, 2] + (w_l) * (SMR[i] - AMR[i]) * m[i, 0] * m[i, 1]
    return SxAll, SyAll


def compute_gmr(Rp: float, Rap: float, m1: np.ndarray, m2: np.ndarray):
    """Computes the GMR using parallel and antiparallel resistance.
    :param Rp: parallel resistance
    :param Rap: antiparallel resistance
    :param m1: magnetisation of layer 1
    :param m2: magnetisation of layer 2"""
    return Rp + 0.5 * (Rap - Rp) * np.sum(m1 * m2, axis=0)


def calculate_magnetoresistance(Rp: float, Rap: float, m: np.ndarray):
    """Computes the magnetoresistance using parallel and antiparallel resistance.
    :param Rp: parallel resistance
    :param Rap: antiparallel resistance
    :param m: magnetisation, 2 layers of shape [2, 3, T] where T is the time component
    """
    if not isinstance(m, np.ndarray):
        m = np.asarray(m)
    if m.shape[0] != 2:
        raise ValueError("The magnetoresistance can only be computed for 2 layers" f". Current shape {m.shape}")
    return Rp + 0.5 * (Rap - Rp) * np.sum(m[0] * m[1], axis=0)


def calculate_resistance_parallel(
    Rx0: list[float],
    Ry0: list[float],
    AMR: list[float],
    AHE: list[float],
    SMR: list[float],
    m: list[float],
    l: list[float],
    w: list[float],
):
    """Calculates the resistance of the system in parallel.
    If you want to compute the resistance for an entire time series, pass m as a 3D array.
    [number_of_layers, 3, T] where T is the time component.
    Uses Kim's formula from the paper:
    https://link.aps.org/doi/10.1103/PhysRevLett.116.097201

    :param Rx0: resistance offset in longitudinal direction
    :param Ry0: resistance offset in transverse direction
    :param AMR: anisotropic magnetoresistance
    :param AHE: anomalous Hall effect
    :param SMR: spin Hall magnetoresistance
    :param m: magnetisation of the layers. Shape [number_of_layers, 3, T]
    :param l: length of the layers
    :param w: width of the layers
    """
    SxAll, SyAll = compute_resistance(Rx0, Ry0, AMR, AHE, SMR, m, l, w)
    Rx = 1 / np.sum(1.0 / SxAll, axis=0)
    Ry = 1 / np.sum(1.0 / SyAll, axis=0)
    return Rx, Ry


def calculate_resistance_series(
    Rx0: list[float],
    Ry0: list[float],
    AMR: list[float],
    AHE: list[float],
    SMR: list[float],
    m: list[float],
    l: list[float],
    w: list[float],
):
    """Calculates the resistance of the system in series.
    If you want to compute the resistance for an entire time series, pass m as a 3D array.
    [number_of_layers, 3, T] where T is the time component.
    Uses Kim's formula from the paper:
    https://link.aps.org/doi/10.1103/PhysRevLett.116.097201

    :param Rx0: resistance offset in longitudinal direction
    :param Ry0: resistance offset in transverse direction
    :param AMR: anisotropic magnetoresistance
    :param AHE: anomalous Hall effect
    :param SMR: spin Hall magnetoresistance
    :param m: magnetisation of the layers. Shape [number_of_layers, 3, T]
    :param l: length of the layers
    :param w: width of the layers
    """
    SxAll, SyAll = compute_resistance(Rx0, Ry0, AMR, AHE, SMR, m, l, w)
    Rx = np.sum(SxAll, axis=0)
    Ry = np.sum(SyAll, axis=0)
    return Rx, Ry
