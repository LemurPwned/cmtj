import numpy as np

from .filters import Filters
from .linear import FieldScan

# constants
OetoAm = 79.57747
AmtoOe = 1.0 / OetoAm
TtoAm = 795774.715459
AmtoT = 1.0 / TtoAm
echarge = 1.602e-19
mi0 = 12.566e-7
hplanck = 6.6260e-34
hbar = hplanck / (2 * np.pi)


def compute_sd(dynamic_r: np.ndarray, dynamic_i: np.ndarray,
               integration_step: float) -> np.ndarray:
    """Computes the SD voltage
    :param dynamic_r: magnetoresistance from log
    :param dynamic_i: excitation current
    :param integration_step: integration paramemter from run_simulation
    """
    SD = -dynamic_i * dynamic_r
    fs = 1.0 / integration_step
    SD_dc = Filters.butter_lowpass_filter(SD, cutoff=10e6, fs=fs, order=3)
    return np.mean(SD_dc)


def calculate_resistance(Rx0, Ry0, AMR, AHE, SMR, m, number_of_layers, l, w):
    if m.ndim == 2:
        SxAll = np.zeros((number_of_layers, ))
        SyAll = np.zeros((number_of_layers, ))

    elif m.ndim == 3:
        SxAll = np.zeros((number_of_layers, m.shape[2]))
        SyAll = np.zeros((number_of_layers, m.shape[2]))

    for i in range(0, number_of_layers):
        w_l = w[i] / l[i]
        SxAll[i] = 1 / (Rx0[i] + (AMR[i] * m[i, 0]**2 + SMR[i] * m[i, 1]**2))
        SyAll[i] = 1 / (Ry0[i] + 0.5 * AHE[i] * m[i, 2] + (w_l) *
                        (SMR[i] - AMR[i]) * m[i, 0] * m[i, 1])

    Rx = 1 / np.sum(SxAll, axis=0)
    Ry = 1 / np.sum(SyAll, axis=0)
    return Rx, Ry


__all__ = ["Filters", "FieldScan"]
