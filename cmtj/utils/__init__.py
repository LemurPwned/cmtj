import numpy as np
from .filters import Filters
from .linear import FieldScan
# constants
OetoAm = 79.57747
AmtoOe = 1. / OetoAm
TtoAm = 795774.715459
AmtoT = 1.0 / TtoAm


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


__all__ = ["Filters", "FieldScan"]
