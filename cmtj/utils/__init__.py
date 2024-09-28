import math

from .filters import Filters
from .general import VectorObj, box_muller_random, perturb_position
from .linear import FieldScan
from .resistance import (
    calculate_magnetoresistance,
    calculate_resistance_parallel,
    calculate_resistance_series,
    compute_gmr,
    compute_resistance,
    compute_sd,
)

# constants
OetoAm = 79.57747
AmtoOe = 1.0 / OetoAm
TtoAm = 795774.715459
AmtoT = 1.0 / TtoAm
echarge = -1.602e-19
mu0 = 12.566e-7
hplanck = 6.6260e-34
hbar = hplanck / (2 * math.pi)
gyromagnetic_ratio = 2.211e5  # m/As
gamma = 28024e6  # Hz/T
gamma_rad = 1.76e11  # rad / (s * T)
me = 9.109e-31
bohr_magneton = echarge * hbar / (2 * me)

__all__ = [
    "Filters",
    "FieldScan",
    "compute_sd",
    "compute_resistance",
    "calculate_magnetoresistance",
    "calculate_resistance_series",
    "calculate_resistance_parallel",
    "VectorObj",
    "box_muller_random",
    "perturb_position",
    "compute_gmr",
]
