import math

from .filters import Filters
from .linear import FieldScan
from .resistance import *

# constants
OetoAm = 79.57747
AmtoOe = 1.0 / OetoAm
TtoAm = 795774.715459
AmtoT = 1.0 / TtoAm
echarge = 1.602e-19
mu0 = 12.566e-7
hplanck = 6.6260e-34
hbar = hplanck / (2 * math.pi)
gyromagnetic_ratio = 2.211e5

__all__ = [
    "Filters",
    "FieldScan",
    "compute_sd",
    "compute_resistance",
    "calculate_magnetoresistance",
    "calculate_resistance_series",
    "calculate_resistance_parallel",
]
