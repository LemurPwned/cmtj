from junction import Layer, Junction
from constants import Constants
from workers import frequency_analysis, voltage_spin_diode
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

constant = Constants()

l1 = Layer(id_="free",
           start_mag=[0.0, 0.0, 1.0],
           start_anisotropy=[0.0, 0.0, 1.0],
           K=900e3,
           Ms=1200e3,
           coupling=3e-6,
           thickness=2e-9)


dipole_tensor = [[0.00706885975425371, 0.0, 0., ],
                 [0.0, 0.00706885975425371, 0., ],
                 [0., 0., -0.0141377195085074]]
demag_tensor = [[0.00832453381627329, 0., 0.],
                [0., 0.00832453381627329, 0.],
                [0.0, 0.0, 0.765750932367453]]
l1.dipole_tensor = dipole_tensor
l1.demagnetisation_tensor = demag_tensor

# l1.dipole_tensor = np.array([[0.00706885975425371, 0.0, 0., ],
#                              [0.0, 0.00706885975425371, 0., ],
#                              [0., 0., -0.0141377195085074]])
# l1.demagnetisation_tensor = np.array([[0.00832453381627329, 0., 0.],
#                                       [0., 0.00832453381627329, 0.],
#                                       [0.0, 0.0, 0.765750932367453]])
l2 = Layer(id_="bottom",
           start_mag=[0., 0.0, 1.0],
           start_anisotropy=[0., 0., 1.0],
           K=1000e3,
           Ms=1000e3,
           coupling=3e-6,
           thickness=2e-9)
l2.dipole_tensor = dipole_tensor
l2.demagnetisation_tensor = demag_tensor

# l2.dipole_tensor = np.array([[0.00706885975425371, 0.0, 0., ],
#                              [0.0, 0.00706885975425371, 0., ],
#                              [0., 0., -0.0141377195085074]])
# l2.demagnetisation_tensor = np.array([[0.00832453381627329, 0., 0.],
#                                       [0., 0.00832453381627329, 0.],
#                                       [0.0, 0.0, 0.765750932367453]])
junction = Junction('MTJ',
                    layers=[l1, l2], couplings=[[2], [1]], persist=True)


# def step_field(time, step_start=7e-9, step_stop=7.01e-9):
#     Hval = np.zeros((3,))
#     if time <= step_stop and time >= step_start:
#         Hval[0] = 10e-3*constant.TtoAm
#     return Hval


# def anisotropy_update(time):
#     frequency = 6.84e9  # 10 Ghz
#     omega = 2 * math.pi * frequency
#     return 1e3*math.sin(2*omega*time)

def coupling_update(time):
    frequency = 6.93e9  # 10 Ghz
    omega = 2 * math.pi * frequency
    return 8e-7 * math.sin(omega * time)


voltage_spin_diode(junction, 0, 400e-3)


# junction.set_junction_global_external_field(
#     200e-3*constant.TtoAm, axis='x')
# junction.set_global_anisotropy_function(anisotropy_update)
# junction.set_global_field_function(step_field)
# junction.run_simulation(10e-9)
