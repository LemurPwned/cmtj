from junction import Layer, Junction
from constants import Constants
from workers import frequency_analysis, voltage_spin_diode
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

constant = Constants()

dipole_tensor = np.array([
    [6.8353909454237E-4, 0.,0.],
    [0., 0.00150694452305927, 0.],
    [0., 0., 0.99780951638608]
], dtype=np.double)
demag_tensor = np.array([
    [5.57049776248663E-4,0., 0.],
    [0., 0.00125355500286346, 0.],
    [0., 0.0, -0.00181060482770131]
], dtype=np.double)
l1 = Layer(id_="free",
           start_mag=[0.0, 0.0, 1.0],
           start_anisotropy=[0.0, 0.0, 1.0],
           K=900e3,
           Ms=1200e3,
           coupling=-3e-6,
           thickness=1.4e-9)

l1.dipole_tensor = dipole_tensor
l1.demagnetisation_tensor = demag_tensor
l2 = Layer(id_="bottom",
           start_mag=[0., 0.0, 1.0],
           start_anisotropy=[0., 0., 1.0],
           K=1000e3,
           Ms=1000e3,
           coupling=-3e-6,
           thickness=7e-10)
l2.dipole_tensor = dipole_tensor
l2.demagnetisation_tensor = demag_tensor
junction = Junction('MTJ', layers=[l1, l2], couplings=[[2], [1]], persist=True)


def step_field(time, step_start=7e-9, step_stop=7.001e-9):
    Hval = np.zeros((3,))
    if time <= step_stop and time >= step_start:
        Hval[0] = 0.001254*constant.TtoAm
    return Hval
# def step_field(time, step_start=7e-9, step_stop=7.01e-9):
#     Hval = np.zeros((3,))
#     if time <= step_stop and time >= step_start:
#         Hval[0] = 10e-3*constant.TtoAm
#     return Hval


# def anisotropy_update(time):
#     frequency = 6.84e9  # 10 Ghz
#     omega = 2 * math.pi * frequency
#     return 1e3*math.sin(2*omega*time)


def get_resonance_frequency(junction):
    junction.set_junction_global_external_field(
                    250e-3*constant.TtoAm, axis='x')
    junction.run_simulation(10e-9)
    print(frequency_analysis(junction))


def perform_vsd(junction):
    voltage_spin_diode(junction, 0, 400e-3)

# perform_vsd(junction)
get_resonance_frequency(junction)
# # junction.set_global_anisotropy_function(anisotropy_update)
# junction.set_global_field_function(step_field)
# voltage_spin_diode(junction, 0, 400e-3)
# junction.set_junction_global_external_field(
                # 200e-3*constant.TtoAm, axis='x')
# junction.set_global_anisotropy_function(anisotropy_update)
# junction.set_global_field_function(step_field)
# junction.run_simulation(10e-9)
# junction.junction_result[['K_log_free']].plot()
# plt.show()
# junction.junction_result[['Hext_const_x_free',
                        #   'Hext_const_y_free', 'Hext_const_x_bottom', 'Hext_const_y_bottom']].plot()
# junction.junction_result[['Hext_x_free',
                        #   'Hext_y_free', 'Hext_x_bottom', 'Hext_y_bottom']].plot()

# print(frequency_analysis(junction))
# junction.junction_result[['m_x_free',
#                           'm_y_free', 'm_z_free', 
#                           'm_x_bottom', 'm_y_bottom', 'm_z_bottom']].plot()
# junction.junction_result['R_free_bottom'].plot()
# plt.show()
# # 