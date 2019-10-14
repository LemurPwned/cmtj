from pymtj.junction import Layer, Junction
from pymtj.constants import Constants
from pymtj.automation.workers import frequency_analysis, voltage_spin_diode
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

constant = Constants()

dipole_tensor = [[6.8353909454237E-4, 0., 0.], [0., 0.00150694452305927, 0.],
                 [0., 0., 0.99780951638608]]
demag_tensor = [[5.57049776248663E-4, 0., 0.], [0., 0.00125355500286346, 0.],
                [0., 0.0, -0.00181060482770131]]

l1 = Layer(id_="free",
           start_mag=[0.0, 0.0, 1.0],
           anisotropy=[0.0, 0.0, 1.0],
           K=900e3,
           Ms=1200e3,
           coupling=-3e-6,
           thickness=1.4e-9,
           demag_tensor=demag_tensor,
           dipole_tensor=dipole_tensor)

l2 = Layer(id_="bottom",
           start_mag=[0., 0.0, 1.0],
           anisotropy=[0., 0., 1.0],
           K=1000e3,
           Ms=1000e3,
           coupling=-3e-6,
           thickness=7e-10,
           demag_tensor=demag_tensor,
           dipole_tensor=dipole_tensor)

junction = Junction('MTJ', layers=[l1, l2], couplings=[[2], [1]], persist=True)


def step_field(time, step_start=5e-9, step_stop=5.001e-9):
    Hval = np.zeros((3, ))
    if time <= step_stop and time >= step_start:
        Hval[0] = 0.001254 * constant.TtoAm
    return Hval


def coupling_update(time):
    frequency = 6.93e9  # 10 Ghz
    omega = 2 * math.pi * frequency
    return 8e-7 * math.sin(omega * time)


def get_resonance_frequency(junction: Junction):
    # junction.set_global_field_function(step_field)
    junction.set_junction_global_external_field(250e-3 * constant.TtoAm,
                                                axis='x')
    junction.run_simulation(15e-9)
    print(frequency_analysis(junction))


def display_results():
    df = pd.read_csv('results.csv')
    plt.plot(df['time'], df[['m_x_free', 'm_y_free', 'm_z_free']])
    plt.show()


def perform_vsd(junction):
    voltage_spin_diode(junction, 0, 400e-3)


# perform_vsd(junction)
get_resonance_frequency(junction)