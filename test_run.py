from junction import Layer, Junction, plot_results
from constants import Constants
from workers import frequency_analysis, voltage_spin_diode
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
constant = Constants()

l1 = Layer(id_="free",
           start_mag=[0.0, 0.0, 1.0],
           start_anisotropy=[0.0, 0.0, 1.],
           K=900e3,
           Ms=1200e3,
           thickness=2e-9)

l1.dipole_tensor = np.array([[0.00706885975425371, 0.0, 0., ],
                             [0.0, 0.00706885975425371, 0., ],
                             [0., 0., -0.0141377195085074]])
l1.demagnetisation_tensor = np.array([[0.00832453381627329, 0., 0.],
                                      [0., 0.00832453381627329, 0.],
                                      [0.0, 0.0, 0.765750932367453]])
l2 = Layer(id_="bottom",
           start_mag=[0., 0.0, 1.0],
           start_anisotropy=[0., 0., 1.0],
           K=1500e3,
           Ms=1000e3,
           thickness=2e-9)
l2.dipole_tensor = np.array([[0.00706885975425371, 0.0, 0., ],
                             [0.0, 0.00706885975425371, 0., ],
                             [0., 0., -0.0141377195085074]])
l2.demagnetisation_tensor = np.array([[0.00832453381627329, 0., 0.],
                                      [0., 0.00832453381627329, 0.],
                                      [0.0, 0.0, 0.765750932367453]])
junction = Junction('MTJ', layers=[l1, l2], couplings=None, persist=True)


def step_field(time, step_start=5e-9, step_stop=5.01e-9):
    Hval = np.zeros((3,))
    if time < step_stop and time > step_start:
        Hval[0] = 10e-3*constant.TtoAm
    return Hval


# voltage_spin_diode(junction, 0, 500e-3)
df = pd.read_csv('voltage_spin-diode.csv')
plt.plot(df['H'], df['Vmix'], '.')
plt.show()
# junction.set_global_field_function(step_field)
# junction.set_junction_global_external_field(300e-3*constant.TtoAm, axis='x')
# junction.run_simulation(8e-9)
# junction.junction_result[['Hext_const_x_free',
#                           'Hext_const_y_free', 'Hext_const_x_bottom', 'Hext_const_y_bottom']].plot()
# junction.junction_result[['Hext_x_free',
#                           'Hext_y_free', 'Hext_x_bottom', 'Hext_y_bottom']].plot()

# print(frequency_analysis(junction))
# junction.junction_result[['m_x_free',
#                           'm_y_free', 'm_z_free', 'm_x_bottom', 'm_y_bottom', 'm_z_bottom']].plot()
# plt.show()
