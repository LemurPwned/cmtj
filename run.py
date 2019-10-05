import pyximport
from junction import Layer, Junction
from constants import Constants
import cProfile
import numpy as np
import math

constant = Constants()

pyximport.install()

if __name__ == "__main__":
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

    l2 = Layer(id_="bottom",
               start_mag=[0., 0.0, 1.0],
               start_anisotropy=[0., 0., 1.0],
               K=1000e3,
               Ms=1000e3,
               coupling=3e-6,
               thickness=2e-9)
    l2.dipole_tensor = dipole_tensor
    l2.demagnetisation_tensor = demag_tensor
    junction = Junction('MTJ', layers=[l1, l2], couplings=[
                        [2], [1]], persist=True)

    def coupling_update(time):
        frequency = 6.93e9  # 10 Ghz
        omega = 2 * math.pi * frequency
        return 8e-7 * math.sin(omega * time)
    junction.set_junction_global_external_field(400e-3*constant.TtoAm,
                                                axis='x')
    junction.set_global_coupling_function(coupling_update)

    cProfile.runctx('junction.run_simulation(15e-9)', globals(), locals(),
                    "c_test_run.prof")
