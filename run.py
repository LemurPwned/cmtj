import pyximport
from junction import Layer, Junction
from constants import Constants
import cProfile
import numpy as np

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

    l1.dipole_tensor = np.array([[0.00706885975425371, 0.0, 0., ],
                                 [0.0, 0.00706885975425371, 0., ],
                                 [0., 0., -0.0141377195085074]], dtype=np.double
                                )
    l1.demagnetisation_tensor = np.array([[0.00832453381627329, 0., 0.],
                                          [0., 0.00832453381627329, 0.],
                                          [0.0, 0.0, 0.765750932367453]], dtype=np.double)
    l2 = Layer(id_="bottom",
               start_mag=[0., 0.0, 1.0],
               start_anisotropy=[0., 0., 1.0],
               K=1000e3,
               Ms=1000e3,
               coupling=3e-6,
               thickness=2e-9)
    l2.dipole_tensor = np.array([[0.00706885975425371, 0.0, 0., ],
                                 [0.0, 0.00706885975425371, 0., ],
                                 [0., 0., -0.0141377195085074]], dtype=np.double)
    l2.demagnetisation_tensor = np.array([[0.00832453381627329, 0., 0.],
                                          [0., 0.00832453381627329, 0.],
                                          [0.0, 0.0, 0.765750932367453]], dtype=np.double)
    junction = Junction('MTJ', layers=[l1, l2], couplings=[
                        [2], [1]], persist=True)

    cProfile.runctx('junction.run_simulation(10e-9)', globals(), locals(),
                    "c_test_run.prof")
