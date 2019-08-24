from junction import Layer, Junction
from constants import Constants
import cProfile
import numpy as np

constant = Constants()

import pyximport
pyximport.install()

if __name__ == "__main__":
    # j = Junction.from_json('junction.json', persist=True)
    # j.run_simulation(2e-9)
    # print(np.mean(j.junction_result['time']))
    l1 = Layer(id_="free",
               start_mag=[0.8, 0.1, 0.],
               start_anisotropy=[0., 1., 0],
               K=800e3,
               Ms=1200e3,
               thickness=2e-9)
    l2 = Layer(id_="bottom",
               start_mag=[0., 1.0, 0.],
               start_anisotropy=[0., 1., 0],
               K=1500e3,
               Ms=1200e3,
               thickness=2e-9)
    junction = Junction('MTJ', layers=[l1, l2], couplings=None)

    # def field_change(Hext, t):
    #     return np.array([
    #         200e-3 * constant.TtoAm +
    #         20e-3 * constant.TtoAm * np.sin(5e9 * 2 * np.pi * t), 0., 0.
    #     ])

    # l1.update_external_field = field_change
    # junction.set_junction_global_external_field(200)
    # junction.overwrite_snapshot()
    # junction.to_json('junction.json')
    # junction.set_junction_global_external_field(100)
    # # junction.restart()
    # junction.to_json('junction2.json')
    cProfile.runctx('junction.run_simulation(10e-9)', globals(), locals(),
                    "c_test_run.prof")
    # (10e-9)÷
    # junction.run_simulation(5e-9)
    # plot_results()
