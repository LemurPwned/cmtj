import numpy as np
from junction import Junction


def voltage_spin_diode(junction: Junction, start_h, stop_h, step_h):

    for h_value in np.linspace(start_h, stop_h, step_h):
        # set the field
        junction.layers_id['free'].set_global_external_field_value = h_value
        # restart simualtion
        junction.restart()
        junction.run_simulation(10e-9, persist=True)
        # extract the magnetisation value
        # calcualte magnetoresistance and get the current
        junction.
