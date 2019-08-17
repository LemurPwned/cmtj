import numpy as np
from junction import Junction


def voltage_spin_diode(junction: Junction, start_h, stop_h, step_h):
    """
    Field scan
        - scan with the field in range [min_field, max_field]
        - collect the magnetoresistance 
        - calculate average magnetoresistance assuming some current fed into the circuit
        - calculate constant DC voltage
        - aggregate
    """
    for h_value in np.linspace(start_h, stop_h, step_h):
        # set the field
        junction.set_junction_global_external_field(h_value)
        # restart simualtion
        junction.restart()
        junction.run_simulation(10e-9, persist=True)
        # extract the magnetisation value
        # calcualte magnetoresistance and get the current
        # junction.
