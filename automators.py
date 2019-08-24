import numpy as np


def voltage_spin_diode_ex(junction,
                          min_field,
                          max_field,
                          simulation_time=10e9,
                          samples=30):
    """
    Field scan
        - scan with the field in range [min_field, max_field]
        - collect the magnetoresistance 
        - calculate average magnetoresistance assuming some current fed into the circuit
        - calculate constant DC voltage
        - aggregate
    """
    field_values = np.linspace(min_field, max_field, samples)
    f = junction.f 
    I = np.sqrt(10/1e3)*np.sin(2*np.pi*f*t)
    Vdcs = []
    for H_ext in field_values:
        junction.set_external_field(H_ext)
        junction.run_simulation(simulation_time)
        R = junction.get_magnetoresistance()
        Vdc = np.mean(R*I)
        Vdcs.append(Vdc)


def fourier_scan():
    pass
