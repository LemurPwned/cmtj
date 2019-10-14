from pymtj.junction import Layer, Junction
from pymtj.constants import Constants
from pymtj.automation.workers import frequency_analysis, voltage_spin_diode
from multiprocessing import Pool
from itertools import repeat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import os

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


def calculate_single_voltage(h_value, junction: Junction, frequency, 
                             constant_coupling, coupling_excitation):

    # set the field
    phase_shift = 0
    power = 10e-6
    omega = 2 * np.pi * frequency

    def coupling_update(time):
        return coupling_excitation * np.sin(omega * time)

    junction.restart()
    junction.save = False
    junction.set_global_coupling(constant_coupling)
    junction.set_junction_global_external_field(float(h_value *
                                                      constant.TtoAm),
                                                axis='x')
    junction.set_global_coupling_function(coupling_update)
    # restart simulation
    junction.run_simulation(20e-9)
    # extract the magnetisation value
    # wait for 5ns
    limited_res = junction.junction_result[
        junction.junction_result['time'] >= 15e-9]
    avg_resistance = np.mean(limited_res['R_free_bottom'])
    # print(f"Avg resistance {avg_resistance}")
    amplitude = np.sqrt(power / avg_resistance)
    current = amplitude * np.sin(omega * limited_res['time'] + phase_shift)
    voltage = limited_res['R_free_bottom'] * current

    return h_value, np.mean(voltage)


def parameter_scan(junction: Junction):
    """
    Scan a range of parameters and 
    obtain the function FFT amplitude ()
    """

    # define the step function
    def step_field(time, step_start=5e-9, step_stop=5.001e-9):
        Hval = [0., 0. ,0.]
        if time <= step_stop and time >= step_start:
            Hval[1] = 0.001254 * constant.TtoAm
        return Hval

    couplings = [-8e-6, -5e-6, -3e-6, -1e-6,
                 -8e-5, -5e-5, -3e-5, -1e-5,
                 -8e-4, -5e-4, -3e-4, -1e-4]
    experiment_folder = "experiment2"
    try:
        os.mkdir(experiment_folder)
    except:
        # folder exists
        pass
    for K_top in [800e3, 850e3, 900e3, 950e3]:
        print(f"Calcuating {K_top}...")
        x_f, y_f, z_f = [], [], []
        x_amp, y_amp, z_amp = [], [], []
        for coupling in couplings:
            curr_coupling = float(coupling)
            junction.restart()
            junction.save = False
            junction.persist = True
            junction.set_layer_K("free", K_top)
            junction.set_global_coupling(curr_coupling)
            junction.set_global_field_function(step_field)
            junction.set_junction_global_external_field(250e-3 * constant.TtoAm,
                                                        axis='x')
            junction.run_simulation(15e-9)
            frequencies = np.around(frequency_analysis(junction), 3)
            xf, yf, zf = frequencies[:, 0].tolist()  # only freq, no amp
            xamp, yamp, zamp = frequencies[:, 1].tolist()
            x_f.append(xf)
            y_f.append(yf)
            z_f.append(zf)
            x_amp.append(xamp)
            y_amp.append(yamp)
            z_amp.append(zamp)

        # draw dispersion plots
        for name, freq in zip(['x', 'y', 'z'], [x_f, y_f, z_f]):
            fig = plt.figure()
            plt.plot(couplings, freq, 'o')
            title = f"Amp_plot_{name}_K_{K_top}"
            fig.suptitle(title, fontsize=12)
            plt.xlabel("Coupling [J/m^2]")
            plt.ylabel("Amplitude of resonance frequency")
            plt.savefig(os.path.join(experiment_folder, title + '.png'))

        # save all parameters
        df = pd.DataFrame.from_dict({
            'coupling': couplings,
            'x_f': x_f,
            'y_f': y_f,
            'z_f': z_f,
            'x_amp': x_amp,
            'y_amp': y_amp,
            'z_amp': z_amp
        })
        df.to_csv(os.path.join(experiment_folder, f"K_{K_top}.csv"))
        junction.to_json(os.path.join(experiment_folder, f"junction_base.json"))

def resonance_curve_scan(junction, experiment_folder):
    try:
        os.mkdir(experiment_folder)
    except:
        # folder exists
        pass
    h_vals = np.linspace(0, 600e-3, 60)
    for resonance_frequency, constant_coupling in zip([7e9,7e9, 7e9],
                                                    [-3e-6, -5e-6, -8e-5]):
        subfolder_name = os.path.join(experiment_folder, 
                    f"Freq_{resonance_frequency}_coupling_{constant_coupling}")
        print(subfolder_name + "...")
        os.makedirs(subfolder_name,exist_ok=True)
        for coupling_excitation in [1e-7, 5e-7, 8e-7, 2e-6]:
            with Pool(4) as pool:
                hvals_voltages = pool.starmap(
                    calculate_single_voltage,
                    zip(h_vals, repeat(junction), repeat(resonance_frequency),
                        repeat(constant_coupling), repeat(coupling_excitation)))
            h_vals, voltages = zip(*hvals_voltages)
            h_vals = list(h_vals)
            voltages = list(voltages)
            fig = plt.figure()
            plt.plot(h_vals, voltages, 'o')
            title = f"FMR_" + '{:.2e}'.format(coupling_excitation) + '.png'
            fig.suptitle(title, fontsize=12)
            plt.xlabel("H[mT]")
            plt.ylabel("Vmix[V]")
            plt.savefig(os.path.join(subfolder_name, title + '.png'))
            df = pd.DataFrame.from_dict(
                {
                    'Vmix': voltages,
                    'H': h_vals
                })
            df.to_csv(os.path.join(subfolder_name,
                    f"Coupling_ex_"+ '{:.2e}'.format(coupling_excitation) + ".csv"))
            junction.to_json(os.path.join(subfolder_name, f"junction_base.json"))

if __name__ == "__main__":
    # parameter_scan(junction)
    resonance_curve_scan(junction, "FMR_test")