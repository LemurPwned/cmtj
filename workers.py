import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
from itertools import repeat


from junction import Junction
from constants import Constants
constant = Constants()


def step_field(time, step_start=5e-9, step_stop=5.01e-9):
    Hval = np.zeros((3,))
    if time < step_stop and time > step_start:
        Hval[0] = 10e-3*constant.TtoAm
    return Hval


def calculate_single_voltage(h_value, junction: Junction, frequency):
    # set the field
    phase_shift = 0
    power = 10e-6
    omega = 2 * np.pi * frequency
    print(f"Simulation for {h_value}")
    junction.restart()
    junction.set_junction_global_external_field(
        h_value*constant.TtoAm, axis='x')
    junction.set_global_anisotropy_function(anisotropy_update)
    # restart simualtion
    junction.run_simulation(9e-9)
    # extract the magnetisation value
    # wait for 5ns
    limited_res = junction.junction_result[
        junction.junction_result['time'] >= 5e-9]
    avg_resistance = np.mean(limited_res['R_free_bottom'])
    print(f"Avg resistance {avg_resistance}")
    amplitude = np.sqrt(power / avg_resistance)
    current = amplitude * np.sin(omega * limited_res['time']*1e-9 +
                                    phase_shift)
    voltage = limited_res['R_free_bottom'] * current

    return np.mean(voltage)


def anisotropy_update(time):
    frequency = 5.4e9  # 10 Ghz
    omega = 2 * np.pi * frequency
    anisotropy = np.zeros((3,))
    anisotropy[0] = 100e3*np.sin(2*omega*time)
    return anisotropy

def voltage_spin_diode(junction: Junction, start_h, stop_h, multiprocess=True):
    """
    Field scan
            - scan with the field in range [min_field, max_field]
            - collect the magnetoresistance
            - calculate average magnetoresistance assuming some current fed into the circuit
            - calculate constant DC voltage
            - aggregate
    """
    phase_shift = 0
    power = 10e-6
    frequency = 5.4e9  # 10 Ghz
    omega = 2 * np.pi * frequency
    voltages = []

    h_vals = np.linspace(start_h, stop_h, 30)
    if multiprocess:
        junction.set_global_anisotropy_function(anisotropy_update)
        with Pool() as pool:
            voltages = pool.starmap(
                calculate_single_voltage,
                zip(h_vals, repeat(junction), repeat(frequency))
            )
    else:
        for h_value in h_vals:
            # set the field
            print(f"Simulation for {h_value}")
            junction.restart()
            junction.set_junction_global_external_field(
                h_value*constant.TtoAm, axis='x')
            junction.set_global_anisotropy_function(anisotropy_update)
            # restart simualtion
            junction.run_simulation(9e-9)
            # extract the magnetisation value
            # wait for 5ns
            limited_res = junction.junction_result[
                junction.junction_result['time'] >= 5e-9]
            avg_resistance = np.mean(limited_res['R_free_bottom'])
            print(f"Avg resistance {avg_resistance}")
            amplitude = np.sqrt(power / avg_resistance)
            current = amplitude * np.sin(omega * limited_res['time']*1e-9 +
                                        phase_shift)
            voltage = limited_res['R_free_bottom'] * current
            dc_component = np.mean(voltage)
            voltages.append(dc_component)
        # calcualte magnetoresistance and get the current
    # plt.plot(h_vals, voltages)
    # plt.show()
    df = pd.DataFrame.from_dict({'H': h_vals, 'Vmix': voltages})
    df.to_csv('voltage_spin-diode.csv')


def find_resonant_frequency(junction: Junction):
    def step_field(time, step_start=2e-9, step_stop=3e-9):
        Hval = np.zeros((3,))
        if time < step_stop and time > step_start:
            Hval[1] = 10e-3*constant.TtoAm
        return Hval
    junction.set_global_field_function(step_field)
    junction.run_simulation(10e-9)


def frequency_analysis(junction: Junction, time_step=1e-13):
    """
    Returns the values of frequency analysis
    For each magnetisation axis the function returns the values of 
    1. maximum-amplitude frequency (resonant frequency)
    2. the amplitude of resonant frequency
    :param results
        filepath to csv file
    :param time_step
        time step for fft Fourier steps
    """
    # send a step pulse to excite the system
    print(
        f"Calculating the resonant frequencies for the system..., step size {time_step}")    # measure the response
    limited_res = junction.junction_result[
        junction.junction_result['time'] >= 5.1e-9]
    mx_fft = np.fft.fft(limited_res['m_x_free'], axis=0)
    my_fft = np.fft.fft(limited_res['m_y_free'], axis=0)
    mz_fft = np.fft.fft(limited_res['m_z_free'], axis=0)
    frequency_steps = np.fft.fftfreq(
        mx_fft.size, d=time_step)
    # print(frequency_steps)
    max_freq_set = []
    for freq_data in [mx_fft, my_fft, mz_fft]:
        freq_data = abs(freq_data)
        max_val = 0
        max_freq = 0
        for frequency, amp in zip(frequency_steps, freq_data):
            if np.abs(amp) > max_val and frequency > 0:
                max_val = amp
                max_freq = frequency
        max_freq_set.append([max_freq / 1e9, max_val])
    # display Fourier
    return np.array(max_freq_set, dtype=float)


def frequency_analysis_csv(results, time_step=1e-13):
    """
    Returns the values of frequency analysis
    For each magnetisation axis the function returns the values of 
    1. maximum-amplitude frequency (resonant frequency)
    2. the amplitude of resonant frequency
    :param results
        filepath to csv file
    :param time_step
        time step for fft Fourier steps
    """
    # send a step pulse to excite the system
    df = pd.read_csv(results)
    print(
        f"Calculating the resonant frequencies for the system..., step size {time_step}")
    # measure the response
    limited_res = df[df['time'] >= 5.1e-9]
    mx_fft = np.fft.fft(limited_res['m_x_free'], axis=0)
    my_fft = np.fft.fft(limited_res['m_y_free'], axis=0)
    mz_fft = np.fft.fft(limited_res['m_z_free'], axis=0)
    frequency_steps = np.fft.fftfreq(
        mx_fft.size, d=time_step)
    # print(frequency_steps)
    max_freq_set = []
    for freq_data in [mx_fft, my_fft, mz_fft]:
        freq_data = abs(freq_data)
        max_val = 0
        max_freq = 0
        for frequency, amp in zip(frequency_steps, freq_data):
            if np.abs(amp) > max_val and frequency > 0:
                max_val = amp
                max_freq = frequency
        max_freq_set.append([max_freq / 1e9, max_val])
    # display Fourier
    return np.array(max_freq_set, dtype=float)


if __name__ == "__main__":
    # j = Junction.from_json('junction.json', persist=True)
    voltage_spin_diode(j, 0, 400e-3 * constant.TtoAm, 20e-3 * constant.TtoAm)
    print(frequency_analysis_csv('results.csv'))
