import numpy as np
import matplotlib.pyplot as plt

from junction import Junction
from constants import Constants
constant = Constants()


def voltage_spin_diode(junction: Junction, start_h, stop_h, step_h):
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
    frequency = 10e9  # 10 Ghz
    voltages = []
    h_vals = np.linspace(start_h, stop_h, 10)
    for h_value in h_vals:
        # set the field
        print(f"Simulation for {h_value}")
        junction.restart()
        junction.set_junction_global_external_field(
            np.array([h_value, 0.0, 0.0]))
        # restart simualtion
        junction.run_simulation(4.5e-9)
        # extract the magnetisation value
        # wait for 5ns
        limited_res = junction.junction_result[
            junction.junction_result['time'] >= 3e-9]
        avg_resistance = np.mean(limited_res['R_free_bottom'])
        omega = 2 * np.pi * frequency
        amplitude = np.sqrt(power / avg_resistance)
        current = amplitude * np.sin(omega * limited_res['time'] +
                                     phase_shift)
        voltage = limited_res['R_free_bottom'] * current
        dc_component = np.mean(voltage)
        voltages.append(dc_component)
        # calcualte magnetoresistance and get the current
    plt.plot(h_vals, voltages)
    plt.show()


def find_resonant_frequency(junction: Junction):
    def step_field(time, step_start=2e-9, step_stop=3e-9):
        Hval = np.zeros((3,))
        if time < step_stop and time > step_start:
            Hval[1] = 10e-3*constant.TtoAm
        return Hval
    junction.set_global_field_function(step_field)
    junction.run_simulation(10e-9)


def frequency_analysis(junction: Junction, time_step=1e-11):
    # send a step pulse to excite the system

    # measure the response
    limited_res = junction.junction_result[
        junction.junction_result['time'] >= 3e-9]
    mx_fft = np.fft.fft(limited_res['m_x_free'], axis=0)
    my_fft = np.fft.fft(limited_res['m_y_free'], axis=0)
    mz_fft = np.fft.fft(limited_res['m_z_free'], axis=0)
    frequency_steps = np.fft.fftfreq(
        mx_fft[0].size, d=time_step)
    max_freq_set = []
    for freq_data in [*mx_fft, *my_fft, *mz_fft]:
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
    j = Junction.from_json('junction.json', persist=True)
    # voltage_spin_diode(j, 0, 400e-3 * constant.TtoAm, 20e-3 * constant.TtoAm)
