import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from constants import Constants

constant = Constants()


class Layer:
    def __init__(self, id_, start_mag,
                 start_anisotropy, Ms, thickness):
        self.id = id_
        if start_mag is not np.array:
            start_mag = np.array(start_mag, dtype=float)
        self.m = start_mag
        self.m = self.m/np.linalg.norm(self.m)
        self.m_history = []
        if start_anisotropy is not np.array:
            start_anisotropy = np.array(start_anisotropy, dtype=float)
        self.anisotropy = start_anisotropy
        self.Ms = Ms
        self.thickness = thickness
        self.coupling = 1e-5
        self.demagnetisation_tensor = np.array(
            [[0.0005, 0, 0], [0, 0.005, 0], [0, 0, 0.93]], dtype=float)

        self.dipole_tensor = np.array(
            [[0.0173, 0.0, 0.0], [0.0, 0.0173, 0.0], [0.0, 0.0, -0.0345]], dtype=float)

        self.params = ['m', 'anisotropy', 'Hext', 'coupling']
        self.Hext = np.array([0.0, 0.0, 0.0])
        self.update_external_field = None
        self.update_anisotropy = None
        self.update_coupling = None
        self.log = []

    def Heff(self, time, coupled_layers):
        if self.update_anisotropy:
            self.anisotropy = self.update_anisotropy(time)
        if self.update_external_field:
            self.Hext = self.update_external_field(self.Hext, time)
        if self.update_coupling:
            self.coupling = self.update_coupling(time)

        heff = \
            self.Hext + \
            self.calculate_demagnetisation_field() + \
            self.calculate_anisotropy(time) + \
            self.calculate_interlayer_exchange_coupling(time, coupled_layers)
        return heff

    def calculate_zeeman(self, time):
        return -1 * constant.MAGNETIC_PERMEABILITY*np.dot(self.m, self.Hext)

    def calculate_dipole_interaction(self, time):
        return self.dipole_tensor@np.array([0, 0, 1.])

    def calculate_anisotropy(self, time):
        nom = 2*np.dot(self.m, self.anisotropy)*self.anisotropy
        return nom/(constant.MAGNETIC_PERMEABILITY*self.Ms)

    def calculate_interlayer_exchange_coupling(self, time, coupled_layers):
        heff_iec = np.array([0, 0, 0])
        if not coupled_layers:
            return heff_iec
        for layer in coupled_layers:
            heff_iec += self.coupling*layer.m / \
                (constant.MAGNETIC_PERMEABILITY*self.Ms*self.thickness)
        return heff_iec

    def calculate_demagnetisation_field(self):
        return -1.*self.demagnetisation_tensor@self.m*self.Ms

    def log_layer_parameters(self, time):
        value_list = [time/1e-9]  # nanoseconds
        for param in self.params:
            pval = getattr(self, param)
            if type(pval) is np.ndarray:
                value_list.extend(pval)
            else:
                value_list.append(getattr(self, param))
        self.log.append(value_list)

    def llg(self, time, m, coupled_layers):
        heff = self.Heff(time, coupled_layers)
        dmdt = -1*constant.GYRO * \
            np.cross(m, heff) - 1*constant.GYRO * \
            constant.DAMPING*np.cross(m, np.cross(m, heff))
        return dmdt

    def run_simulation(self, stop_time, dt=1e-12, time_step=1e-13):
        # LLG#
        for t in np.arange(0, stop_time, dt):
            time = t
            stopping_cond = (t + dt - time_step)
            while time < stopping_cond:
                self.rk4_step(time, time_step)
                time += time_step
            self.log_layer_parameters(t)
        pd.DataFrame(data=self.log, columns=['time', 'mx', 'my', 'mz',
                                             'ax', 'ay', 'az',
                                             'Hext_x', 'Hext_y', 'Hext_z']).to_csv('results.csv')

    def rk4_step(self, time, time_step, coupled_layers=None):
        m = self.m
        k1 = time_step*self.llg(time, self.m, coupled_layers)
        k2 = time_step*self.llg(time + 0.5*time_step,
                                self.m+0.5*k1, coupled_layers)
        k3 = time_step*self.llg(time + 0.5*time_step,
                                self.m+0.5*k1, coupled_layers)
        k4 = time_step*self.llg(time + time_step, self.m+k3, coupled_layers)
        m = m + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        m = m/np.linalg.norm(m)
        # update m
        self.m_history.append(self.m)
        self.m = m


def plot_results():
    df = pd.read_csv('results.csv')
    # df[['mx', 'my', 'mz']].plot()
    df[['mx', 'my', 'mz']].plot()
    plt.show()


class Junction():
    def __init__(self):
        self.layers = []
        self.couplings = []
        self.log = []
        self.params = ['m', 'anisotropy', 'Hext', 'coupling']

    def log_layer_parameters(self, time):
        value_list = [time/1e-9]  # nanoseconds
        for param in self.params:
            for layer in self.layers:
                pval = getattr(layer, param)
                value_list.extend(pval)
        self.log.append(value_list)

    def run_simulation(self, stop_time, dt=1e-12, time_step=1e-13):
        # LLG#
        if not self.layers:
            raise ValueError("No layers in junction!")
        for t in np.arange(0, stop_time, dt):
            time = t
            stopping_cond = (t + dt - time_step)
            while time < stopping_cond:
                for i, layer in enumerate(self.layers):
                    print(i, layer.id)
                    coupled_layers = [self.layers[j]
                                      for j in self.couplings[i]]
                    layer.rk4_step(time, time_step, coupled_layers)
            self.log_layer_parameters(t)

        cols = []
        for param in self.params:
            for layer in self.layers:
                for suffix in ['x', 'y', 'z']:
                    cols.append(param+'_'+suffix+'_'+layer.id)
        pd.DataFrame(data=self.log, columns=[
                     'time', *cols]).to_csv('results.csv')


if __name__ == "__main__":

    junction = Junction()
    l1 = Layer(
        "free", [0, 0.5, 0.5], [0, 100e6, 0], 1200e3, 3e-9
    )
    l2 = Layer(
        "reference", [0, 0., 1.], [0, 0, 0], 1200e3, 3e-9
    )
    juncton.layers = [l1, l2]
    junction.couplings = [[1], [0]]
    l2.update_external_field = lambda Hext, t: np.array(
        [0, 20*np.sin(5e9*2*np.pi*t), 0.])
    l1.update_external_field = lambda Hext, t: np.array(
        [0, 20*np.sin(5e9*2*np.pi*t), 0.])
    junction.run_simulation(1e-9)
    # plot_results()
