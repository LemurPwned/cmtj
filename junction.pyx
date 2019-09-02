# cython: profile=False

# cython: linetrace=False
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False

import math
import json
import time as tm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt

cimport numpy as np


MAGNETIC_PERMEABILITY = 12.57e-7
GYRO = 2.21e5
DAMPING = 0.01
TtoAm = 795774.715459
HBAR = 6.62607015e-34/(2*pi)
ELECTRON_CHARGE = 1.60217662e-19

cdef np.ndarray calculate_tensor_interaction(m, tensor, Ms):
    return -1. * tensor @ m * Ms

cdef np.ndarray c_cross(np.ndarray[double, ndim=1] a, np.ndarray[double, ndim=1]b):
    return np.array([
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ], dtype=np.double)


cdef double c_dot(np.ndarray[double, ndim=1] a, np.ndarray[double, ndim=1] b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


class Layer:
    def __init__(self, id_, start_mag, start_anisotropy, K, Ms, thickness):
        self.id = id_
        # magnetisation options
        if start_mag is not np.array:
            start_mag = np.array(start_mag, dtype=float)
        self.m = start_mag
        self.m = self.m / np.linalg.norm(self.m)
        self.m_history = []
        if start_anisotropy is not np.array:
            start_anisotropy = np.array(start_anisotropy, dtype=float)

        # other parameters section
        self.anisotropy = start_anisotropy
        self.K = K
        self.K_log = K


        self.Ms = Ms
        self.thickness = thickness
        self.coupling = 1e-5
        self.demagnetisation_tensor = np.array(
            [[0.000, 0, 0], [0, 0.00, 0], [0, 0, 0.93]], dtype=float)

        self.dipole_tensor = np.array(
            [[0.0173, 0.0, 0.0], [0.0, 0.0173, 0.0], [0.0, 0.0, -0.0345]],
            dtype=float)

        # utility function
        self.params = ['m', 'anisotropy', 'Hext', 'coupling']
        self.Hext = np.array([0.0, 0.0, 0.0])
        self.Hext_const = np.array([0.0, 0.0, 0.0])
        self.update_external_field = None
        self.update_anisotropy = None
        self.update_coupling = None
        self.dmdt = np.array([0., 0., 0.])
        self.log = []
        self.junction_result = None

    def to_dict(self):
        d = {}
        for param in [
                'id', 'anisotropy', 'K', 'm', 'coupling', 'thickness', 'Hext',
                'Hext_const', 'demagnetisation_tensor', 'dipole_tensor', 'Ms'
        ]:
            param_val = getattr(self, param)
            if type(param_val) is np.ndarray:
                param_val = param_val.tolist()
            d[param] = param_val
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d['id'], d['m'], d['anisotropy'], d['K'], d['Ms'],
                   d['thickness'])

    def Heff(self, time, coupled_layers):
        cdef double K_var =0.0
        if self.update_anisotropy:
            K_var = self.update_anisotropy(time)
        if self.update_external_field:
            self.Hext = self.update_external_field(time)
        if self.update_coupling:
            self.coupling = self.update_coupling(time)

        self.K_log = self.K + K_var

        heff = \
            self.Hext_const + self.Hext +\
            self.calculate_anisotropy() +\
            calculate_tensor_interaction(self.m, self.demagnetisation_tensor, self.Ms) +\
            self.calculate_interlayer_exchange_coupling_field(coupled_layers) +\
            calculate_tensor_interaction(self.m, self.dipole_tensor, self.Ms)
        return heff

    def set_global_external_field_value(self, hval):
        self.Hext_const = hval

    def calculate_dipole_interaction(self):
        return -1. * self.dipole_tensor @ self.m * self.Ms

    def calculate_demagnetisation_field(self):
        return -1. * self.demagnetisation_tensor @ self.m * self.Ms

    def calculate_anisotropy(self):
        nom = 2 * self.K_log * c_dot(self.m, self.anisotropy) * self.anisotropy
        return nom / (MAGNETIC_PERMEABILITY * self.Ms)

    """
    # this is the old formula that yields non-zero field 
    # for parallel spins
    def calculate_interlayer_exchange_coupling(self, time, coupled_layers):
        heff_iec = np.array([0, 0., 0], dtype=float)
        if not coupled_layers:
            return heff_iec
        for layer in coupled_layers:
            heff_iec += self.coupling*layer.m / \
                (MAGNETIC_PERMEABILITY*self.Ms*self.thickness)
        return heff_iec
    """

    def calculate_interlayer_exchange_coupling_field(self, coupled_layers):
        heff_iec = np.array([0, 0., 0], dtype=float)
        if not coupled_layers:
            return heff_iec
        for layer in coupled_layers:
            heff_iec += self.coupling*(layer.m - self.m) / \
                ( MAGNETIC_PERMEABILITY*self.Ms*self.thickness)
        return heff_iec

    def llg(self, time, m, coupled_layers):
        heff = self.Heff(time, coupled_layers)
        prod = c_cross(m, heff)
        dmdt = - GYRO * \
            prod - GYRO * \
             DAMPING*c_cross(m, prod)
        return dmdt

    def llg_STT(self, time, m, coupled_layers, spin_polarized_layer):
        """
        :param spin_polarized_layer is a vector of a layer in which 
            the current is polarised
        """
        heff = self.Heff(time, coupled_layers)
        aj = ( HBAR * self.current_density *
              self.spin_polarisation_efficiency) / (
                  2 *  ELECTRON_CHARGE *
                   MAGNETIC_PERMEABILITY * self.Ms * self.thickness)
        bj = np.array([0., 0., 0.])
        p_cross = c_cross(m, spin_polarized_layer)
        stt_term =  GYRO * (aj * c_cross(m, p_cross) + bj * p_cross)

        dmdt = -1.0* GYRO * \
            c_cross(m, heff) - 1.0* GYRO * \
             DAMPING*c_cross(m, c_cross(m, heff)) +\
            stt_term
        return dmdt

    def rk4_step(self, time, time_step, coupled_layers=None):
        cdef np.ndarray[double, ndim=1] k1, k2, k3, k4, m
        m = self.m
        k1 = time_step * self.llg(time, m, coupled_layers)
        k2 = time_step * self.llg(time + 0.5 * time_step, m + 0.5 * k1,
                                  coupled_layers)
        k3 = time_step * self.llg(time + 0.5 * time_step, m + 0.5 * k2,
                                  coupled_layers)
        k4 = time_step * self.llg(time + time_step, m + k3, coupled_layers)
        m = m + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        m = m / sqrt(m[0]**2 + m[1]**2 + m[2]**2)
        # update m
        self.m = m


class Junction():
    def __init__(self, _id, layers, couplings, persist=False, save=True):
        self.id = _id
        self.layers = layers
        self.layers_id = {layer.id: layer for layer in layers}
        self.couplings = couplings
        self.log = []
        self.vec_params = ['m', 'Hext', 'Hext_const']
        self.scalar_params = ['K_log', 'K']

        self.save = False
        self.persist = persist  # read file only or memory leftover?

        self.R_labs = []
        self.Rp = 100   
        self.Rap = 200
        self.avg_RpRap = 0.5 * (self.Rap - self.Rp)
        self.snapshot = self.to_dict()

    def overwrite_snapshot(self):
        self.snapshot = self.to_dict()

    def to_json(self, filename):
        json.dump(self.to_dict(), open(filename, 'w'), indent=4)

    @classmethod
    def from_json(cls, filename, persist=False):
        d = json.load(open(filename, 'r'))
        return Junction.from_dict(d, persist)

    @classmethod
    def from_dict(cls, d, persist=False):
        layers = [Layer.from_dict(layer) for layer in d['layers']]
        return cls(d['id'], layers, d['couplings'], persist=persist)

    def to_dict(self):
        d = {}
        layers_dict = [layer.to_dict() for layer in self.layers]
        d['layers'] = layers_dict
        for param in ['id', 'couplings', 'Rap', 'Rp']:
            d[param] = getattr(self, param)
        return d

    def magnetoresistance(self, costheta):
        """
        :param theta
            the angle between two magnetisation layers
        """
        return self.Rp + self.avg_RpRap * (1. - costheta)

    def log_layer_parameters(self, time, R):
        value_list = [time , *R]  # nanoseconds
        for param in self.scalar_params:
            for layer in self.layers:
                pval = getattr(layer, param)
                value_list.append(pval)
        for param in self.vec_params:
            for layer in self.layers:
                pval = getattr(layer, param)
                value_list.extend(pval)
        self.log.append(value_list)

    def set_junction_global_external_field(self, constant_field_value, axis='x'):
        Hext = np.zeros((3,))
        if axis == 'x':
            Hext[0] = constant_field_value
        elif axis == 'y':
            Hext[1] = constant_field_value
        elif axis == 'z':
            Hext[2] = constant_field_value
        for layer in self.layers:
            layer.Hext_const = Hext

    def set_global_field_function(self, field_function):
        for layer in self.layers:
            layer.update_external_field = field_function

    def set_global_anisotropy_function(self, anisotropy_function):
        for layer in self.layers:
            layer.update_anisotropy = anisotropy_function

    def restart(self):
        # just revert to initial parameters
        # that will involve reading the saved file and overrite the state
        # restore the snapshot
        self.R_labs = []
        self.log = []
        for key in self.snapshot:
            if key != 'layers':
                setattr(self, key, self.snapshot[key])
        self.layers = [
            Layer.from_dict(layer_dict)
            for layer_dict in self.snapshot['layers']
        ]

    def run_simulation(self, stop_time, time_step=1e-13):
        # LLG#
        cdef double tstart, tend, t
        cdef int iterations, i
        if not self.layers:
            raise ValueError("No layers in junction!")
        tstart = tm.time()
        iterations = int(stop_time / time_step)

        # just for labels
        for layer1 in self.layers:
            for layer2 in self.layers[1:]:
                if layer1.id != layer2.id:
                    self.R_labs.append(f'R_{layer1.id}_{layer2.id}')
        # start the simulation
        for i in range(0, iterations):
            t = i * time_step
            for layer in self.layers:
                layer.rk4_step(t, time_step)
            # calculate magnetoresistance
            R = []
            for layer1 in self.layers:
                for layer2 in self.layers[1:]:
                    if layer1.id != layer2.id:
                        cosTheta = c_dot(layer1.m, layer2.m)
                        R.append(self.magnetoresistance(cosTheta))
            self.log_layer_parameters(t, R)

        tend = tm.time()
        print(
            f"Simulation has finished in {tend-tstart}s. Writting the results..."
        )
        cols = []
        for param in self.scalar_params:
            for layer in self.layers:
                cols.append(param + '_' + layer.id)
        for param in self.vec_params:
            for layer in self.layers:
                for suffix in ['x', 'y', 'z']:
                    cols.append(param + '_' + suffix + '_' + layer.id)
        if self.persist:
            self.junction_result = pd.DataFrame(
                data=self.log, columns=['time', *self.R_labs,
                                        *cols])
            if self.save:
                self.junction_result.to_csv('results.csv', index=False)
        else:
            if self.save:
                pd.DataFrame(data=self.log, columns=['time', *self.R_labs,
                                                 *cols]).to_csv('results.csv',
                                                 index=False)


def plot_results():
    df = pd.read_csv('results.csv')
    # df[['m_x_free', 'm_y_free', 'm_z_free']].plot()
    df['R_free_bottom'].plot()
    plt.show()
