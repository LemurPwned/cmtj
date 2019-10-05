# cython: profile=False

# cython: linetrace=False
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False

import math
import json
import time as tm
import pandas as pd
from math import pi, sqrt



MAGNETIC_PERMEABILITY = 12.57e-7
GYRO = 2.21e5
DAMPING = 0.011
TtoAm = 795774.715459
HBAR = 6.62607015e-34/(2*pi)
ELECTRON_CHARGE = 1.60217662e-19

cdef list calculate_tensor_interaction(list m, list tensor, double Ms):
    return [
        -Ms*tensor[0][0]*m[0] -Ms* tensor[0][1]*m[1] -Ms* tensor[0][2]*m[2],
        -Ms*tensor[1][0]*m[0] -Ms* tensor[1][1]*m[1] -Ms* tensor[1][2]*m[2],
        -Ms*tensor[2][0]*m[0] -Ms* tensor[2][1]*m[1] -Ms* tensor[2][2]*m[2]
        ]

cdef list c_cross(list a, list b):
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ]


cdef double c_dot(list a, list b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


class Layer:
    def __init__(self, id_, start_mag, start_anisotropy, K, Ms, coupling, thickness):
        self.id = id_
        # magnetisation options

        self.m = start_mag
        norm = sqrt(self.m[0]**2 + self.m[1]**2 + self.m[2]**2)
        self.m = [self.m[i] / norm for i in range(3)]

        # other parameters section
        self.anisotropy = start_anisotropy
        self.K = K
        self.K_log = K

        self.Ms = Ms
        self.thickness = thickness
        self.coupling = coupling
        self.coupling_log = self.coupling

        self.demagnetisation_tensor = \
            [[0.000, 0, 0], [0, 0.00, 0], [0, 0, 0.93]]

        self.dipole_tensor = \
            [[0.0173, 0.0, 0.0], [0.0, 0.0173, 0.0], [0.0, 0.0, -0.0345]]

        # utility function
        self.params = ['m', 'anisotropy', 'Hext', 'coupling']
        self.Hext = [0.0, 0.0, 0.0]
        self.Hext_const = [0.0, 0.0, 0.0]

        self.update_external_field = None
        self.update_anisotropy = None
        self.update_coupling = None
        self.dmdt = [0., 0., 0.]
        self.log = []
        self.junction_result = None

    def to_dict(self):
        d = {}
        for param in [
                'id', 'anisotropy', 'K', 'm', 'coupling', 'thickness', 'Hext',
                'Hext_const', 'demagnetisation_tensor', 'dipole_tensor', 'Ms'
        ]:
            param_val = getattr(self, param)
            d[param] = param_val
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d['id'], d['m'], d['anisotropy'], d['K'], d['Ms'],
                   d['thickness'], d['coupling'])

    def Heff(self, time, coupled_layers):
        cdef double K_var, coupling_var = 0.0
        if self.update_anisotropy:
            K_var = self.update_anisotropy(time)
        if self.update_external_field:
            self.Hext = self.update_external_field(time)
        if self.update_coupling:
            coupling_var = self.update_coupling(time)

        self.K_log = self.K + K_var
        self.coupling_log = self.coupling + coupling_var

        anis = self.calculate_anisotropy()
        iec = self.calculate_interlayer_exchange_coupling_field(coupled_layers)
        ti_demag = calculate_tensor_interaction(self.m, self.demagnetisation_tensor, self.Ms) 
        ti_dip = calculate_tensor_interaction(self.m, self.dipole_tensor, self.Ms)

        heff = [
            self.Hext_const[i] +\
            self.Hext[i] +\
            anis[i] +\
            ti_demag[i] +\
            ti_dip[i] +\
            iec[i]
            for i in range(3)
        ]
        return heff

    def set_global_external_field_value(self, hval):
        self.Hext_const = hval

    def calculate_anisotropy(self):
        nom = 2 * self.K_log * c_dot(self.m, self.anisotropy) / (MAGNETIC_PERMEABILITY * self.Ms)
        return [nom * anis_v for anis_v in self.anisotropy]

    def calculate_interlayer_exchange_coupling_field(self, coupled_layers_m=None):
        if not coupled_layers_m:
            return [0., 0., 0.]
        return [self.coupling_log*(coupled_layers_m.m[i] - self.m[i]) /(MAGNETIC_PERMEABILITY*self.Ms*self.thickness)
                for i in range(3)]

    def llg(self, time, m, coupled_layers):
        heff = self.Heff(time, coupled_layers)
        prod = c_cross(m, heff)
        prod2 = c_cross(m, prod)
        return [- GYRO * \
            prod[i] - GYRO * \
            DAMPING*prod2[i] for i in range(3)]

    def rk4_step(self, time, time_step, coupled_layers=None):
        m = self.m
        k1 =  self.llg(time, m, coupled_layers)
        k2 =  self.llg(time + 0.5 * time_step, [m[i] + 0.5 * k1[i] for i in range(3)],
                       coupled_layers)
        k3 =  self.llg(time + 0.5 * time_step, [m[i] + 0.5 * k2[i] for i in range(3)],
                       coupled_layers)
        k4 =  self.llg(time + time_step, [m[i] + k3[i] for i in range(3)], coupled_layers)

        m[0] += time_step*(k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0
        m[1] += time_step*(k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0
        m[2] += time_step*(k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0
        norm = sqrt(m[0]**2 + m[1]**2 + m[2]**2)
        # update m
        self.m = [mv/norm for mv in m]


class Junction():
    def __init__(self, _id, layers, couplings, persist=False, save=True):
        self.id = _id
        self.layers = layers
        self.layers_id = {layer.id: layer for layer in layers}
        self.couplings = couplings
        self.log = []
        self.vec_params = ['m', 'Hext', 'Hext_const']
        self.scalar_params = ['K_log', 'K', 'coupling_log', 'coupling']

        self.save = save
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
        Hext = [0., 0., 0.]
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

    def set_global_coupling_function(self, coupling_function):
        for layer in self.layers:
            layer.update_coupling = coupling_function

    def restart(self):
        # just revert to initial parameters
        # that will involve reading the saved file and override the state
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
        write_steps = int(stop_time/1e-9)-1
        save_step = 0 
        for i in range(0, iterations):
            t = i * time_step
            for layer_no, layer in enumerate(self.layers):
                layer.rk4_step(t, time_step, # -1 because id "1" to logical id "0"
                    self.layers[0])
            # calculate magnetoresistance
            if save_step == i:
            # the write time is quite limited
                R = []
                for layer1 in self.layers:
                    for layer2 in self.layers[1:]:
                        if layer1.id != layer2.id:
                            cosTheta = c_dot(layer1.m, layer2.m)
                            R.append(self.magnetoresistance(cosTheta))
                self.log_layer_parameters(t, R)
                save_step += 1

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

