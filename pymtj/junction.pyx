# cython: profile=False

# cython: linetrace=False
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False

import math
import json
import time as tm
import os 
import pandas as pd
from math import pi, sqrt

MAGNETIC_PERMEABILITY = 12.57e-7
GYRO = 221000.0
DAMPING = 0.011
TtoAm = 795774.715459
HBAR = 6.62607015e-34 / (2 * pi)
ELECTRON_CHARGE = 1.60217662e-19

cdef list calculate_tensor_interaction(list m, list tensor, float Ms):
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
    def __init__(self, id_, start_mag, anisotropy, K, Ms, coupling, thickness,
                 demag_tensor, dipole_tensor):
        self.id = id_
        # magnetisation options
        norm = sqrt(start_mag[0]**2 + start_mag[1]**2 + start_mag[2]**2)
        self.m = [start_mag[i] / norm for i in range(3)]
        self.start_mag = start_mag
        # other parameters section
        self.anisotropy = anisotropy
        self.K = K
        self.K_log = self.K

        self.Ms = Ms
        self.thickness = thickness
        self.coupling = coupling
        self.coupling_log = self.coupling

        # tensors
        self.demagnetisation_tensor = demag_tensor
        self.dipole_tensor = dipole_tensor

        # utility function
        self.params = ['m', 'anisotropy', 'Hext', 'coupling']
        self.Hext = [0.0, 0.0, 0.0]
        self.Hext_const = [0.0, 0.0, 0.0]

        # updating functions
        self.update_external_field = None
        self.update_anisotropy = None
        self.update_coupling = None

        self.junction_result = None

    def to_dict(self):
        d = {}
        for param in [
                'id', 'anisotropy', 'K', 'm', 'coupling', 'thickness', 'Hext',
                'Hext_const', 'demagnetisation_tensor', 'dipole_tensor', 'Ms',
                'start_mag'
        ]:
            param_val = getattr(self, param)
            d[param] = param_val
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d['id'],
                   start_mag=d['start_mag'],
                   anisotropy=d['anisotropy'],
                   K=d['K'],
                   Ms=d['Ms'],
                   coupling=d['coupling'],
                   thickness=d['thickness'],
                   dipole_tensor=d['dipole_tensor'],
                   demag_tensor=d['demagnetisation_tensor'])

    def Heff(self, m, time, coupled_layers, Ms):
        cdef float K_var, coupling_var = 0.0
        if self.update_anisotropy:
            K_var = self.update_anisotropy(time)
        if self.update_external_field:
            self.Hext = self.update_external_field(time)
        if self.update_coupling:
            coupling_var = self.update_coupling(time)

        self.K_log = self.K + K_var
        self.coupling_log = self.coupling + coupling_var

        anis = self.calculate_anisotropy(m)
        iec = self.calculate_interlayer_exchange_coupling_field(
            m, coupled_layers)
        ti_demag = calculate_tensor_interaction(coupled_layers, 
                                                self.demagnetisation_tensor,
                                                Ms)
        ti_dip = calculate_tensor_interaction(m, self.dipole_tensor, self.Ms)

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
        """
        :param hval
            Value of the field in A/m
        Units:
            returns the value in A/m
        """
        self.Hext_const = hval

    def calculate_anisotropy(self, m):
        """
        :param m
            magnetisation value A/m, normalized
        :returns
            anisotropy field value
        Units:
            returns the value in A/m
        """
        const = (2 * self.K_log) / (MAGNETIC_PERMEABILITY * self.Ms)
        nom = const * c_dot(
            m, self.anisotropy)
        return [nom * anis_v for anis_v in self.anisotropy]

    def calculate_interlayer_exchange_coupling_field(self,
                                                     m,
                                                     coupled_layers_m=None):
        """
        :param m
            magnetisation value A/m, normalized
        :param coupled_layers_m
            the coupled layer's magnetisation A/m, normalized
        :returns
            IEC field value
        Units:
            returns the value in A/m
        """
        if not coupled_layers_m:
            return [0., 0., 0.]
        const = self.coupling_log / (MAGNETIC_PERMEABILITY * self.Ms *
                                     self.thickness)
        return [const * (coupled_layers_m[i] - m[i]) for i in range(3)]

    def llg(self, heff, time, m, coupled_layers):
        """
        Calculate LLG section in form
        dmdt = -gamma*(m x heff) - gamma*damping*(m x (m x heff))
        :param time
            time of the simulation
        :param m 
            magnetisation value A/m, normalized
        :param coupled_layers_m
            the coupled layer's magnetisation A/m, normalized
        """
        # normalize -- needed for rk4 equations
        norm = sqrt(m[0]**2 + m[1]**2 + m[2]**2)
        m[0] = m[0] / norm
        m[1] = m[1] / norm
        m[2] = m[2] / norm
        # heff = self.Heff(m, time, coupled_layers)
        prod = c_cross(m, heff)
        prod2 = c_cross(m, prod)
        dmdt = [-GYRO * prod[i] - GYRO * DAMPING * prod2[i] for i in range(3)]
        return dmdt

    def rk4_step(self, time, time_step, coupled_layers, Ms):
        """
        Solves RK4 for the LLG equation
        :param time
            time of the simulation
        :param m 
            magnetisation value A/m, normalized
        :param coupled_layers_m
            the coupled layer's magnetisation A/m, normalized
        """
        m = self.m
        heff = self.Heff(m, time, coupled_layers, Ms)
        # heff=0
        k1 = self.llg(heff, time, m, coupled_layers)
        k1 = [time_step * k1[i] for i in range(3)]

        k2 = self.llg(heff, time + (0.5 * time_step),
                      [m[i] + (0.5 * k1[i]) for i in range(3)], coupled_layers)
        k2 = [time_step * k2[i] for i in range(3)]

        k3 = self.llg(heff, time + (0.5 * time_step),
                      [m[i] + (0.5 * k2[i]) for i in range(3)], coupled_layers)
        k3 = [time_step * k3[i] for i in range(3)]

        k4 = self.llg(heff, time + time_step, [m[i] + k3[i] for i in range(3)],
                      coupled_layers)
        k4 = [time_step * k4[i] for i in range(3)]

        m[0] += ((k1[0] + (2.0 * k2[0]) + (2.0 * k3[0]) + k4[0]) / 6.0)
        m[1] += ((k1[1] + (2.0 * k2[1]) + (2.0 * k3[1]) + k4[1]) / 6.0)
        m[2] += ((k1[2] + (2.0 * k2[2]) + (2.0 * k3[2]) + k4[2]) / 6.0)
        norm = sqrt(m[0]**2 + m[1]**2 + m[2]**2)
        # # update m
        m[0] = m[0] / norm
        m[1] = m[1] / norm
        m[2] = m[2] / norm
        self.m = m


class Junction():
    def __init__(self,
                 _id,
                 layers,
                 couplings,
                 persist=True,
                 save=False,
                 save_location='results.csv'):
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
        self.save_location = save_location
        self.avg_RpRap = 0.5 * (self.Rap - self.Rp)
        self.snapshot = self.to_dict()

    def overwrite_snapshot(self):
        self.snapshot = self.to_dict()

    def to_json(self, filename):
        json.dump(self.to_dict(), open(filename, 'w'), indent=4)

    @classmethod
    def from_json(cls, filename, persist=True):
        d = json.load(open(filename, 'r'))
        return Junction.from_dict(d, persist)

    @classmethod
    def from_dict(cls, d, persist=True, save=False):
        layers = [Layer.from_dict(layer) for layer in d['layers']]
        return cls(d['id'], layers, d['couplings'], persist=persist, save=save)

    def to_dict(self):
        d = {}
        layers_dict = [layer.to_dict() for layer in self.layers]
        d['layers'] = layers_dict
        for param in ['id', 'couplings', 'Rap', 'Rp']:
            d[param] = getattr(self, param)
        return d

    def magnetoresistance(self, costheta):
        """
        :param costheta
            the cosine(angle) between two magnetisation layers
        """
        return self.Rp + (self.avg_RpRap * (1. - costheta))

    def log_layer_parameters(self, time, R):
        value_list = [time, *R]  # nanoseconds
        for param in self.scalar_params:
            for layer in self.layers:
                pval = getattr(layer, param)
                value_list.append(pval)
        for param in self.vec_params:
            for layer in self.layers:
                pval = getattr(layer, param)
                value_list.extend(pval)
        self.log.append(value_list)

    def set_junction_global_external_field(self,
                                           constant_field_value,
                                           axis='x'):
        """
        Set the external magnetic field for all layers in the junction
        Can only be set along one axis
        :param constant_field_value
            value of the constant field in a given axis
        :param axis
            x, y or z axis of the field (positive)
        """
        Hext = [0., 0., 0.]
        if axis == 'x':
            Hext[0] = constant_field_value
        elif axis == 'y':
            Hext[1] = constant_field_value
        elif axis == 'z':
            Hext[2] = constant_field_value
        for layer in self.layers:
            layer.Hext_const = Hext


    def set_global_coupling(self, coupling_val):
        """
        Sets coupling to the coupling val for all layers
        in the junction
        :param coupling_val
            Coupling value, unit [J/m^2]
        """
        for layer in self.layers:
            layer.coupling = coupling_val

    def set_layer_coupling(self, layer_id, coupling_val):
        """
        Sets coupling to the coupling val for a specific layer
        in the junction
        :param layer_id
            layer id, for which the value will be set
        :param coupling_val
            Coupling value, unit [J/m^2]
        """
        found = False 
        for layer in self.layers:
            if layer.id == layer_id:
                found = True 
                layer.coupling = coupling_val
                return
        if not found:
            raise AttributeError(f"Layer with id: {layer_id} not found!")

    def set_layer_K(self, layer_id, K):
        """
        Sets K (anisotropy) to the coupling val for a specific layer
        in the junction
        :param layer_id
            layer id, for which the value will be set
        :param coupling_val
            Coupling value, unit [J/m^3]
        """
        found = False 
        for layer in self.layers:
            if layer.id == layer_id:
                found = True 
                layer.K = K
                return 
        if not found:
            raise AttributeError(f"Layer with id: {layer_id} not found!")


    def set_global_field_function(self, field_function):
        """
        Updates the external field value at each step
        :param field_function
            function that updates the external field at each
            time step
        Units:
            returns the value in A/m
        """
        for layer in self.layers:
            layer.update_external_field = field_function

    def set_global_anisotropy_function(self, anisotropy_function):
        """
        Updates the anisotropy at each step
        :param field_function
            function that updates the external field at each
            time step
        Units:
            returns the value in A/m
        """
        for layer in self.layers:
            layer.update_anisotropy = anisotropy_function

    def set_layer_anisotropy_function(self, layer_id, anisotropy_function):
        """
        Updates the anisotropy for a given layer at each step
        :param field_function
            function that updates the external field at each
            time step
        :param layer_id
            layer id, for which the value will be set
        Units:
            returns the value in A/m
        """
        found = False 
        for layer in self.layers:
            if layer.id == layer_id:
                found = True 
                layer.update_anisotropy = anisotropy_function
                return 
        if not found:
            raise AttributeError(f"Layer with id: {layer_id} not found!")

    def set_global_coupling_function(self, coupling_function):
        """
        Updates the coupling value at each step
        :param field_function
            function that updates the coupling at each
            time step
        Units:
            returns the value in A/m
        """
        for layer in self.layers:
            layer.update_coupling = coupling_function

    def soft_log_reset(self):
        self.log = []

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
        cdef double t
        cdef int iterations, i, STEP, ITERATIONS
        if not self.layers:
            raise ValueError("No layers in junction!")
        tstart = os.times().elapsed
        iterations = int(stop_time / time_step)
        # just for labels
        if not self.R_labs:
            for layer1 in self.layers:
                for layer2 in self.layers[1:]:
                    if layer1.id != layer2.id:
                        self.R_labs.append(f'R_{layer1.id}_{layer2.id}')

        # start the simulation
        save_step = 0
        write_steps = int(stop_time/(0.01e-9))
        mod = int(iterations/write_steps)
        STEP_TIME = 0.01*1e-9
        ITERATIONS = int(stop_time / time_step)
        STEP = int(STEP_TIME/ time_step)
        writeEvery = STEP-1; 
        # print(STEP, writeEvery, ITERATIONS)
        for i in range(0, ITERATIONS):
            t = i * time_step
            for layer_no, layer in enumerate(self.layers):
                other_m = self.layers[self.couplings[layer_no][0] - 1].m
                other_Ms = self.layers[self.couplings[layer_no][0] - 1].Ms
                layer.rk4_step(
                    t,
                    time_step,  # -1 because id "1" to logical id "0"
                    other_m,
                    other_Ms)
            # calculate magnetoresistance
            # the write time is quite limited
            if i == writeEvery:
                R = []
                for layer1 in self.layers:
                    for layer2 in self.layers[1:]:
                        if layer1.id != layer2.id:
                            cosTheta = c_dot(layer1.m, layer2.m)
                            # norm1 = sqrt(layer1.m[0]**2 + layer1.m[1]**2 + layer1.m[2]**2)
                            # norm2 = sqrt(layer2.m[0]**2 + layer2.m[1]**2 + layer2.m[2]**2)
                            # cosTheta = cosTheta/(norm1*norm2)
                            R.append(self.magnetoresistance(cosTheta))
                writeEvery += STEP
                self.log_layer_parameters(t, R)

        tend = os.times().elapsed
        print(
            f"{iterations} finished in {tend-tstart}s. Writting the results..."
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
                data=self.log, columns=['time', *self.R_labs, *cols])
            if self.save:
                self.junction_result.to_csv(self.save_location, index=False)
        else:
            if self.save:
                pd.DataFrame(data=self.log,
                             columns=['time', *self.R_labs,
                                      *cols]).to_csv(self.save_location,
                                                     index=False)
