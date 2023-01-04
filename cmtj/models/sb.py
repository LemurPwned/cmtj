import math
import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import List

import numpy as np
from tqdm import tqdm

from ..utils import gamma_rad, mu0


@dataclass
class VectorObj:
    """Vector object for standard manipulation.
    :param theta: positive z-axis angle (in xz plane) in radians.
    :param phi: positive x-axis (in xy plane) angle in radians
    :param mag: magnitude of the vector, if not set defaults to 1 *unit vector*
    """
    theta: float  # in radians
    phi: float  # rad
    mag: float = 1

    def get_cartesian(self):
        """Returns the vector in Cartesian coordinates with (x, y, z) compnents"""
        return VectorObj.from_spherical(self.theta, self.phi, self.mag)

    @staticmethod
    def from_spherical(theta, phi, mag=1):
        """Creates a Cartesian vector from spherical components"""
        return [
            mag * math.sin(theta) * math.cos(phi),
            mag * math.sin(theta) * math.sin(phi), mag * math.cos(theta)
        ]


def box_muller_random(mean, std):
    """
    Generates Gaussian noise with mean and standard deviation
    using the Box-Muller transform.
    https://en.wikipedia.org/wiki/Box–Muller_transform
    :param mean: mean of the Gaussian.
    :param std: standard deviation of the Gaussian.
    """
    u1 = np.random.uniform(0, 1)
    u2 = np.random.uniform(0, 1)
    mag = std * math.sqrt(-2.0 * math.log(u1))
    z0 = mag * math.cos(2 * math.pi * u2) + mean
    z1 = mag * math.sin(2 * math.pi * u2) + mean
    return z0, z1


@dataclass
class LayerSB:
    """Basic Layer for Smit-Beljers model.
    :param thickness: thickness of the FM layer (effective).
    :param m: initial magnetisation vector (Cauchy condition).
    :param Kv: volumetric (in-plane) anisotropy. Only phi and mag count [J/m^3].
    :param Ks: surface anisotropy (out-of plane, or perpendicular) value [J/m^3].
    :param Ms: magnetisation saturation value in [A/m].
    :param thermal_noise: if !=0 then introduces a small random disturbance to init m.
    """
    thickness: float
    m: VectorObj
    Kv: VectorObj
    Ks: float
    Ms: float
    thermal_noise: float = 0

    def __post_init__(self):
        if self.thermal_noise:
            self.m = self.add_thermal_noise(self.m, self.thermal_noise)

    def add_thermal_noise(self, m, thermal) -> VectorObj:
        """
        Adds small themal noise to the magnetization vector
        """
        theta, phi = m.theta, m.phi
        z0, z1 = box_muller_random(0, thermal)
        theta += z0
        phi += z1
        return VectorObj(theta, phi)

    @property
    def stheta(self):
        return math.sin(self.m.theta)

    @property
    def ctheta(self):
        return math.cos(self.m.theta)

    @property
    def sphi(self):
        return math.sin(self.m.phi)

    @property
    def cphi(self):
        return math.cos(self.m.phi)

    def ext_field(self, Hinplane: VectorObj):
        hx, hy, hz = Hinplane.get_cartesian()
        return -self.Ms * mu0 * (hx * self.stheta * self.cphi + hy *
                                 self.sphi * self.ctheta + hz * self.ctheta)

    def grad_ext_field(self, Hinplane: VectorObj):
        hx, hy, hz = Hinplane.get_cartesian()
        pref = -self.Ms * mu0
        dEdtheta = pref * (hx * self.cphi * self.ctheta +
                           hy * self.sphi * self.ctheta - hz * self.stheta)

        dEdphi = pref * (-hx * self.sphi * self.stheta +
                         hy * self.stheta * self.cphi)

        d2Edtheta2 = pref * (-hx * self.stheta * self.cphi -
                             hy * self.sphi * self.stheta - hz * self.ctheta)
        d2Edphi2 = pref * (-hx * self.stheta * self.cphi -
                           hy * self.sphi * self.stheta)
        d2Edphidtheta = pref * (-hx * self.sphi * self.ctheta +
                                hy * self.cphi * self.ctheta)
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def iec_interaction(self, J, layer: "LayerSB"):
        """
        Energy is symmetric and computed only once
        """
        if layer is None:
            return 0
        [m1x, m1y, m1z] = self.m.get_cartesian()
        [m2x, m2y, m2z] = layer.m.get_cartesian()
        constJ = J / self.thickness
        return -constJ * (m1x * m2x + m1y * m2y + m1z * m2z)

    def grad_iec_interaction(self, J_top, J_bottom, top_layer: "LayerSB",
                             bottom_layer: "LayerSB"):
        """
        IEC gradient is not symmetric and is computed per layer
        Compute joint interaction on the layers from the top and bottom
        """
        grad1 = self.__any_grad_iec_interaction(J_top, top_layer)
        grad2 = self.__any_grad_iec_interaction(J_bottom, bottom_layer)
        return [
            grad1[0] + grad2[0], grad1[1] + grad2[1], grad1[2] + grad2[2],
            grad1[3] + grad2[3], grad1[4] + grad2[4]
        ]

    def __any_grad_iec_interaction(self, J, layer: "LayerSB"):
        # from any of the layers
        if (layer is None) or J == 0:
            return [0, 0, 0, 0, 0]
        constJ = J / self.thickness

        dEdtheta = -constJ * (self.sphi * math.sin(layer.m.phi) *
                              math.sin(layer.m.theta) * self.stheta -
                              self.stheta * math.cos(layer.m.theta) +
                              math.sin(layer.m.theta) * self.cphi *
                              math.cos(layer.m.phi) * self.ctheta)

        dEdphi = -constJ * (-self.sphi * self.stheta * math.sin(
            layer.m.theta) * math.cos(layer.m.phi) + math.sin(layer.m.phi) *
                            self.stheta * math.sin(layer.m.theta) * self.cphi)

        d2Edtheta2 = -constJ * (
            -self.sphi * layer.sphi * self.stheta * layer.stheta -
            self.stheta * layer.stheta * self.cphi * layer.cphi -
            self.ctheta * layer.ctheta)

        d2Edphi2 = -constJ * (
            -self.sphi * layer.sphi * self.stheta * layer.stheta -
            self.stheta * layer.stheta**self.cphi * layer.cphi)

        d2Edphidtheta = -constJ * (
            -self.sphi * layer.stheta * layer.cphi * self.ctheta +
            layer.sphi * layer.stheta * self.cphi * self.ctheta)

        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def surface_anisotropy(self):
        return (-self.Ks + 0.5 * self.Ms * mu0) * self.ctheta**2

    def grad_surface_anisotropy(self):
        msq = math.pow(self.Ms, 2) * mu0
        dEdphi = 0
        d2Edphi2 = 0
        d2Edphidtheta = 0
        dEdtheta = (2 * self.Ks - msq) * self.stheta * self.ctheta
        d2Edtheta2 = (2 * self.Ks - msq) * (self.ctheta**2 - self.stheta**2)
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def volume_anisotropy(self):
        """
        E_k = K (m o a)^2
        """
        ax = math.cos(self.Kv.phi)
        ay = math.sin(self.Kv.phi)
        mx, my, _ = self.m.get_cartesian()
        return -self.Kv.mag * (mx * ax + my * ay)**2

    def grad_volume_anisotropy(self):
        """
        E_k/dtheta = ~ (a_x * m_x) * m_x *dm_x/dtheta
        """
        Kv = self.Kv.mag
        angdiff = self.Kv.phi - self.m.phi
        dEdtheta = -2 * Kv * self.stheta * self.ctheta * (math.cos(angdiff)**2)
        dEdphi = -2 * Kv * (self.stheta**
                            2) * math.sin(angdiff) * math.cos(angdiff)
        d2Edtheta2 = -2 * Kv * math.cos(2 * self.m.theta) * (math.cos(angdiff)
                                                             **2)
        d2Edphi2 = 2 * Kv * (self.stheta**2) * math.cos(2 * angdiff)
        d2Edphidtheta = -0.5 * Kv * (
            math.cos(-2 * self.Kv.phi + 2 * self.m.phi + 2 * self.m.theta) -
            math.cos(2 * self.Kv.phi - 2 * self.m.phi + 2 * self.m.theta))
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def compute_grad_energy(self, Hinplane: VectorObj, Jtop: float,
                            Jbottom: float, top_layer: "LayerSB",
                            bottom_layer: "LayerSB"):
        g1 = self.grad_ext_field(Hinplane)
        g2 = self.grad_surface_anisotropy()
        g3 = self.grad_volume_anisotropy()
        g4 = self.grad_iec_interaction(Jtop, Jbottom, top_layer, bottom_layer)
        return [g1[i] + g2[i] + g3[i] + g4[i] for i in range(5)]

    def compute_energy(self, Hinplane: VectorObj, Jtop: float, Jbottom: float,
                       top_layer: "LayerSB", bottom_layer: "LayerSB"):
        e1 = self.ext_field(Hinplane)
        e2 = self.surface_anisotropy()
        e3 = self.volume_anisotropy()
        e4 = self.iec_interaction(
            Jbottom, bottom_layer) + self.iec_interaction(Jtop, top_layer)
        evec = [e1, e2, e3, e4]
        return evec

    def get_current_position(self):
        return [self.m.theta, self.m.phi]

    def set_current_position(self, pos):
        self.m.theta = pos[0]
        self.m.phi = pos[1]

    def compute_frequency_at_equilibrum(self, Hinplane: VectorObj, Jtop: float,
                                        Jbottom: float, top_layer: "LayerSB",
                                        bottom_layer: "LayerSB"):
        """Computes the resonance frequency (FMR) of the layers
        :param Hinplance: vector that describes the applied H.
        :param Jtop: IEC constant from the layer above the current one.
        :param Jbottom: IEC constant from the layer below the current one.
        :param top_layer: LayerSB definition of the layer above the current one.
        :param bottom layer: LayerSB definition of the layer below the current one."""
        (_, _, d2Edtheta2, d2Edphi2,
         d2Edphidtheta) = self.compute_grad_energy(Hinplane, Jbottom, Jtop,
                                                   top_layer, bottom_layer)

        if self.stheta != 0.:
            fmr = (d2Edtheta2 * d2Edphi2 - d2Edphidtheta**2) / math.pow(
                self.stheta * self.Ms, 2)
            if fmr > 0:
                fmr = math.sqrt(fmr) * gamma_rad / (2 * math.pi)
        else:
            return -4
        return fmr


class SmitBeljersModel:
    """
    Smits-Beljers model for the energy of a multilayered system
    It allows to compute the FMR of the system.
    """

    def __init__(self,
                 layers: List[LayerSB],
                 Hext: VectorObj,
                 J: List[float] = [],
                 silent: bool = False) -> None:
        """
        :param layers: list of LayerSB, layer definitions.
        :param Hext: applied external magnetic field vector. Defined for all the layers.
        :param J: the list of IEC constants. 0 element defines coupling between 0 and 1 layer, etc.
        :param silent: debug mode? defaults to False.
        """
        if len(layers) != (len(J) + 1):
            raise ValueError("Number of layers must be equal to len(J) + 1")
        self.J = J
        self.layers: List[LayerSB] = layers
        self.energy = np.zeros((len(self.layers), ), dtype=np.float32)
        self.current_position = np.zeros((len(self.layers), 2),
                                         dtype=np.float32)
        self.grad_energy = np.zeros((len(self.layers), 5), dtype=np.float32)
        self.Hext = Hext
        self.silent = silent
        self.ener_labels = ['external', 'surface', 'volume', 'iec']
        self.history = defaultdict(list)

    def compute_energy_step(self):
        """
        Note that symmetric coupling is only computed once and not per layer
        """
        for i, layer in enumerate(self.layers):
            if i > 0:
                bottom_layer_handle = self.layers[i - 1]
                Jbottom = self.J[i - 1]
            else:
                bottom_layer_handle = None
                Jbottom = 0
            if i < len(self.layers) - 1:
                top_layer_handle = self.layers[i + 1]
                Jtop = self.J[i]
            else:
                top_layer_handle = None
                Jtop = 0

            evec = layer.compute_energy(self.Hext, Jtop, Jbottom,
                                        top_layer_handle, bottom_layer_handle)
            self.energy[i] = sum(evec)
            self.current_position[i, :] = layer.get_current_position()
            self.grad_energy[i, :] = layer.compute_grad_energy(
                self.Hext, Jtop, Jbottom, top_layer_handle,
                bottom_layer_handle)
            for ener_val, ener_label in zip(evec, self.ener_labels):
                self.history[f"energy_{ener_label}_{i}"].append(ener_val)
            self.history[f"energy_{i}"].append(self.energy[i])
            self.history[f"theta_{i}"].append(layer.m.theta)
            self.history[f"phi_{i}"].append(layer.m.phi)

    def update_layers(self, position_vector):
        # update layers
        for l, layer in enumerate(self.layers):
            layer.set_current_position(position_vector[l])

    def position_difference(self, position_vector):
        # compute difference between current position and new position
        return np.linalg.norm(self.current_position - position_vector)

    def gradient_descent(self,
                         max_steps: int,
                         tol: float = 1e-8,
                         learning_rate=1e-3):
        """
        Main gradient descent algorithm. Currently implements a basic version.
        TODO: implement a more advanced version -- conjugate gradient descent
        or Nesterov's accelerated gradient descent
        """
        warnings.warn("Use ADAM gradient descent for better results!")
        for i in tqdm(range(int(max_steps)), mininterval=0.5):
            self.compute_energy_step()
            # compute gradient update
            new_position = self.current_position - learning_rate * self.grad_energy[:, :
                                                                                    2]
            if self.position_difference(new_position) < tol:
                break
            self.update_layers(new_position)
        self.print_summary(i)

    def get_fmr(self, layer_index: int) -> float:
        """Computes FMR values for each of the layers in the model."""
        if layer_index > 0:
            bottom_layer_handle = self.layers[layer_index - 1]
            Jbottom = self.J[layer_index - 1]
        else:
            bottom_layer_handle = None
            Jbottom = 0
        if layer_index < len(self.layers) - 1:
            top_layer_handle = self.layers[layer_index + 1]
            Jtop = self.J[layer_index]
        else:
            top_layer_handle = None
            Jtop = 0
        return self.layers[layer_index].compute_frequency_at_equilibrum(
            self.Hext, Jtop, Jbottom, top_layer_handle, bottom_layer_handle)

    def print_summary(self, steps):
        """Prints summary of the solution if silent param is true."""
        print(f"Gradient descent finished in {steps} steps")
        print(f"Final energy: {self.energy.sum()*1e3:.4f} mJ/m^2")
        print(f"Final position: {np.rad2deg(self.current_position)}")
        for i, _ in enumerate(self.layers):
            fmr = self.get_fmr(i)
            print(f"FMR: {fmr/1e9:.4f} GHz")

    def adam_gradient_descent(self,
                              max_steps: int,
                              tol: float = 1e-8,
                              learning_rate: float = 1e-4,
                              first_momentum_decay: float = 0.9,
                              second_momentum_decay: float = 0.999):
        """
        A naive implementation of Adam gradient descent.
        See: ADAM: A METHOD FOR STOCHASTIC OPTIMIZATION, Kingma et Ba, 2015
        :param max_steps: maximum number of gradient steps.
        :param tol: tolerance of the solution.
        :param learning_rate: the learning rate (descent speed).
        :param first_momentum_decay: constant for the first momentum.
        :param second_momentum_decay: constant for the second momentum.
        """
        step = 0
        m = np.zeros_like(self.grad_energy[:, :2])  # first momentum
        v = np.zeros_like(self.grad_energy[:, :2])  # second momentum
        eps = 1e-12
        while True:
            step += 1
            self.compute_energy_step()
            m = first_momentum_decay * m + (
                1. - first_momentum_decay) * self.grad_energy[:, :2]
            v = second_momentum_decay * v + (
                1. - second_momentum_decay) * self.grad_energy[:, :2]**2
            m_hat = m / (1. - first_momentum_decay**step)
            v_hat = v / (1. - second_momentum_decay**step)
            new_position = self.current_position - learning_rate * m_hat / (
                np.sqrt(v_hat) + eps)
            if step > max_steps:
                break
            if self.position_difference(new_position) < tol:
                break
            self.update_layers(new_position)
        if not self.silent:
            self.print_summary(step)
