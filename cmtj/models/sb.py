import math
from dataclasses import dataclass
from functools import lru_cache
from typing import List

import numpy as np
from tqdm import tqdm

from ..utils import TtoAm, gyromagnetic_ratio, mu0


@dataclass
class VectorObj:
    theta: float  # in radians
    phi: float  # rad
    mag: float = 1

    @lru_cache
    def get_spherical(self):
        return [
            self.mag * math.sin(self.theta) * math.cos(self.phi),
            self.mag * math.sin(self.theta) * math.sin(self.phi),
            self.mag * math.cos(self.theta),
        ]


@dataclass
class LayerSB:
    thickness: float
    m: VectorObj
    Kv: VectorObj
    Ks: float
    Ms: float
    energy: float = 0

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
        hx, hy, hz = Hinplane.get_spherical()
        return -self.Ms * (hx * self.stheta * self.cphi +
                           hy * self.sphi * self.ctheta + hz * self.ctheta)

    def grad_ext_field(self, Hinplane: VectorObj):
        hx, hy, hz = Hinplane.get_spherical()
        dEdtheta = -self.Ms * (hx * self.cphi * self.ctheta +
                               hy * self.sphi * self.ctheta - hz * self.stheta)

        dEdphi = -self.Ms * (-hx * self.sphi * self.stheta +
                             hy * self.stheta * self.cphi)

        d2Edtheta2 = -self.Ms * (hx * self.stheta * self.cphi - hy *
                                 self.sphi * self.stheta - hz * self.ctheta)
        d2Edphi2 = -self.Ms * (-hx * self.stheta * self.cphi -
                               hy * self.sphi * self.stheta)
        d2Edphidtheta = self.Ms * (-hx * self.sphi * self.ctheta +
                                   hy * self.cphi * self.ctheta)
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def iec_interaction(self, J, layer: "LayerSB"):
        """
        Energy is symmetric and computed only once
        """
        if layer is None:
            return 0
        [m1x, m1y, m1z] = self.m.get_spherical()
        [m2x, m2y, m2z] = layer.m.get_spherical()
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
            -self.sphi(layer.stheta * layer.cphi * self.ctheta +
                       layer.sphi * layer.stheta * self.cphi * self.ctheta))

        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def surface_anisotropy(self):
        return (-self.Ks + 0.5 * (self.Ms**2) / mu0) * self.ctheta

    def grad_surface_anisotropy(self):
        pref = (-self.Ks + 0.5 * (self.Ms**2) / mu0)
        dEdtheta = -pref * self.stheta
        dEdphi = 0
        d2Edtheta2 = -pref * self.ctheta
        d2Edphi2 = 0
        d2Edphidtheta = 0
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def volume_anisotropy(self):
        ax = math.cos(self.Kv.phi)
        ay = math.sin(self.Kv.phi)
        mx, my, _ = self.m.get_cartesian()
        return -self.Kv.mag * (mx * ax + my * ay)

    def grad_volume_anisotropy(self):
        ax = math.cos(self.Kv.phi)
        ay = math.sin(self.Kv.phi)
        Kv = self.Kv.mag
        dEdtheta = -Kv * (self.cphi * self.ctheta * ax +
                          self.sphi * self.ctheta * ay)
        dEdphi = -Kv * (-self.sphi * self.stheta * ax +
                        self.cphi * self.stheta * ay)
        d2Edtheta2 = -Kv * (-self.cphi * self.stheta * ax -
                            self.sphi * self.stheta * ay)

        d2Edphi2 = -Kv * (-self.cphi * self.stheta * ax -
                          self.sphi * self.stheta * ay)
        d2Edphidtheta = -Kv * (-self.sphi * self.ctheta * ax +
                               self.cphi * self.ctheta * ay)
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
        return e1 + e2 + e3 + e4

    def get_current_position(self):
        return [self.m.theta, self.m.phi]

    def set_current_position(self, pos):
        self.m.theta = pos[0]
        self.m.phi = pos[1]

    def compute_frequency_at_equilibrum(self, Hinplane: VectorObj, Jtop: float,
                                        Jbottom: float, top_layer: "LayerSB",
                                        bottom_layer: "LayerSB"):
        (_, _, d2Edtheta2, d2Edphi2,
         d2Edphidtheta) = self.compute_grad_energy(Hinplane, Jbottom, Jtop,
                                                   top_layer, bottom_layer)
        t = self.m.theta
        if t != 0.:
            fmr = d2Edtheta2 * d2Edphi2 - math.sqrt(d2Edphidtheta) / t
        else:
            fmr = -4
        frequency = gyromagnetic_ratio * math.sqrt(fmr) / (
            TtoAm * self.Ms * 2 * math.pi * self.thickness)
        return frequency


class SmitBeljersModel:
    """
    Smits-Beljers model for the energy of a multilayered system
    It allows to compute the FMR of the system.
    """

    def __init__(self,
                 layers: List[LayerSB],
                 Hext: VectorObj,
                 J: List[float] = [],
                 eta: float = 1e-2) -> None:

        if len(layers) != (len(J) + 1):
            raise ValueError(
                "Number of layers must be equal to number of J + 1")
        self.J = J
        self.layers: List[LayerSB] = layers
        self.energy = np.zeros((len(self.layers), 1))
        self.current_position = np.zeros((len(self.layers), 2))
        self.grad_energy = np.zeros((len(self.layers), 5))
        self.eta = eta
        self.Hext = Hext

    def compute_coupling_energy(self, J, layerA: "LayerSB", layerB: "LayerSB"):
        """
        Compute the coupling energy between two layers.
        This should be symmetric and computed only once.
        Gradients are computed per layer, they are not symmetric.
        For now IEC is the only symmetric contribution.
        """
        iec_interaction = -J * (layerA.stheta * layerB.stheta *
                                math.cos(layerA.m.phi - layerB.m.phi) +
                                layerA.ctheta * layerB.ctheta)
        return iec_interaction

    def compute_energy_step(self):
        """
        Note that symmetric coupling is only computed once and not per layer
        """
        for i, layer in enumerate(self.layers):
            if i > 0:
                self.compute_coupling_energy(self.layers[i - 1],
                                             self.layers[i])
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

            self.energy[i] = layer.compute_energy(self.Hext, Jtop, Jbottom,
                                                  top_layer_handle,
                                                  bottom_layer_handle)
            self.current_position[i, :] = layer.get_current_position()
            self.grad_energy[i, :] = layer.compute_grad_energy(
                self.Hext, Jtop, Jbottom, top_layer_handle,
                bottom_layer_handle)

    def update_layers(self, position_vector):
        # update layers
        for l, layer in enumerate(self.layers):
            layer.set_current_position(position_vector[l])

    def position_difference(self, position_vector):
        # compute difference between current position and new position
        return np.linalg.norm(self.current_position - position_vector)

    def gradient_descent(self, max_steps: int, tol: float = 1e-8):
        """
        Main gradient descent algorithm. Currently implements a basic version.
        TODO: implement a more advanced version -- conjugate gradient descent
        or Nesterov's accelerated gradient descent
        """
        for i in tqdm(range(int(max_steps)), mininterval=0.5):
            self.compute_energy_step()
            # compute gradient update
            new_position = self.current_position - self.eta * self.grad_energy[:, :
                                                                               2]
            if self.position_difference(new_position) < tol:
                break
            self.update_layers(new_position)
        self.print_summary(i)

    def print_summary(self, steps):
        print(f"Gradient descent finished in {steps} steps")
        print(f"Final energy: {self.energy.sum()*1e6:.4f} uJ")
        print(f"Final position: {np.rad2deg(self.current_position)}")
        # for layer in self.layers:
        #     print(
        #         f"FMR: {layer.compute_frequency_at_equilibrum(self.Hext, 0, 0, None, None):.2f} GHz"
        #     )

    def adam_gradient_descent(self,
                              max_steps: int,
                              tol: float = 1e-8,
                              learning_rate: float = 1e-4,
                              first_momentum_decay: float = 0.9,
                              second_momentum_decay: float = 0.999):
        """
        A naive implementation of Adam gradient descent.
        """
        step = 0
        m = 0  # first momentum
        v = 0  # second momentum
        eps = 1e-8
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
        self.print_summary(step)
