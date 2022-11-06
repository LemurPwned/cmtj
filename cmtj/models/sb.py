import math
from dataclasses import dataclass
from typing import List

import numpy as np
from tqdm import tqdm

from ..utils import TtoAm
"""
3 gradients:
dE/dtheta
dE/dphi
dE/dt => dE^2/dt^2
dE^2/dtdf
dE^2/df^2
"""


@dataclass
class VectorObj:
    theta: float  # in radians
    phi: float  # rad
    mag: float = 1


@dataclass
class LayerSB:
    thickness: float
    m: VectorObj
    K: VectorObj
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

    def inplane_anisotropy(self):
        return -self.K.mag * self.thickness * math.pow(
            math.cos(self.m.theta - self.K.phi), 2) * math.pow(
                math.sin(self.m.theta - self.K.phi), 2)

    def grad_inplane_anisotropy(self):
        dEdtheta = self.thickness * 2 * self.stheta * self.ctheta * (
            -self.K.mag * math.pow(math.cos(self.m.phi - self.K.phi), 2) -
            TtoAm * math.pow(self.Ms, 2) / 2)
        dEdphi = self.thickness * math.pow(
            self.stheta, 2) * self.K.mag * math.sin(2 *
                                                    (self.m.phi - self.K.phi))

        d2Edtheta2 = (-math.pow(math.cos(self.m.phi - self.K.phi), 2)) - (
            (self.Ms**2) * TtoAm / 2)(
                self.thickness * 2 * math.cos(2 * self.m.theta))
        d2Edphi2 = (2 * self.K.mag * math.pow(self.stheta)) * math.cos(
            2 * (self.m.phi - self.K.phi)) * self.thickness
        d2Edphittheta = (self.K.mag * math.sin(2 * self.m.theta) *
                         math.sin(2 *
                                  (self.m.phi - self.K.phi))) * self.thickness
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphittheta]

    def perpendicular_anisotropy(self):
        return -self.K.mag * self.thickness * math.pow(self.ctheta, 2)

    def grad_perpendicular_anisotropy(self):
        dEdtheta = self.thickness * 2 * self.stheta * self.ctheta * (
            self.K.mag - TtoAm * math.pow(self.Ms, 2) / 2)
        d2Edtheta2 = (self.K.mag - (self.Ms**2) * TtoAm /
                      2) * self.thickness * 2 * math.cos(2 * self.m.theta)
        return [dEdtheta, 0, d2Edtheta2, 0, 0]

    def ext_field(self, Hinplane: VectorObj):
        return -self.Ms * self.thickness * Hinplane.mag * (
            self.stheta * math.sin(Hinplane.phi) * self.cphi +
            self.ctheta * math.cos(Hinplane.phi))

    def grad_ext_field(self, Hinplane: VectorObj):
        pre = -Hinplane.mag * self.Ms * self.thickness
        sH = math.sin(Hinplane.phi)
        cH = math.cos(Hinplane.phi)
        dEdtheta = pre * (self.ctheta * sH * self.cphi - self.stheta * cH)
        dEdphi = pre * (self.stheta * sH * self.sphi)

        # second energy derivative wrt theta
        d2Edtheta2 = pre * (self.stheta * math.phi(Hinplane.phi) * self.cphi +
                            self.ctheta * math.cos(Hinplane.phi))
        # derivative wrt phi
        d2Edphi2 = pre * (self.stheta(math.sin(Hinplane.phi) * self.cphi))
        # mixed derivative
        d2Edphittheta = pre * (self.ctheta * math.sin(Hinplane.phi) *
                               self.cphi)

        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphittheta]

    def demagnetisation(self):
        return -TtoAm / 2 * self.Ms * self.Ms * self.thickness * (self.ctheta**
                                                                  2)

    def grad_demagnetisation(self):
        pre = TtoAm * self.Ms * self.Ms * self.thickness
        dEdtheta = pre * (2 * self.stheta * self.ctheta)
        d2Edtheta2 = pre * (-2 * math.cos(2 * self.m.theta) * self.stheta +
                            2 * math.sin(2 * self.m.theta) * self.ctheta)
        return [dEdtheta, 0, d2Edtheta2, 0, 0]

    def compute_energy(self, Hinplane: VectorObj):
        energy = self.ext_field(
            Hinplane=Hinplane) + self.perpendicular_anisotropy(
            ) + self.inplane_anisotropy()
        return energy

    def iec_interaction(self, J, layer: "LayerSB"):
        """
        Energy is symmetric and computed only once
        """
        return -J * (
            self.stheta * layer.stheta * math.cos(self.m.phi - layer.m.phi) +
            self.ctheta * layer.ctheta)

    def grad_iec_interaction(self, J, top_layer: "LayerSB",
                             bottom_layer: "LayerSB"):
        """
        IEC gradient is not symmetric and is computed per layer
        Compute joint interaction on the layers from the top and bottom
        """
        grad1 = self.__any_grad_iec_interaction(J, top_layer)
        grad2 = self.__any_grad_iec_interaction(J, bottom_layer)
        return [
            grad1[0] + grad2[0], grad1[1] + grad2[1], grad1[2] + grad2[2],
            grad1[3] + grad2[3], grad1[4] + grad2[4]
        ]

    def __any_grad_iec_interaction(self, J, any_layer: "LayerSB"):
        # from any of the layers
        if any_layer is None:
            return [0, 0, 0, 0, 0]
        dEdtheta = -J * (
            self.ctheta * self.stheta * math.cos(self.m.phi - any_layer.m.phi)
            - self.stheta * math.cos(self.m.theta - any_layer.m.theta))
        dEdphi = -J * math.sin(
            any_layer.m.theta) * self.stheta * math.sin(self.m.phi -
                                                        any_layer.m.phi)

        d2Edtheta2 = -J * (-self.stheta * any_layer.stheta *
                           math.cos(self.m.phi - any_layer.m.phi) -
                           self.ctheta * any_layer.ctheta)
        d2Edphi2 = -J * (-self.stheta * any_layer.stheta *
                         math.cos(self.m.phi - any_layer.m.phi))
        d2Edfdt = -J * (-self.ctheta * any_layer.stheta *
                        math.sin(self.m.phi - any_layer.m.phi))
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edfdt]

    def compute_grad_energy(self, Hinplance: VectorObj, Jbottom: float,
                            Jtop: float, top_layer: "LayerSB",
                            bottom_layer: "LayerSB"):
        g1 = self.grad_ext_field(Hinplance)
        g2 = self.grad_perpendicular_anisotropy()
        g3 = self.grad_inplane_anisotropy()
        g4 = self.grad_demagnetisation()
        g5 = self.grad_iec_interaction(Jbottom, Jtop, top_layer, bottom_layer)
        return [g1[i] + g2[i] + g3[i] + g4[i] + g5[i] for i in range(5)]

    def get_current_position(self):
        return [self.m.theta, self.m.phi]

    def set_current_position(self, pos):
        self.m.theta = pos[0]
        self.m.phi = pos[1]


class SmitBeljersModel:
    """
    Smits-Beljers model for the energy of a multilayered system
    It allows to compute the FMR of the system.
    """

    def __init__(self,
                 layers: List[LayerSB],
                 J: List[float],
                 Hext: VectorObj,
                 eta: float = 1e-2) -> None:

        if len(layers) != (len(J) + 1):
            raise ValueError(
                "Number of layers must be equal to number of J + 1")
        self.J = J
        self.layers: List[LayerSB] = layers
        self.energy = np.zeros((len(self.layers), 1))
        self.current_position = np.zeros((len(self.layers), 2))
        self.grad_energy = np.zeros((len(self.layers), 2))
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

            self.energy[i] = layer.compute_energy(self.Hext, Jbottom, Jtop,
                                                  top_layer_handle,
                                                  bottom_layer_handle)
            self.current_position[i, :] = layer.get_current_position()
            self.grad_energy[i, :] = layer.compute_grad_energy(self.Hext)

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
            new_position = self.current_position - self.eta * self.grad_energy
            self.update_layers(new_position)
            pos_difference = np.abs(self.current_position - new_position).sum()
            if ((pos_difference < tol)
                    or (np.linalg.norm(self.grad_energy) < tol)):
                break
            self.current_position = new_position

        self.print_summary(i)

    def print_summary(self, steps):
        print(f"Gradient descent finished in {steps} steps")
        print(f"Final energy: {self.energy.sum()*1e6:.4f} uJ")
        print(f"Final position: {np.rad2deg(self.current_position)}")

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
                1 - first_momentum_decay) * self.grad_energy
            v = second_momentum_decay * v + (
                1 - second_momentum_decay) * self.grad_energy**2
            m_hat = m / (1 - first_momentum_decay**step)
            v_hat = v / (1 - second_momentum_decay**step)
            new_position = self.current_position - learning_rate * m_hat / (
                np.sqrt(v_hat) + eps)
            if step > max_steps:
                break
            if self.position_difference(new_position) < tol:
                break
            self.update_layers(new_position)
        self.print_summary(step)
