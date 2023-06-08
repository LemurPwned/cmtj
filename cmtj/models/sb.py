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

    def grad_ext_field(self, Hinplane: VectorObj, full_grad: bool):
        hx, hy, hz = Hinplane.get_cartesian()
        pref = -self.Ms * mu0
        dEdtheta = pref * (hx * self.cphi * self.ctheta +
                           hy * self.sphi * self.ctheta - hz * self.stheta)

        dEdphi = pref * (-hx * self.sphi * self.stheta +
                         hy * self.stheta * self.cphi)
        if not full_grad:
            return [dEdtheta, dEdphi, 0, 0, 0]
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
                             bottom_layer: "LayerSB", full_grad: bool):
        """
        IEC gradient is not symmetric and is computed per layer
        Compute joint interaction on the layers from the top and bottom
        """
        grad1 = self.__any_grad_iec_interaction(J_top,
                                                top_layer,
                                                full_grad=full_grad)
        grad2 = self.__any_grad_iec_interaction(J_bottom,
                                                bottom_layer,
                                                full_grad=full_grad)
        return [
            grad1[0] + grad2[0], grad1[1] + grad2[1], grad1[2] + grad2[2],
            grad1[3] + grad2[3], grad1[4] + grad2[4]
        ]

    def __any_grad_iec_interaction(self, J, layer: "LayerSB", full_grad: bool):
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
        if not full_grad:
            return [dEdtheta, dEdphi, 0, 0, 0]
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

    def dmi_interaction(self, D, layer: "LayerSB"):
        """
        DMi energy is Edmi = D z(m1 x m2)
        """
        if layer is None:
            return 0
        [m1x, m1y, _] = self.m.get_cartesian()
        [m2x, m2y, _] = layer.m.get_cartesian()
        return D * (m1x * m2y - m1y * m2x)

    def grad_dmi_interaction(self, D_top, D_bottom, top_layer: "LayerSB",
                             bottom_layer: "LayerSB", full_grad: bool):
        """
        DMI energy and its gradient are both not symmetric and is computed per layer
        Compute joint interaction on the layers from the top and bottom
        """
        grad1 = self.__any_grad_dmi_interaction(D_top,
                                                top_layer,
                                                full_grad=full_grad)
        grad2 = self.__any_grad_dmi_interaction(D_bottom,
                                                bottom_layer,
                                                full_grad=full_grad)
        return [
            grad1[0] + grad2[0], grad1[1] + grad2[1], grad1[2] + grad2[2],
            grad1[3] + grad2[3], grad1[4] + grad2[4]
        ]

    def __any_grad_dmi_interaction(self, D, layer: "LayerSB", full_grad: bool):
        if (layer is None) or D == 0:
            return [0, 0, 0, 0, 0]
        constD = D

        dEdtheta = -constD * math.sin(
            layer.m.theta) * math.sin(self.m.phi - layer.m.phi) * self.ctheta

        dEdphi = -constD * self.stheta * math.sin(
            layer.m.theta) * math.cos(self.m.phi - layer.m.phi)
        if not full_grad:
            return [dEdtheta, dEdphi, 0, 0, 0]
        d2Edtheta2 = constD * self.stheta * math.sin(self.m.phi - layer.m.phi)

        d2Edphi2 = d2Edtheta2

        d2Edphidtheta = -D * math.sin(
            layer.m.phi) * self.ctheta * math.cos(self.m.phi - layer.m.phi)

        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def surface_anisotropy(self):
        return (-self.Ks + 0.5 * self.Ms * mu0) * self.ctheta**2

    def grad_surface_anisotropy(self, full_grad: bool):
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

    def grad_volume_anisotropy(self, full_grad: bool):
        """
        E_k/dtheta = ~ (a_x * m_x) * m_x *dm_x/dtheta
        """
        Kv = self.Kv.mag
        angdiff = self.Kv.phi - self.m.phi
        dEdtheta = -2 * Kv * self.stheta * self.ctheta * (math.cos(angdiff)**2)
        dEdphi = -2 * Kv * (self.stheta**
                            2) * math.sin(angdiff) * math.cos(angdiff)
        if not full_grad:
            return [dEdtheta, dEdphi, 0, 0, 0]
        d2Edtheta2 = -2 * Kv * math.cos(2 * self.m.theta) * (math.cos(angdiff)
                                                             **2)
        d2Edphi2 = 2 * Kv * (self.stheta**2) * math.cos(2 * angdiff)
        d2Edphidtheta = -0.5 * Kv * (
            math.cos(-2 * self.Kv.phi + 2 * self.m.phi + 2 * self.m.theta) -
            math.cos(2 * self.Kv.phi - 2 * self.m.phi + 2 * self.m.theta))
        return [dEdtheta, dEdphi, d2Edtheta2, d2Edphi2, d2Edphidtheta]

    def compute_grad_energy(self,
                            Hinplane: VectorObj,
                            Jtop: float,
                            Jbottom: float,
                            Dtop: float,
                            Dbottom: float,
                            top_layer: "LayerSB",
                            bottom_layer: "LayerSB",
                            full_grad: bool = False):
        g1 = self.grad_ext_field(Hinplane, full_grad)
        g2 = self.grad_surface_anisotropy(full_grad)
        g3 = self.grad_volume_anisotropy(full_grad)
        g4 = self.grad_iec_interaction(Jtop, Jbottom, top_layer, bottom_layer,
                                       full_grad)
        g5 = self.grad_dmi_interaction(Dtop, Dbottom, top_layer, bottom_layer,
                                       full_grad)
        return [g1[i] + g2[i] + g3[i] + g4[i] + g5[i] for i in range(5)]

    def compute_energy(self, Hinplane: VectorObj, Jtop: float, Jbottom: float,
                       Dtop: float, Dbottom: float, top_layer: "LayerSB",
                       bottom_layer: "LayerSB"):
        e1 = self.ext_field(Hinplane)
        e2 = self.surface_anisotropy()
        e3 = self.volume_anisotropy()
        e4 = self.iec_interaction(
            Jbottom, bottom_layer) + self.iec_interaction(Jtop, top_layer)
        e5 = self.dmi_interaction(Dtop, top_layer) + self.dmi_interaction(
            Dbottom, bottom_layer)
        evec = [e1, e2, e3, e4, e5]
        return evec

    def get_current_position(self):
        return [self.m.theta, self.m.phi]

    def set_current_position(self, pos):
        self.m.theta = pos[0]
        self.m.phi = pos[1]

    def compute_frequency_at_equilibrum(self, Hinplane: VectorObj, Jtop: float,
                                        Jbottom: float, Dtop: float,
                                        Dbottom: float, top_layer: "LayerSB",
                                        bottom_layer: "LayerSB"):
        """Computes the resonance frequency (FMR) of the layers
        :param Hinplance: vector that describes the applied H.
        :param Jtop: IEC constant from the layer above the current one.
        :param Jbottom: IEC constant from the layer below the current one.
        :param top_layer: LayerSB definition of the layer above the current one.
        :param bottom layer: LayerSB definition of the layer below the current one."""
        (_, _, d2Edtheta2, d2Edphi2,
         d2Edphidtheta) = self.compute_grad_energy(Hinplane,
                                                   Jtop,
                                                   Jbottom,
                                                   Dtop,
                                                   Dbottom,
                                                   top_layer,
                                                   bottom_layer,
                                                   full_grad=True)

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
                 D: List[float] = [],
                 silent: bool = False) -> None:
        """
        :param layers: list of LayerSB, layer definitions.
        :param Hext: applied external magnetic field vector. Defined for all the layers.
        :param J: list of IEC constants. 0 element defines coupling between 0 and 1 layer, etc.
        :param D: list of DMI constants.
        :param silent: debug mode? defaults to False.
        """
        if len(layers) != (len(J) + 1):
            raise ValueError("Number of layers must be equal to len(J) + 1")
        self.J = J
        self.D = D
        self.layers: List[LayerSB] = layers
        self.energy = np.zeros((len(self.layers), ), dtype=np.float32)
        self.current_position = np.zeros((len(self.layers), 2),
                                         dtype=np.float32)
        self.grad_energy = np.zeros((len(self.layers), 5), dtype=np.float32)
        self.Hext = Hext
        self.silent = silent
        self.ener_labels = ['external', 'surface', 'volume', 'iec']
        self.history = defaultdict(list)

    def __get_interaction_constant(self, layers: List[LayerSB],
                                   layer_indx: int,
                                   interaction_constant_list: List[float]):
        """Simple function to figure out top and bottom handles and interaction constants"""
        if layer_indx > 0:
            bottom_handle = layers[layer_indx - layer_indx]
            interaction_const_bottom = interaction_constant_list[layer_indx -
                                                                 1]
        else:
            bottom_handle = None
            interaction_const_bottom = 0
        if layer_indx < len(layers) - 1:
            top_handle = layers[layer_indx + 1]
            interaction_const_top = interaction_constant_list[layer_indx]
        else:
            top_handle = None
            interaction_const_top = 0
        return interaction_const_top, interaction_const_bottom, top_handle, bottom_handle

    def compute_energy_step(self):
        """
        Note that symmetric coupling is only computed once and not per layer
        """
        for i, layer in enumerate(self.layers):
            Jtop, Jbottom, top_layer_handle, bottom_layer_handle = self.__get_interaction_constant(
                self.layers, i, self.J)
            Dtop, Dbottom, _, _ = self.__get_interaction_constant(
                self.layers, i, self.D)
            evec = layer.compute_energy(self.Hext, Jtop, Jbottom, Dtop,
                                        Dbottom, top_layer_handle,
                                        bottom_layer_handle)
            self.energy[i] = sum(evec)
            self.current_position[i, :] = layer.get_current_position()
            self.grad_energy[i, :] = layer.compute_grad_energy(
                self.Hext, Jtop, Jbottom, Dtop, Dbottom, top_layer_handle,
                bottom_layer_handle)
            for ener_val, ener_label in zip(evec, self.ener_labels):
                self.history[f"energy_{ener_label}_{i}"].append(ener_val)
            self.history[f"energy_{i}"].append(self.energy[i])
            self.history[f"theta_{i}"].append(layer.m.theta)
            self.history[f"phi_{i}"].append(layer.m.phi)

    def transform_flat_position_to_layer_position(self, position_vector):
        position_ = [
            position_vector[n:n + 2] for n in range(0, len(position_vector), 2)
        ]
        return position_

    def compute_energy(self, position_vector):
        position_ = self.transform_flat_position_to_layer_position(
            position_vector=position_vector)
        self.update_layers(position_)
        total_energy = 0
        for i, layer in enumerate(self.layers):
            Jtop, Jbottom, top_layer_handle, bottom_layer_handle = self.__get_interaction_constant(
                self.layers, i, self.J)
            Dtop, Dbottom, _, _ = self.__get_interaction_constant(
                self.layers, i, self.D)
            evec = layer.compute_energy(self.Hext, Jtop, Jbottom, Dtop,
                                        Dbottom, top_layer_handle,
                                        bottom_layer_handle)
            total_energy += sum(evec)
        return total_energy

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
        Jtop, Jbottom, top_layer_handle, bottom_layer_handle = self.__get_interaction_constant(
            self.layers, layer_index, self.J)
        Dtop, Dbottom, _, _ = self.__get_interaction_constant(
            self.layers, layer_index, self.D)
        return self.layers[layer_index].compute_frequency_at_equilibrum(
            self.Hext, Jtop, Jbottom, Dtop, Dbottom, top_layer_handle,
            bottom_layer_handle)

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

    def nelder_mead_method(self,
                           max_steps: int,
                           tol: float = 1e-8,
                           alpha: float = 1.,
                           gamma: float = 2.,
                           rho: float = -0.5,
                           sigma: float = 0.5):
        """
        A naive implementation of the Nedler-Mead method.
        See: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
        Defaults taken from the wikipedia page.
        :param max_steps: maximum number of steps.
        :param tol: tolerance of the solution.
        """

        def simplex_centroid(points):
            return np.mean(points, axis=0)

        def reflection(point, centroid):
            return centroid + alpha * (centroid - point)

        def expansion(point, centroid):
            return centroid + gamma * (point - centroid)

        def contraction(point, centroid):
            return centroid + rho * (point - centroid)

        def shrink(points):
            return np.asarray(
                [points[0], *(points[0] + sigma * (points[1:] - points[0]))])

        def point_energies(points):
            return np.asarray([self.compute_energy(p) for p in points])

        def sort_points(points, energies):
            return points[np.argsort(energies)]

        steps = 0
        flat_coords = self.current_position.flatten()
        points = [flat_coords]
        for i in range(2):
            points.append(flat_coords +
                          np.random.normal(size=flat_coords.shape))

        points = np.asarray(points)

        while True:
            # sort points by energy
            energies = point_energies(points)
            # sort points by energy, lowest first
            points = sort_points(points, energies=energies)
            # print("Best point", points[0])
            position_ = self.transform_flat_position_to_layer_position(
                position_vector=points[0])
            self.update_layers(position_)
            if self.position_difference(points[0]) < tol:
                break
            # do not allow more than max_steps
            if steps > max_steps:
                break
            steps += 1
            # compute excluding the worst point
            centroid = simplex_centroid(points[:-1])
            # reflection
            reflected_point = reflection(points[-1], centroid)
            reflected_point_energy = self.compute_energy(reflected_point)
            # if reflected is better than the second worst point, replace worst
            if (reflected_point_energy <
                    energies[-2]) and (reflected_point_energy >= energies[0]):
                points[-1] = reflected_point
                continue

            # if reflected point is better than the best point, expand
            if reflected_point_energy < energies[0]:
                expansion_point = expansion(reflected_point, centroid)
                expansion_point_energy = self.compute_energy(expansion_point)
                if expansion_point_energy < reflected_point_energy:
                    points[-1] = expansion_point
                else:
                    points[-1] = reflected_point
                continue

            # if reflected point is better than the worst point, contract
            if reflected_point_energy < energies[-1]:
                contracted_point = contraction(reflected_point, centroid)
                contracted_point_energy = self.compute_energy(contracted_point)
                if contracted_point_energy < reflected_point_energy:
                    points[-1] = contracted_point
                else:
                    points = shrink(points)
            # if reflected point is worse than the worst point, contract
            else:
                contracted_point = contraction(points[-1], centroid)
                contracted_point_energy = self.compute_energy(contracted_point)
                if contracted_point_energy < energies[-1]:
                    points[-1] = contracted_point
                else:
                    points = shrink(points)

        # print("Best point", points[0])
        position_ = self.transform_flat_position_to_layer_position(
            position_vector=points[0])
        self.update_layers(position_)
        for i in range(len(self.layers)):
            self.history[f"energy_{i}"].append(energies[0])
            self.history[f"theta_{i}"].append(points[0][2 * i])
            self.history[f"phi_{i}"].append(points[0][2 * i + 1])

        if not self.silent:
            self.print_summary(steps)
