import math
from dataclasses import dataclass
from functools import lru_cache
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

        return [dEdtheta, dEdphi]

    def perpendicular_anisotropy(self):
        return -self.K.mag * self.thickness * math.pow(self.ctheta, 2)

    def grad_perpendicular_anisotropy(self):
        dEdtheta = self.thickness * 2 * self.stheta * self.ctheta * (
            self.K.mag - TtoAm * math.pow(self.Ms, 2) / 2)
        dEdphi = 0
        return [dEdtheta, dEdphi]

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
        return [dEdtheta, dEdphi]

    def demagnetisation(self):
        return -TtoAm / 2 * self.Ms * self.Ms * self.thickness * (self.ctheta**
                                                                  2)

    def grad_demagnetisation(self):
        dEdtheta = TtoAm * self.Ms * self.Ms * self.thickness * (
            2 * self.stheta * self.ctheta)
        dEdphi = 0
        return [dEdtheta, dEdphi]

    def compute_energy(self, Hinplane: VectorObj):
        energy = self.ext_field(
            Hinplane=Hinplane) + self.perpendicular_anisotropy(
            ) + self.inplane_anisotropy()
        return energy

    def iec_interaction(self, J, layer: "LayerSB"):
        return -J * (
            self.stheta * layer.stheta * math.cos(self.m.phi - layer.m.phi) +
            self.ctheta * layer.ctheta)

    def grad_iec_interaction(self, J, layer: "LayerSB"):
        dEdtheta = -J * (
            self.ctheta * self.stheta * math.cos(self.m.phi - layer.m.phi) -
            self.stheta * math.cos(self.m.theta - layer.m.theta))
        dEdphi = -J * math.sin(
            layer.m.theta) * self.stheta * math.sin(self.m.phi - layer.m.phi)
        return [dEdtheta, dEdphi]

    def compute_grad_energy(self, Hinplance: VectorObj, Jbottom: float,
                            Jtop: float, top_layer: "LayerSB",
                            bottom_layer: "LayerSB"):
        g1 = self.grad_ext_field(Hinplance)
        g2 = self.grad_perpendicular_anisotropy()
        g3 = self.grad_inplane_anisotropy()
        g4 = self.grad_demagnetisation()
        g5 = self.grad_iec_interaction(Jbottom, Jtop, top_layer, bottom_layer)

        return [
            g1[0] + g2[0] + g3[0] + g4[0] + g5[0],
            g1[1] + g2[1] + g3[1] + g4[1] + g5[1]
        ]

    def other_gradients(self, Htotal: VectorObj):
        # second energy derivative wrt theta
        pref = self.Ms * self.thickness * Htotal.mag
        d2edt2 = pref * (self.stheta * math.phi(Htotal.phi) * self.cphi +
                         self.ctheta * math.cos(Htotal.phi))
        # derivative wrt phi
        d2edp2 = pref * (self.stheta(math.sin(Htotal.phi) * self.cphi))
        # mixed derivative
        d2edft = pref * (self.ctheta * math.sin(Htotal.phi) * self.cphi)

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

    def compute_energy_step(self):
        """
        # TODO: how to main
        """
        for i, layer in enumerate(self.layers):
            if (i - 1) < 0:
                Jbottom = 0
                top_layer_handle = None
            else:
                Jbottom = self.J[i - 1]
                bottom_layer_handle = self.layers[i - 1]
            if (i + 1) >= len(self.layers):
                Jtop = 0
                top_layer_handle = None
            else:
                Jtop = self.J[i]
                top_layer_handle = self.layers[i + 1]
            self.energy[i] = layer.compute_energy(self.Hext, Jbottom, Jtop,
                                                  top_layer_handle,
                                                  bottom_layer_handle)
            self.current_position[i, :] = layer.get_current_position()
            self.grad_energy[i, :] = layer.compute_grad_energy(self.Hext)

    def gradient_descent(self, max_steps: int, tol: float = 1e-8):
        for i in tqdm(range(int(max_steps)), mininterval=0.5):
            self.compute_energy_step()
            # compute gradient update
            new_position = self.current_position - self.eta * self.grad_energy
            # update layers
            for l, layer in enumerate(self.layers):
                layer.set_current_position(new_position[l])
            pos_difference = np.abs(self.current_position - new_position).sum()
            if ((pos_difference < tol)
                    or (np.linalg.norm(self.grad_energy) < tol)):
                break
            self.current_position = new_position

        print(f"Gradient descent finished in {i} steps")
        print(f"Final energy: {self.energy.sum()*1e6:.4f} uJ")
        print(f"Final position: {np.rad2deg(self.current_position)}")


# # Equilibrum
# import sympy as sym
# import mpmath

# mpmath.mp.dps = 10
# # define symbols

# theta, phi, K1, K2, Ms = sym.symbols(r'\theta \phi K_1 K_2, M_s')
# Hx, Hy, Hz = sym.symbols(r'H_x H_y H_z')
# Nx, Ny, Nz = sym.symbols(r'N_x N_y N_z')
# ctheta = sym.cos(theta)
# stheta = sym.sin(theta)
# cphi = sym.cos(phi)
# sphi = sym.sin(phi)

# eq1 = -Hx * ctheta * cphi - Hy * ctheta * sphi + Hz * sphi + (
#     Ms * Nx * cphi**2 + Ms * Ny * sphi**2 - Ms * Nz) * ctheta * stheta + (
#         (2 * K1 / Ms) *
#         (2 * ctheta * (stheta**3) * (cphi**2) * (sphi**2) + ctheta * stheta *
#          (ctheta**2 - stheta**2))) + (2 * (K2 / Ms) * ctheta * (stheta**3) *
#                                       (cphi**2) * (sphi**2) *
#                                       (2 * ctheta**2 - stheta**2))

# eq2 = Hx * sphi - Hy * cphi + (Ms * Ny - Ms * Nx) * stheta * cphi * sphi + (
#     (2 * K1 / Ms) * (stheta**3) * cphi * sphi * (cphi**2 - sphi**2) +
#     (2 * K2 / Ms) * (ctheta**2) * (stheta**2) * cphi * sphi *
#     (cphi**2 - sphi**2))

# m_subs = {
#     K1: 1.05e3,
#     K2: 10,
#     Ms: 1.65,
#     Hx: 0,
#     Hy: 300e3,
#     Hz: 0,
#     Nx: 0,
#     Ny: 0,
#     Nz: 1
# }
# theta_init = np.pi / 2
# phi_init = 0

# seq1 = eq1.subs(m_subs)
# seq2 = eq2.subs(m_subs)

# result = sym.solvers.nsolve((seq1, seq2), (theta, phi), (theta_init, phi_init))

# resonance_equation1 = (Hx * stheta * ctheta + Hy * stheta * sphi +
#                        Hz * ctheta +
#                        (Ms * Nx * ctheta**2 + Ms * Ny * sphi**2 - Ms * Nz) *
#                        (ctheta**2 - stheta**2) + (2 * K1 / Ms) *
#                        (ctheta**4 + (stheta**4) * (cphi**4 + sphi**4) - 3 *
#                         (ctheta**2) * (stheta**2) * (1 + cphi**4 + sphi**4)) +
#                        (2 * K2 / Ms) * ((stheta**2) * (cphi**2) * (sphi**2) *
#                                         (6 * ctheta**4 + stheta**4 - 11 *
#                                          (ctheta**2) * (stheta**2))))
# resonance_equation2 = (Hx * stheta * cphi + Hy * stheta * sphi + Hz * ctheta +
#                        Ms * Nx * (sphi**2 - (stheta**2) * (cphi**2)) +
#                        Ms * Ny * (cphi**2 - (stheta**2) * (sphi**2)) -
#                        Ms * Nz * ctheta**2 + (2 * K1 / Ms) *
#                        (ctheta**4 + (stheta**4) * (cphi**4 - sphi**4) - 6 *
#                         (stheta**2) * (cphi**2) * (sphi**2)) + (2 * K2 / Ms) *
#                        ((ctheta**2) * (stheta**2) *
#                         (cphi**4 + sphi**4 - (4 + 3 * stheta**2) * (cphi**2) *
#                          (sphi**2))))
# main_resonance_equation = resonance_equation1 * resonance_equation2
# resonance = main_resonance_equation.subs(m_subs).subs({
#     theta: float(result[0]),
#     phi: float(result[1])
# })
# # omega = 2*pi*f
# gamma = 28024e6
# f = np.sqrt(float(resonance))*gamma
