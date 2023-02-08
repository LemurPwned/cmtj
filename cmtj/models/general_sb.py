from dataclasses import dataclass
from typing import List

import numpy as np
import sympy as sym

from cmtj.utils import VectorObj, gamma, mu0


@dataclass
class LayerSB:
    """Basic Layer for Smit-Beljers model.
    :param thickness: thickness of the FM layer (effective).
    :param m: initial magnetisation vector (Cauchy condition).
    :param Kv: volumetric (in-plane) anisotropy. Only phi and mag count [J/m^3].
    :param Ks: surface anisotropy (out-of plane, or perpendicular) value [J/m^3].
    :param Ms: magnetisation saturation value in [A/m].
    """
    _id: int
    thickness: float
    Kv: VectorObj
    Ks: float
    Ms: float

    def __post_init__(self):
        if self._id > 9:
            raise ValueError("Only up to 10 layers supported.")
        self.theta, self.phi = sym.symbols(rf'theta_{self._id} phi_{self._id}')
        self.m = sym.Matrix([
            sym.sin(self.theta) * sym.cos(self.phi),
            sym.sin(self.theta) * sym.sin(self.phi),
            sym.cos(self.theta)
        ])

    def get_coord_sym(self):
        """Returns the symbolic coordinates of the layer."""
        return self.theta, self.phi

    def get_m_sym(self):
        """Returns the magnetisation vector."""
        return self.m

    def symbolic_layer_energy(self, H: sym.Matrix, J1: float, J2: float,
                              down_layer: "LayerSB"):
        """Returns the symbolic expression for the energy of the layer.
        Only from the bottom layer (top-down crawl)"""
        m = self.get_m_sym()

        alpha = sym.Matrix([sym.cos(self.Kv.phi), sym.sin(self.Kv.phi), 0])

        field_energy = -mu0 * self.Ms * m.dot(H)
        surface_anistropy = (-self.Ks +
                             (1. / 2.) * mu0 * self.Ms**2) * (m[-1]**2)
        volume_anisotropy = -self.Kv.mag * m.dot(alpha)**2

        if down_layer is None:
            iec_energy = 0
        else:
            other_m = down_layer.get_m_sym()
            iec_energy = -(J1 / self.thickness) * m.dot(other_m) - (
                J2 / self.thickness) * m.dot(other_m)**2

        return field_energy + surface_anistropy + volume_anisotropy + iec_energy

    def sb_correction(self):
        omega = sym.Symbol(r'\omega')
        Z = (omega / gamma) * self.Ms * sym.sin(self.theta)
        return Z


@dataclass
class Solver:
    layers: List[LayerSB]
    H: VectorObj
    J1: List[float]
    J2: List[float]

    def create_energy(self):
        """Creates the symbolic energy expression."""
        h = self.H.get_cartesian()
        H = sym.Matrix(h)
        energy = 0
        for i, layer in enumerate(self.layers):
            if i == 0:
                energy += layer.symbolic_layer_energy(H, 0, 0, None)
            else:
                energy += layer.symbolic_layer_energy(H, self.J1[i - 1],
                                                      self.J2[i - 1],
                                                      self.layers[i - 1])
        return energy

    def create_energy_hessian(self, equilibrium_position: List[float]):
        """Creates the symbolic hessian of the energy expression."""
        energy = sym.simplify(self.create_energy())
        symbols = [
            symbol for layer in self.layers
            for symbol in layer.get_coord_sym()
        ]
        subs = {}
        for i in range(len(self.layers)):
            theta, phi = self.layers[i].get_coord_sym()
            subs[theta] = equilibrium_position[2 * i]
            subs[phi] = equilibrium_position[2 * i + 1]

        N = len(self.layers)
        hessian = [[0 for _ in range(2 * N)] for _ in range(2 * N)]
        # energy = sym.Function('E')(*symbols)
        for i in range(N):
            z = self.layers[i].sb_correction()
            # z = sym.Symbol("Z")
            for j in range(i, N):
                # dtheta dtheta
                theta_i, phi_i = self.layers[i].get_coord_sym()
                theta_j, phi_j = self.layers[j].get_coord_sym()

                expr = sym.diff(sym.diff(energy, theta_i), theta_j)
                hessian[2 * i][2 * j] = expr
                hessian[2 * j][2 * i] = expr

                # dphi dphi
                expr = sym.diff(sym.diff(energy, phi_i), phi_j)
                hessian[2 * i + 1][2 * j + 1] = expr
                hessian[2 * j + 1][2 * i + 1] = expr

                # mixed terms
                if i == j:
                    expr = sym.diff(sym.diff(energy, theta_i), phi_j)
                    hessian[2 * i + 1][2 * j] = expr + sym.I * z
                    hessian[2 * i][2 * j + 1] = expr - sym.I * z
                else:
                    s1 = r"\theta_" + str(i + 1)
                    s2 = r"\phi_" + str(j + 1)
                    expr = sym.diff(sym.diff(energy, theta_i), phi_j)
                    hessian[2 * i][2 * j + 1] = expr
                    hessian[2 * j + 1][2 * i] = expr

                    s1 = r"\phi_" + str(i + 1)
                    s2 = r"\theta_" + str(j + 1)
                    expr = sym.diff(sym.diff(energy, phi_i), theta_j)
                    hessian[2 * i + 1][2 * j] = expr
                    hessian[2 * j][2 * i + 1] = expr

        return sym.Matrix(hessian).subs(subs)

    def get_gradient_expr(self, accel="numpy"):
        """Returns the symbolic gradient of the energy expression."""
        energy = self.create_energy()
        grad_vector = []
        symbols = []
        for layer in self.layers:
            (theta, phi) = layer.get_coord_sym()
            grad_vector.append(sym.diff(energy, theta))
            grad_vector.append(sym.diff(energy, phi))
            symbols.append(theta)
            symbols.append(phi)
        return sym.lambdify(symbols, grad_vector, accel)

    def adam_gradient_descent(self,
                              init_position: np.ndarray,
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
        gradfn = self.get_gradient_expr()
        current_position = init_position

        m = np.zeros_like(current_position)  # first momentum
        v = np.zeros_like(current_position)  # second momentum
        eps = 1e-12
        while True:
            step += 1
            grad = np.asarray(gradfn(*current_position))
            m = first_momentum_decay * m + (1. - first_momentum_decay) * grad
            v = second_momentum_decay * v + (1. -
                                             second_momentum_decay) * grad**2
            m_hat = m / (1. - first_momentum_decay**step)
            v_hat = v / (1. - second_momentum_decay**step)
            new_position = current_position - learning_rate * m_hat / (
                np.sqrt(v_hat) + eps)
            if step > max_steps:
                break
            if np.linalg.norm(current_position - new_position) < tol:
                break
            current_position = new_position
        print(f"Converged in {step} steps")
        return np.asarray(current_position)
