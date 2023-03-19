import math
import time
import warnings
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, List, Tuple, Union

import numpy as np
import sympy as sym
from numba import njit
from tqdm import tqdm

from cmtj.utils import VectorObj, gamma, gamma_rad, mu0

from ..utils.solvers import RootFinder


def real_deocrator(fn):

    def wrap_fn(*args):
        return np.real(fn(*args))

    return wrap_fn


@njit
def fast_norm(x):
    """Fast norm function for !D arrays.
    """
    sum_ = 0
    for x_ in x:
        sum_ += x_**2
    return math.sqrt(sum_)


@lru_cache
def general_hessian_functional(N: int):
    """Create a generalised hessian functional for N layers.
    :param N: number of layers.
    ! WARNING: remember Ms Symbols must match!!!
    """
    all_symbols = []
    for i in range(N):
        # indx_i = str(i + 1) # for display purposes
        indx_i = str(i)
        all_symbols.append(sym.Symbol(r"\theta_" + indx_i))
        all_symbols.append(sym.Symbol(r"\phi_" + indx_i))
    energy_functional_expr = sym.Function("E")(*all_symbols)
    return get_hessian_from_energy_expr(
        N, energy_functional_expr), energy_functional_expr


@lru_cache
def get_hessian_from_energy_expr(N: int, energy_functional_expr: sym.Expr):
    hessian = [[0 for _ in range(2 * N)] for _ in range(2 * N)]
    for i in range(N):
        # indx_i = str(i + 1) # for display purposes
        indx_i = str(i)
        # z = sym.Symbol("Z")

        # these here must match the Ms symbols!
        z = sym.Symbol(r"\omega") * sym.Symbol(
            r"M_{" + indx_i + "}") * sym.sin(sym.Symbol(r"\theta_" + indx_i))
        for j in range(i, N):
            # indx_j = str(j + 1) # for display purposes
            indx_j = str(j)
            s1 = sym.Symbol(r"\theta_" + indx_i)
            s2 = sym.Symbol(r"\theta_" + indx_j)

            expr = sym.diff(energy_functional_expr, s1, s2)
            hessian[2 * i][2 * j] = expr
            hessian[2 * j][2 * i] = expr
            s1 = sym.Symbol(r"\phi_" + indx_i)
            s2 = sym.Symbol(r"\phi_" + indx_j)
            expr = sym.diff(energy_functional_expr, s1, s2)
            hessian[2 * i + 1][2 * j + 1] = expr
            hessian[2 * j + 1][2 * i + 1] = expr

            if i == j:
                s1 = sym.Symbol(r"\theta_" + indx_i)
                s2 = sym.Symbol(r"\phi_" + indx_j)
                expr = sym.diff(energy_functional_expr, s1, s2)
                hessian[2 * i + 1][2 * j] = expr + sym.I * z
                hessian[2 * i][2 * j + 1] = expr - sym.I * z
            else:
                s1 = sym.Symbol(r"\theta_" + indx_i)
                s2 = sym.Symbol(r"\phi_" + indx_j)
                expr = sym.diff(energy_functional_expr, s1, s2)
                hessian[2 * i][2 * j + 1] = expr
                hessian[2 * j + 1][2 * i] = expr

                s1 = sym.Symbol(r"\phi_" + indx_i)
                s2 = sym.Symbol(r"\theta_" + indx_j)
                expr = sym.diff(energy_functional_expr, s1, s2)
                hessian[2 * i + 1][2 * j] = expr
                hessian[2 * j][2 * i + 1] = expr
    return sym.Matrix(hessian)


@lru_cache()
def solve_for_determinant(N: int):
    """Solve for the determinant of the hessian functional.
    :param N: number of layers.
    """
    hessian, energy_functional_expr = general_hessian_functional(N)
    if N == 1:
        return hessian.det(), energy_functional_expr

    # LU decomposition
    _, U, _ = hessian.LUdecomposition()
    return U.det(), energy_functional_expr


@lru_cache
def find_analytical_roots(N: int):
    det_expr, energy_functional_expr = solve_for_determinant(N)
    solutions = sym.solve(sym.simplify(det_expr), sym.Symbol(r"\omega"))
    return solutions, energy_functional_expr


def get_all_second_derivatives(energy_functional_expr,
                               energy_expression,
                               subs={}):
    """Get all second derivatives of the energy expression.
    :param energy_functional_expr: symbolic energy_functional expression
    :param energy_expression: symbolic energy expression (from solver)
    :param subs: substitutions to be made."""
    second_derivatives = subs
    symbols = energy_expression.free_symbols
    for i, s1 in enumerate(symbols):
        for j, s2 in enumerate(symbols):
            if i <= j:
                org_diff = sym.diff(energy_functional_expr, s1, s2)
                if subs is not None:
                    second_derivatives[org_diff] = sym.diff(
                        energy_expression, s1, s2).subs(subs)
                else:
                    second_derivatives[org_diff] = sym.diff(
                        energy_expression, s1, s2)
    return second_derivatives


@dataclass
class LayerSB:
    """Basic Layer for Smit-Beljers model.
    :param thickness: thickness of the FM layer (effective).
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
        self.theta = sym.Symbol(r"\theta_" + str(self._id))
        self.phi = sym.Symbol(r"\phi_" + str(self._id))
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
        Coupling contribution comes only from the bottom layer (top-down crawl)"""
        m = self.get_m_sym()

        alpha = sym.Matrix([sym.cos(self.Kv.phi), sym.sin(self.Kv.phi), 0])

        field_energy = -mu0 * self.Ms * m.dot(H)
        surface_anistropy = (-self.Ks +
                             (1. / 2.) * mu0 * self.Ms**2) * (m[-1]**2)
        volume_anisotropy = -self.Kv.mag * (m.dot(alpha)**2)

        if down_layer is None:
            iec_energy = 0
        else:
            other_m = down_layer.get_m_sym()
            iec_energy = -(J1 / self.thickness) * m.dot(other_m) - (
                J2 / self.thickness) * m.dot(other_m)**2

        return field_energy + surface_anistropy + volume_anisotropy + iec_energy

    def sb_correction(self):
        omega = sym.Symbol(r'\omega')
        # we use unreduced gamma
        # later we don't need to divide omega by 2pi
        # idk but if we use reduced gamma it numerically breaks lol
        Z = (omega / gamma) * self.Ms * sym.sin(self.theta)
        return Z


@dataclass
class SolverSB:
    layers: List[LayerSB]
    J1: List[float]
    J2: List[float]
    H: VectorObj = None

    def __post_init__(self):
        if len(self.layers) != len(self.J1) + 1:
            raise ValueError("Number of layers must be 1 more than J1.")
        if len(self.layers) != len(self.J2) + 1:
            raise ValueError("Number of layers must be 1 more than J2.")

        # TODO: check the ids of the layers
        id_sets = set([layer._id for layer in self.layers])
        ideal_set = set(range(len(self.layers)))
        if id_sets != ideal_set:
            raise ValueError("Layer ids must be 0, 1, 2, ... and unique")

    def create_energy(self, H: Union[VectorObj, sym.Matrix] = None):
        """Creates the symbolic energy expression."""
        if H is None:
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
        energy = self.create_energy()
        subs = self.get_subs(equilibrium_position)
        N = len(self.layers)
        hessian = [[0 for _ in range(2 * N)] for _ in range(2 * N)]
        for i in range(N):
            z = self.layers[i].sb_correction()
            theta_i, phi_i = self.layers[i].get_coord_sym()
            for j in range(i, N):
                # dtheta dtheta
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
                    expr = sym.diff(sym.diff(energy, theta_i), phi_j)
                    hessian[2 * i][2 * j + 1] = expr
                    hessian[2 * j + 1][2 * i] = expr

                    expr = sym.diff(sym.diff(energy, phi_i), theta_j)
                    hessian[2 * i + 1][2 * j] = expr
                    hessian[2 * j][2 * i + 1] = expr

        hes = sym.Matrix(hessian)
        _, U, _ = hes.LUdecomposition()
        return U.det().subs(subs)

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
            # if np.linalg.norm(current_position - new_position) < tol:
            # break
            if fast_norm(current_position - new_position) < tol:
                break
            current_position = new_position
        return np.asarray(current_position)

    def single_layer_equilibrium(self, layer_indx: int,
                                 eq_position: np.ndarray):
        """We can compute the equilibrium position of a single layer directly.
        :param layer_indx: the index of the layer to compute the equilibrium
        :param eq_position: the equilibrium position vector"""
        # compute grads
        theta_eq, phi_eq = eq_position
        layer = self.layers[layer_indx]
        theta, phi = self.layers[layer_indx].get_coord_sym()
        energy = self.create_energy()
        subs = {theta: theta_eq, phi: phi_eq}
        d2Edtheta2 = sym.diff(sym.diff(energy, theta), theta).subs(subs)
        d2Edphi2 = sym.diff(sym.diff(energy, phi), phi).subs(subs)
        # mixed, assuming symmetry
        d2Edthetaphi = sym.diff(sym.diff(energy, theta), phi).subs(subs)
        vareps = 1e-12

        fmr = (d2Edtheta2 * d2Edphi2 - d2Edthetaphi**2) / np.power(
            np.sin(theta_eq + vareps) * layer.Ms, 2)
        fmr = np.sqrt(float(fmr)) * gamma_rad / (2 * np.pi)
        return fmr

    def solve(self,
              init_position: np.ndarray,
              max_steps: int = 1e9,
              learning_rate: float = 1e-4,
              first_momentum_decay: float = 0.9,
              second_momentum_decay: float = 0.999,
              ftol: float = 0.01e9,
              max_freq: float = 80e9):
        """Solves the system.
        1. Computes the energy functional.
        2. Computes the gradient of the energy functional.
        3. Performs a gradient descent to find the equilibrium position.
        Returns the equilibrium position and frequencies in [GHz].
        If there's only one layer, the frequency is computed analytically.
        For full analytical solution, see: `analytical_field_scan`
        :param init_position: initial position for the gradient descent.
                              Must be a 1D array of size 2 * number of layers (theta, phi)
        :param max_steps: maximum number of gradient steps.
        :param learning_rate: the learning rate (descent speed).
        :param first_momentum_decay: constant for the first momentum.
        :param second_momentum_decay: constant for the second momentum.
        :param ftol: tolerance for the frequency search. [numerical only]
        :param max_freq: maximum frequency to search for. [numerical only]
        """
        if self.H is None:
            raise ValueError(
                "H must be set before solving the system numerically.")
        eq = self.adam_gradient_descent(
            init_position=init_position,
            max_steps=max_steps,
            tol=1e-9,
            learning_rate=learning_rate,
            first_momentum_decay=first_momentum_decay,
            second_momentum_decay=second_momentum_decay,
        )
        if len(self.layers) == 1:
            return eq, self.single_layer_equilibrium(0, eq) / 1e9

        return self.num_solve(eq, ftol=ftol, max_freq=max_freq)

    def num_solve(self,
                  eq: List[float],
                  ftol: float = 0.01e9,
                  max_freq: float = 80e9):
        hes = self.create_energy_hessian(eq)
        omega = sym.Symbol(r"\omega")

        y = real_deocrator(njit(sym.lambdify(omega, hes, 'math')))
        r = RootFinder(0, max_freq, step=ftol, xtol=1e-8, root_dtype="float16")
        roots = r.find(y)
        # print("Root finding took {:.2f} seconds".format(t2 - t1))
        # convert to GHz
        # reduce unique solutions to 2 decimal places
        # don't divide by 2pi, we used gamma instead of gamma / 2pi
        f = np.unique(np.around(roots / 1e9, 2))
        return eq, f

    def analytical_roots(self):
        """Find & cache the analytical roots of the system.
        Returns a list of solutions.
        Ineffecient for more than 2 layers (can try though).
        """
        Hsym = sym.Matrix([
            sym.Symbol(r"H_{x}"),
            sym.Symbol(r"H_{y}"),
            sym.Symbol(r"H_{z}"),
        ])
        N = len(self.layers)
        if N > 2:
            warnings.warn(
                "Analytical solutions for over 2 layers may be computationally expensive."
            )
        system_energy = self.create_energy(H=Hsym)
        root_expr, energy_functional_expr = find_analytical_roots(N)
        subs = get_all_second_derivatives(energy_functional_expr,
                                          energy_expression=system_energy,
                                          subs={})
        subs.update(self.get_ms_subs())
        # substitute all known values
        solutions = [s.subs(subs) for s in root_expr]
        return solutions

    def get_subs(self, equilibrium_position: List[float]):
        """Returns the substitution dictionary for the energy expression."""
        subs = {}
        for i in range(len(self.layers)):
            theta, phi = self.layers[i].get_coord_sym()
            subs[theta] = equilibrium_position[2 * i]
            subs[phi] = equilibrium_position[(2 * i) + 1]
        return subs

    def get_ms_subs(self):
        """Returns a dictionary of substitutions for the Ms symbols."""
        return {
            r"M_{" + str(layer._id) + "}": layer.Ms
            for layer in self.layers
        }

    def set_H(self, H: VectorObj):
        """Sets the external field."""
        self.H = H

    def analytical_field_scan(
        self,
        Hrange: List[VectorObj],
        init_position: List[float] = None,
        max_steps: int = 1e9,
        learning_rate: float = 1e-4,
        first_momentum_decay: float = 0.9,
        second_momentum_decay: float = 0.999,
        disable_tqdm: bool = False
    ) -> Iterable[Tuple[List[float], List[float], VectorObj]]:
        """Performs a field scan using the analytical solutions.
        :param Hrange: the range of fields to scan.
        :param init_position: the initial position for the gradient descent.
                              If None, the first field in Hrange will be used.
        :param max_steps: maximum number of gradient steps.
        :param learning_rate: the learning rate (descent speed).
        :param first_momentum_decay: constant for the first momentum.
        :param second_momentum_decay: constant for the second momentum.
        :param disable_tqdm: disable the progress bar.
        :return: an iterable of (equilibrium position, frequencies, field)
        """
        s1 = time.time()
        global_roots = self.analytical_roots()
        s2 = time.time()
        if not disable_tqdm:
            print(f"Analytical roots found in {s2 - s1:.2f} seconds.")
        if init_position is None:
            start = Hrange[0]
            start.mag = 1
            init_position = []
            # align with the first field
            for _ in self.layers:
                init_position.extend([start.theta, start.phi])
        Hsym = sym.Matrix([
            sym.Symbol(r"H_{x}"),
            sym.Symbol(r"H_{y}"),
            sym.Symbol(r"H_{z}"),
        ])
        current_position = init_position
        for Hvalue in tqdm(Hrange, disable=disable_tqdm):
            self.set_H(Hvalue)
            hx, hy, hz = Hvalue.get_cartesian()
            eq = self.adam_gradient_descent(
                init_position=current_position,
                max_steps=max_steps,
                tol=1e-9,
                learning_rate=learning_rate,
                first_momentum_decay=first_momentum_decay,
                second_momentum_decay=second_momentum_decay,
            )
            step_subs = self.get_subs(eq)
            step_subs.update(self.get_ms_subs())
            step_subs.update({Hsym[0]: hx, Hsym[1]: hy, Hsym[2]: hz})
            roots = [s.subs(step_subs) for s in global_roots]
            roots = np.asarray(roots, dtype=np.float32) * gamma / 1e9
            yield eq, roots, Hvalue
            current_position = eq
