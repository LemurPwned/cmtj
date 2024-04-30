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

from ..utils import VectorObj, gamma, gamma_rad, mu0, perturb_position
from ..utils.solvers import RootFinder

EPS = np.finfo("float64").resolution


def real_deocrator(fn):
    """Using numpy real cast is way faster than sympy."""

    def wrap_fn(*args):
        return np.real(fn(*args))

    return wrap_fn


@njit
def fast_norm(x):  # sourcery skip: for-index-underscore, sum-comprehension
    """Fast norm function for 1D arrays."""
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
        all_symbols.extend(
            (sym.Symbol(r"\theta_" + indx_i), sym.Symbol(r"\phi_" + indx_i)))
    energy_functional_expr = sym.Function("E")(*all_symbols)
    return (
        get_hessian_from_energy_expr(N, energy_functional_expr),
        energy_functional_expr,
    )


@lru_cache
def get_hessian_from_energy_expr(N: int, energy_functional_expr: sym.Expr):
    """
    Computes the Hessian matrix of the energy functional expression with respect to the spin angles and phases.

    :param N (int): The number of spins.
    :param energy_functional_expr (sympy.Expr): The energy functional expression.

    returns: sympy.Matrix: The Hessian matrix of the energy functional expression.
    """
    hessian = [[0 for _ in range(2 * N)] for _ in range(2 * N)]
    for i in range(N):
        # indx_i = str(i + 1) # for display purposes
        indx_i = str(i)
        # z = sym.Symbol("Z")
        # these here must match the Ms symbols!
        z = (sym.Symbol(r"\omega") * sym.Symbol(r"M_{" + indx_i + "}") *
             sym.sin(sym.Symbol(r"\theta_" + indx_i)) *
             sym.Symbol(r"t_{" + indx_i + "}"))
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

            s1 = sym.Symbol(r"\theta_" + indx_i)
            s2 = sym.Symbol(r"\phi_" + indx_j)
            expr = sym.diff(energy_functional_expr, s1, s2)
            if i == j:
                hessian[2 * i + 1][2 * j] = expr + sym.I * z
                hessian[2 * i][2 * j + 1] = expr - sym.I * z
            else:
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
                               subs=None):
    """Get all second derivatives of the energy expression.
    :param energy_functional_expr: symbolic energy_functional expression
    :param energy_expression: symbolic energy expression (from solver)
    :param subs: substitutions to be made."""
    if subs is None:
        subs = {}
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
        self.m = sym.ImmutableMatrix([
            sym.sin(self.theta) * sym.cos(self.phi),
            sym.sin(self.theta) * sym.sin(self.phi),
            sym.cos(self.theta),
        ])

    def get_coord_sym(self):
        """Returns the symbolic coordinates of the layer."""
        return self.theta, self.phi

    def get_m_sym(self):
        """Returns the magnetisation vector."""
        return self.m

    @lru_cache(3)
    def symbolic_layer_energy(
        self,
        H: sym.ImmutableMatrix,
        J1top: float,
        J1bottom: float,
        J2top: float,
        J2bottom: float,
        top_layer: "LayerSB",
        down_layer: "LayerSB",
    ):
        """Returns the symbolic expression for the energy of the layer.
        Coupling contribution comes only from the bottom layer (top-down crawl)"""
        m = self.get_m_sym()

        eng_non_interaction = self.no_iec_symbolic_layer_energy(H)

        top_iec_energy = 0
        bottom_iec_energy = 0

        if top_layer is not None:
            other_m = top_layer.get_m_sym()
            top_iec_energy = (-(J1top / self.thickness) * m.dot(other_m) -
                              (J2top / self.thickness) * m.dot(other_m)**2)
        if down_layer is not None:
            other_m = down_layer.get_m_sym()
            bottom_iec_energy = (
                -(J1bottom / self.thickness) * m.dot(other_m) -
                (J2bottom / self.thickness) * m.dot(other_m)**2)
        return eng_non_interaction + top_iec_energy + bottom_iec_energy

    def no_iec_symbolic_layer_energy(self, H: sym.ImmutableMatrix):
        """Returns the symbolic expression for the energy of the layer.
        Coupling contribution comes only from the bottom layer (top-down crawl)"""
        m = self.get_m_sym()

        alpha = sym.ImmutableMatrix(
            [sym.cos(self.Kv.phi),
             sym.sin(self.Kv.phi), 0])

        field_energy = -mu0 * self.Ms * m.dot(H)
        surface_anistropy = (-self.Ks +
                             (1.0 / 2.0) * mu0 * self.Ms**2) * (m[-1]**2)
        volume_anisotropy = -self.Kv.mag * (m.dot(alpha)**2)
        return field_energy + surface_anistropy + volume_anisotropy

    def sb_correction(self):
        omega = sym.Symbol(r"\omega")
        return (omega / gamma) * self.Ms * sym.sin(self.theta) * self.thickness

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, __value: "LayerSB") -> bool:
        return (self._id == __value._id and self.thickness == __value.thickness
                and self.Kv == __value.Kv and self.Ks == __value.Ks
                and self.Ms == __value.Ms)


@dataclass
class LayerDynamic(LayerSB):
    alpha: float

    def rhs_llg(
        self,
        H: sym.Matrix,
        J1top: float,
        J1bottom: float,
        J2top: float,
        J2bottom: float,
        top_layer: "LayerSB",
        down_layer: "LayerSB",
    ):
        """Returns the symbolic expression for the RHS of the spherical LLG equation.
        Coupling contribution comes only from the bottom layer (top-down crawl)"""
        U = self.symbolic_layer_energy(
            H,
            J1top=J1top,
            J1bottom=J1bottom,
            J2top=J2top,
            J2bottom=J2bottom,
            top_layer=top_layer,
            down_layer=down_layer,
        )
        # sum all components
        prefac = gamma_rad / (1.0 + self.alpha)**2
        inv_sin = 1.0 / (sym.sin(self.theta) + EPS)
        dUdtheta = sym.diff(U, self.theta)
        dUdphi = sym.diff(U, self.phi)

        dtheta = -inv_sin * dUdphi - self.alpha * dUdtheta
        dphi = inv_sin * dUdtheta - self.alpha * dUdphi * (inv_sin)**2
        return prefac * sym.ImmutableMatrix([dtheta, dphi]) / self.Ms

    def __eq__(self, __value: "LayerDynamic") -> bool:
        return super().__eq__(__value) and self.alpha == __value.alpha

    def __hash__(self) -> int:
        return super().__hash__()


@dataclass
class Solver:
    """General solver for the system.

    :param layers: list of layers in the system.
    :param J1: list of interlayer exchange constants. Goes (i)-(i+1), i = 0, 1, 2, ...
        with i being the index of the layer.
    :param J2: list of interlayer exchange constants.
    :param H: external field.
    :param Ndipole: list of dipole fields for each layer. Defaults to None.
        Goes (i)-(i+1), i = 0, 1, 2, ... with i being the index of the layer.
    """

    layers: List[Union[LayerSB, LayerDynamic]]
    J1: List[float]
    J2: List[float]
    H: VectorObj = None
    Ndipole: List[List[VectorObj]] = None

    def __post_init__(self):
        if len(self.layers) != len(self.J1) + 1:
            raise ValueError("Number of layers must be 1 more than J1.")
        if len(self.layers) != len(self.J2) + 1:
            raise ValueError("Number of layers must be 1 more than J2.")

        self.dipoleMatrix: list[sym.Matrix] = None
        if self.Ndipole is not None:
            if len(self.layers) != len(self.Ndipole) + 1:
                raise ValueError(
                    "Number of layers must be 1 more than number of tensors.")
            if isinstance(self.layers[0], LayerDynamic):
                raise ValueError(
                    "Dipole coupling is not yet supported for LayerDynamic.")
            self.dipoleMatrix = [
                sym.Matrix([d.get_cartesian() for d in dipole])
                for dipole in self.Ndipole
            ]

        id_sets = {layer._id for layer in self.layers}
        ideal_set = set(range(len(self.layers)))
        if id_sets != ideal_set:
            raise ValueError("Layer ids must be 0, 1, 2, ... and unique."
                             "Ids must start from 0.")

    def get_layer_references(self, layer_indx, interaction_constant):
        """Returns the references to the layers above and below the layer
        with index layer_indx."""
        if len(self.layers) == 1:
            return None, None, 0, 0
        if layer_indx == 0:
            return None, self.layers[layer_indx +
                                     1], 0, interaction_constant[0]
        elif layer_indx == len(self.layers) - 1:
            return self.layers[layer_indx -
                               1], None, interaction_constant[-1], 0
        return (
            self.layers[layer_indx - 1],
            self.layers[layer_indx + 1],
            interaction_constant[layer_indx - 1],
            interaction_constant[layer_indx],
        )

    def compose_llg_jacobian(self, H: VectorObj):
        """Create a symbolic jacobian of the LLG equation in spherical coordinates."""
        # has order theta0, phi0, theta1, phi1, ...
        if isinstance(H, VectorObj):
            H = sym.ImmutableMatrix(H.get_cartesian())

        symbols, fns = [], []
        for i, layer in enumerate(self.layers):
            symbols.extend((layer.theta, layer.phi))
            top_layer, bottom_layer, Jtop, Jbottom = self.get_layer_references(
                i, self.J1)
            _, _, J2top, J2bottom = self.get_layer_references(i, self.J2)
            fns.append(
                layer.rhs_llg(H, Jtop, Jbottom, J2top, J2bottom, top_layer,
                              bottom_layer))
        jac = sym.ImmutableMatrix(fns).jacobian(symbols)
        return jac, symbols

    def create_energy(self,
                      H: Union[VectorObj, sym.ImmutableMatrix] = None,
                      volumetric: bool = False):
        """Creates the symbolic energy expression.

        Due to problematic nature of coupling, there is an issue of
        computing each layer's FMR in the presence of IEC.
        If volumetric = True then we use the thickness of the layer to multiply the
        energy and hence avoid having to divide J by the thickness of a layer.
        If volumetric = False the J constant is divided by weighted thickness
        and included in every layer's energy, correcting FMR automatically.
        """
        if H is None:
            h = self.H.get_cartesian()
            H = sym.ImmutableMatrix(h)
        energy = 0
        if volumetric:
            # volumetric energy -- DO NOT USE IN GENERAL
            for i, layer in enumerate(self.layers):
                top_layer, bottom_layer, Jtop, Jbottom = self.get_layer_references(
                    i, self.J1)
                _, _, J2top, J2bottom = self.get_layer_references(i, self.J2)
                ratio_top, ratio_bottom = 0, 0
                if top_layer:
                    ratio_top = top_layer.thickness / (top_layer.thickness +
                                                       layer.thickness)
                if bottom_layer:
                    ratio_bottom = bottom_layer.thickness / (
                        layer.thickness + bottom_layer.thickness)
                energy += layer.symbolic_layer_energy(
                    H,
                    Jtop * ratio_top,
                    Jbottom * ratio_bottom,
                    J2top,
                    J2bottom,
                    top_layer,
                    bottom_layer,
                )
        else:
            # surface energy for correct angular gradient
            for layer in self.layers:
                # to avoid dividing J by thickness
                energy += layer.no_iec_symbolic_layer_energy(
                    H) * layer.thickness

            for i in range(len(self.layers) - 1):
                l1m = self.layers[i].get_m_sym()
                l2m = self.layers[i + 1].get_m_sym()
                ldot = l1m.dot(l2m)
                energy -= self.J1[i] * ldot
                energy -= self.J2[i] * (ldot)**2

                # dipole fields
                if self.dipoleMatrix is not None:
                    mat = self.dipoleMatrix[i]
                    # is positive, just like demag
                    energy += ((mu0 / 2.0) * l1m.dot(mat * l2m) *
                               self.layers[i].Ms * self.layers[i + 1].Ms *
                               self.layers[i].thickness)
                    energy += ((mu0 / 2.0) * l2m.dot(mat * l1m) *
                               self.layers[i].Ms * self.layers[i + 1].Ms *
                               self.layers[i + 1].thickness)
        return energy

    def create_energy_hessian(self, equilibrium_position: List[float]):
        """Creates the symbolic hessian of the energy expression."""
        energy = self.create_energy(volumetric=False)
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

                expr = sym.diff(sym.diff(energy, theta_i), phi_j)
                # mixed terms
                if i == j:
                    hessian[2 * i + 1][2 * j] = expr + sym.I * z
                    hessian[2 * i][2 * j + 1] = expr - sym.I * z
                else:
                    hessian[2 * i][2 * j + 1] = expr
                    hessian[2 * j + 1][2 * i] = expr

                    expr = sym.diff(sym.diff(energy, phi_i), theta_j)
                    hessian[2 * i + 1][2 * j] = expr
                    hessian[2 * j][2 * i + 1] = expr

        hes = sym.ImmutableMatrix(hessian)
        _, U, _ = hes.LUdecomposition()
        return U.det().subs(subs)

    def get_gradient_expr(self, accel="math"):
        """Returns the symbolic gradient of the energy expression."""
        energy = self.create_energy(volumetric=False)
        grad_vector = []
        symbols = []
        for layer in self.layers:
            (theta, phi) = layer.get_coord_sym()
            grad_vector.extend((sym.diff(energy, theta), sym.diff(energy,
                                                                  phi)))
            symbols.extend((theta, phi))
        return sym.lambdify(symbols, grad_vector, accel)

    def adam_gradient_descent(
        self,
        init_position: np.ndarray,
        max_steps: int,
        tol: float = 1e-8,
        learning_rate: float = 1e-4,
        first_momentum_decay: float = 0.9,
        second_momentum_decay: float = 0.999,
        perturbation: float = 1e-6,
    ):
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
        if perturbation:
            current_position = perturb_position(init_position, perturbation)
        m = np.zeros_like(current_position)  # first momentum
        v = np.zeros_like(current_position)  # second momentum
        eps = 1e-12
        # history = []
        while True:
            step += 1
            grad = np.asarray(gradfn(*current_position))
            m = first_momentum_decay * m + (1.0 - first_momentum_decay) * grad
            v = second_momentum_decay * v + (1.0 -
                                             second_momentum_decay) * grad**2
            m_hat = m / (1.0 - first_momentum_decay**step)
            v_hat = v / (1.0 - second_momentum_decay**step)
            new_position = current_position - learning_rate * m_hat / (
                np.sqrt(v_hat) + eps)
            if step > max_steps:
                break
            if fast_norm(current_position - new_position) < tol:
                break
            current_position = new_position
            # history.append(current_position)
        # return np.asarray(current_position), np.asarray(history)
        return np.asarray(current_position)

    def single_layer_resonance(self, layer_indx: int, eq_position: np.ndarray):
        """We can compute the equilibrium position of a single layer directly.
        :param layer_indx: the index of the layer to compute the equilibrium
        :param eq_position: the equilibrium position vector"""
        layer = self.layers[layer_indx]
        theta_eq = eq_position[2 * layer_indx]
        theta, phi = self.layers[layer_indx].get_coord_sym()
        energy = self.create_energy(volumetric=True)
        subs = self.get_subs(eq_position)
        d2Edtheta2 = sym.diff(sym.diff(energy, theta), theta).subs(subs)
        d2Edphi2 = sym.diff(sym.diff(energy, phi), phi).subs(subs)
        # mixed, assuming symmetry
        d2Edthetaphi = sym.diff(sym.diff(energy, theta), phi).subs(subs)
        vareps = 1e-18

        fmr = (d2Edtheta2 * d2Edphi2 - d2Edthetaphi**2) / np.power(
            np.sin(theta_eq + vareps) * layer.Ms, 2)
        fmr = np.sqrt(float(fmr)) * gamma_rad / (2 * np.pi)
        return fmr

    def solve(
        self,
        init_position: np.ndarray,
        max_steps: int = 1e9,
        learning_rate: float = 1e-4,
        adam_tol: float = 1e-8,
        first_momentum_decay: float = 0.9,
        second_momentum_decay: float = 0.999,
        perturbation: float = 1e-3,
        ftol: float = 0.01e9,
        max_freq: float = 80e9,
        force_single_layer: bool = False,
        force_sb: bool = False,
    ):
        """Solves the system.
        For dynamic LayerDynamic, the return is different, check :return.
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
        :param adam_tol: tolerance for the consecutive Adam minima.
        :param first_momentum_decay: constant for the first momentum.
        :param second_momentum_decay: constant for the second momentum.
        :param perturbation: the perturbation to use for the numerical gradient computation.
        :param ftol: tolerance for the frequency search. [numerical only]
        :param max_freq: maximum frequency to search for. [numerical only]
        :param force_single_layer: whether to force the computation of the frequencies
                                   for each layer individually.
        :param force_sb: whether to force the computation of the frequencies.
                        Takes effect only if the layers are LayerDynamic, not LayerSB.
        :return: equilibrium position and frequencies in [GHz] (and eigenvectors if LayerDynamic instead of LayerSB).
        """
        if self.H is None:
            raise ValueError(
                "H must be set before solving the system numerically.")
        eq = self.adam_gradient_descent(
            init_position=init_position,
            max_steps=max_steps,
            tol=adam_tol,
            learning_rate=learning_rate,
            first_momentum_decay=first_momentum_decay,
            second_momentum_decay=second_momentum_decay,
            perturbation=perturbation,
        )
        if not force_sb and isinstance(self.layers[0], LayerDynamic):
            eigenvalues, eigenvectors = self.dynamic_layer_solve(eq)
            return eq, eigenvalues / 1e9, eigenvectors
        N = len(self.layers)
        if N == 1:
            return eq, [self.single_layer_resonance(0, eq) / 1e9]
        if force_single_layer:
            frequencies = []
            for indx in range(N):
                frequency = self.single_layer_resonance(indx, eq) / 1e9
                frequencies.append(frequency)
            return eq, frequencies
        return self.num_solve(eq, ftol=ftol, max_freq=max_freq)

    def dynamic_layer_solve(self, eq: List[float]):
        """Return the FMR frequencies and modes for N layers using the
        dynamic RHS model
        :param eq: the equilibrium position of the system.
        :return: frequencies and eigenmode vectors."""
        jac, symbols = self.compose_llg_jacobian(self.H)
        subs = {symbols[i]: eq[i] for i in range(len(eq))}
        jac = jac.subs(subs)
        jac = np.asarray(jac, dtype=np.float32)
        eigvals, eigvecs = np.linalg.eig(jac)
        eigvals_im = np.imag(eigvals) / (2 * np.pi)
        indx = np.argwhere(eigvals_im > 0).ravel()
        return eigvals_im[indx], eigvecs[indx]

    def num_solve(self,
                  eq: List[float],
                  ftol: float = 0.01e9,
                  max_freq: float = 80e9):
        hes = self.create_energy_hessian(eq)
        omega = sym.Symbol(r"\omega")
        if len(self.layers) <= 3:
            y = real_deocrator(njit(sym.lambdify(omega, hes, "math")))
        else:
            y = real_deocrator(sym.lambdify(omega, hes, "math"))
        r = RootFinder(0, max_freq, step=ftol, xtol=1e-8, root_dtype="float16")
        roots = r.find(y)
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
        system_energy = self.create_energy(H=Hsym, volumetric=False)
        root_expr, energy_functional_expr = find_analytical_roots(N)
        subs = get_all_second_derivatives(energy_functional_expr,
                                          energy_expression=system_energy,
                                          subs={})
        subs.update(self.get_ms_subs())
        return [s.subs(subs) for s in root_expr]

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
        a = {r"M_{" + str(layer._id) + "}": layer.Ms for layer in self.layers}
        b = {
            r"t_{" + str(layer._id) + r"}": layer.thickness
            for layer in self.layers
        }
        return a | b

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
        disable_tqdm: bool = False,
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
