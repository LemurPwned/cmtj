import math
from abc import ABC
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, List, Literal

from numba import njit
from scipy.integrate import RK45

from ..utils import bohr_magneton, echarge, gyromagnetic_ratio, hbar, mu0
from .sb import VectorObj

gyro = gyromagnetic_ratio
pi2 = math.pi / 2.


class DW:
    """Initial conditions for the phi of DW equation."""
    NEEL_RIGHT = 0
    NEEL_LEFT = math.pi
    BLOCH_UP = math.pi / 2.
    BLOCH_DOWN = 3. * math.pi / 2.


@njit
def get_pinning_field(X, Ms, pinning, Ly, Lz, V0_pin):
    arg = X * math.pi / pinning
    dVdx = 2 * math.pi * V0_pin * math.sin(arg) * math.cos(arg)
    denom = 2 * mu0 * Ms * Lz * Ly
    return -(1. / denom) * dVdx


@njit
def get_edge_field(X, Lx, V0_edge):
    c = Lx / 2
    arg = (X - (Lx / 2)) / c
    p = 6
    return -p * V0_edge * (math.sinh(arg) * math.cosh(arg)**(p - 1) / c)


@njit
def get_field_contribution(X, phi, hx, hy, hz, alpha, dw, Ms, V0_pin, pinning,
                           Ly, Lz):

    pinning = get_pinning_field(X,
                                Ms=Ms,
                                pinning=pinning,
                                Ly=Ly,
                                Lz=Lz,
                                V0_pin=V0_pin)
    dxdt = alpha * gyro * dw * (hz + pinning) + gyro * dw * pi2 * (
        -hy * math.cos(phi) + hx * math.sin(phi))
    dphidt = gyro * (hz + pinning) + alpha * gyro * pi2 * (hy * math.cos(phi) -
                                                           hx * math.sin(phi))
    return dxdt, dphidt


@njit
def get_current_contributions(phi, bj, Hshe, je, alpha, beta, dw):
    bjj = bj * je
    Hshej = Hshe * je * math.cos(phi)
    dxdt = (1 + alpha * beta) * bjj + alpha * gyro * dw * pi2 * Hshej
    dphidt = (beta - alpha) * (bjj / dw) + gyro * pi2 * Hshej
    return dxdt, dphidt


@njit
def compute_dynamics_(X, phi, hx, hy, hz, alpha, dw, bj, beta, Hshe, je, Hdmi,
                      Hk, Ms, V0_pin, pinning, Ly, Lz):
    ag = alpha * gyro
    dg = dw * gyro * (math.pi / 2.)
    field_contr = get_field_contribution(X,
                                         phi,
                                         hx=hx,
                                         hy=hy,
                                         hz=hz,
                                         alpha=alpha,
                                         dw=dw,
                                         Ms=Ms,
                                         V0_pin=V0_pin,
                                         pinning=pinning,
                                         Ly=Ly,
                                         Lz=Lz)
    j_contr = get_current_contributions(phi,
                                        bj=bj,
                                        Hshe=Hshe,
                                        je=je,
                                        alpha=alpha,
                                        beta=beta,
                                        dw=dw)
    dXdt = -gyro * dw * (Hk / 2.) * math.sin(
        2 * phi) + dg * Hdmi * math.sin(phi) + field_contr[0] + j_contr[0]

    dPhidt = ag * (Hk / 2.) * math.sin(2 * phi) - ag * pi2 * Hdmi * math.sin(
        phi) + field_contr[1] + j_contr[1]
    return dXdt, dPhidt


@njit
def compute_gamma_a(X, phi, Q, dw, hk, hx, hy, hdmi, bj, IECterm):
    pi2 = math.pi / 2.
    fact_gamma = -0.5 * hk * math.sin(
        2 * phi) + pi2 * hy * math.cos(phi) - pi2 * hx * math.sin(
            phi) + Q * pi2 * hdmi * math.sin(phi) + IECterm
    fact_stt = bj / dw
    return gyro * fact_gamma + fact_stt


@njit
def compute_gamma_b(X, phi, Q, dw, hshe, hz, hr, beta, bj, Ms, Lx, Ly, Lz,
                    V0_pin, V0_edge, pinning):
    pi2 = math.pi / 2
    hp = get_pinning_field(X,
                           Ms=Ms,
                           pinning=pinning,
                           Ly=Ly,
                           Lz=Lz,
                           V0_pin=V0_pin)
    he = get_edge_field(X, Lx, V0_edge)
    fact_gamma = Q * (he + hz + hp + pi2 * hshe *
                      math.cos(phi)) - beta * pi2 * hr * math.cos(phi)
    fact_stt = beta * bj / dw
    return gyro * fact_gamma + fact_stt


@njit
def compute_dynamics(X, phi, alpha, Q, dw, hx, hy, hz, hk, hdmi, hr, hshe,
                     beta, bj, Ms, Lx, Ly, Lz, V0_pin, V0_edge, pinning,
                     IECterm):
    gamma_a = compute_gamma_a(X, phi, Q, dw, hk, hx, hy, hdmi, bj, IECterm)
    gamma_b = compute_gamma_b(X, phi, Q, dw, hshe, hz, hr, beta, bj, Ms, Lx,
                              Ly, Lz, V0_pin, V0_edge, pinning)
    dXdt = dw * (gamma_a + alpha * gamma_b)
    dPhidt = -alpha * gamma_a + gamma_b
    return dXdt, dPhidt


@dataclass
class DomainWallDynamics:
    """Domain Wall dynamics class.
    :param alpha: Gilbert damping.
    :param domain_width: width of the domain.
    :param H: applied magnetic field vector.
    :param Hk: anisotropy field [A/m].
    :param Ms: magnetisation saturation [A/m].
    :param thickness: thickness of the FM material.
    :param SHE_angle: Spin Hall Effect angle.
    :param D: DMI constant.
    :param beta: STT beta parameter.
    :param p: STT polarisation efficiency.
    :param V0_pin: pinning voltage constant.
    :param pinning: the pinning period.
    :param Lx: z-dimension of the FM block.
    :param Ly: y-dimension of the FM block.
    :param Lz: z-dimension of the FM block.
    :param Q: up-down or down-up wall parameter (either 1 or -1).
    :param Hr: Rashba field [A/m].
    :param kappa: ratio of anisotropies, K = K/K0, (try 0-1)
    :param relativistic: whether to include relativistic effects?
    :param vmax_magnon: m/s, max magnon velocity => 2A/Sd
                        A -- exchange constant.
                        S -- |S1| + |S2| = total spin density of sublattices.
                        d -- interatomic distance.

    For classical formulation see:
    Current-driven dynamics of chiral ferromagnetic domain walls, Emori et al, 2013

    For relativisitc formulation see:
    Relativistic kinematics of a magnetic soliton, Caretta et al., 2020
    """
    H: VectorObj
    alpha: float
    Ms: float
    thickness: float
    SHE_angle: float
    D: float
    Ku: float  # The out-of-plane anisotropy constant
    Kp: float  # The in-plane anisotropy constant
    A: float = 1e-11  # J/m
    beta: float = 1
    p: float = 1
    V0_pin: float = 1.65e-20
    V0_edge: float = 0
    pinning: float = 30e-9
    Lx: float = 120e-9
    Ly: float = 120e-9
    Lz: float = 3e-9
    Q: int = 1
    Hr: float = 0
    kappa: float = 1
    moving_field: Literal["perpendicular", "inplane"] = "perpendicular"

    def __post__init__(self):
        # in post init we already have p
        self.bj = bohr_magneton * self.p / (echarge * self.Ms)
        self.je_driver = lambda t: 0
        denom = (2 * self.Ms * mu0 * echarge * self.thickness)
        self.Hshe = hbar * self.SHE_angle / denom
        self.hx, self.hy, self.hz = self.H.get_cartesian()
        self.dw0 = self.get_unrelaxed_domain_width()
        if self.moving_field == "perpendicular":
            self.Hk = self.get_perpendicular_anisotropy_field()
        elif self.moving_field == "inplane":
            self.Hk = self.get_inplane_anisotropy_field()

    def get_unrelaxed_domain_width(self):
        """Domain width is based off the effective perpendicular anisotropy.
        We reduce the perpendicular anisotropy by demagnetising field"""
        # Keff = self.Ku - 0.5*mu0*(self.Ms/mu0)**2
        Keff = self.Ku - (0.5 / mu0) * (self.Ms**2)
        return math.sqrt(self.A / Keff)

    def set_current_function(self, driver: Callable):
        """
        :param driver: A function of time that returns the current density
        """
        self.je_driver = driver

    def get_Hdmi(self, domain_width):
        """Returns the DMI field"""
        return self.D / (mu0 * self.Ms * domain_width)

    def get_perpendicular_anisotropy_field(self):
        """Returns the perpeanisotropy field"""
        return 2 * self.Ku / (mu0 * self.Ms)

    def get_inplane_anisotropy_field(self):
        """Returns the in-plane anisotropy field"""
        return 2 * self.Kp / (mu0 * self.Ms)


@dataclass
class MultilayerWallDynamics:
    layers: List[DomainWallDynamics]
    J: float = 0
    vector_size: int = 3  # 3 for X, phi, delta

    def __post_init__(self):
        if len(self.layers) > 2:
            raise ValueError(
                "BilayerWallDynamics only supports up to 2 layers")

    def multilayer_dw_llg(self, t, vec):
        """Solve the Thiaville llg equation for LLG.
        :param t: current simulation time.
        :param vec: contains [X, phi, delta], current DW position, its angle and domain width.
        :returns (dXdt, dPhidt, dDeltad): velocity and change of angle and domain width.
        """
        # vector is X1, phi1, X2, phi2, ...
        layer: DomainWallDynamics
        new_vec = []
        for i, layer in enumerate(self.layers):
            je_at_t = layer.je_driver(t=t)
            reduced_alpha = (1. + layer.alpha**2)
            lx = vec[self.vector_size * i]
            lphi = vec[(self.vector_size * i) + 1]
            ldomain_width = vec[(self.vector_size * i) + 2]
            if len(self.layers) == 1:
                Jterm = 0
            else:
                Jterm = 2 * self.J / (layer.Ms * mu0 * layer.thickness)
                otherphi = vec[self.vector_size * (i - 1) + 1]
                Jterm *= math.sin(lphi - otherphi)

            hdmi = layer.get_Hdmi(ldomain_width)
            dXdt, dPhidt = compute_dynamics(lx,
                                            phi=lphi,
                                            Q=layer.Q,
                                            hx=layer.hx,
                                            hy=layer.hy,
                                            hz=layer.hz,
                                            alpha=layer.alpha,
                                            dw=ldomain_width,
                                            bj=layer.bj * je_at_t,
                                            hr=layer.Hr,
                                            beta=layer.beta,
                                            hshe=layer.Hshe * je_at_t,
                                            hdmi=hdmi,
                                            hk=layer.Hk,
                                            Ms=layer.Ms,
                                            IECterm=Jterm,
                                            V0_pin=layer.V0_pin,
                                            V0_edge=layer.V0_edge,
                                            pinning=layer.pinning,
                                            Lx=layer.Lx,
                                            Ly=layer.Ly,
                                            Lz=layer.Lz)
            dXdt = dXdt / reduced_alpha
            dPhidt = dPhidt / reduced_alpha
            pref = gyro / (layer.alpha * mu0 * layer.Ms * layer.thickness)
            # domain width relaxation from Thiaville
            dDeltadt = pref * (layer.A / ldomain_width - ldomain_width *
                               (layer.Ku + layer.Kp * math.sin(lphi)**2))
            new_vec.extend([dXdt, dPhidt, dDeltadt])
        return new_vec

    def run(self,
            sim_time: float,
            starting_conditions: List[float],
            max_step: float = 1e-10):
        """Run simulation of DW dynamics
        :param sim_time: total simulation time (simulation units).
        :param starting_conditions: starting position and angle of the DW.
        :param max_step: maximum allowed step of the RK45 method.
        """
        integrator = RK45(fun=self.multilayer_dw_llg,
                          t0=0.,
                          first_step=1e-16,
                          max_step=max_step,
                          y0=starting_conditions,
                          rtol=1e-12,
                          t_bound=sim_time)
        result = defaultdict(list)
        while True:
            integrator.step()
            if integrator.status == 'failed':
                print("Failed to converge")
                break
            layer_vecs = integrator.y
            result['t'].append(integrator.t)
            for i, _ in enumerate(self.layers):
                x, phi, dw = layer_vecs[self.vector_size * i], layer_vecs[
                    self.vector_size * i +
                    1], layer_vecs[self.vector_size * i + 2]
                vel = (x - integrator.y_old[2 * i]) / integrator.step_size
                result[f'dw_{i}'].append(dw)
                result[f'v_{i}'].append(vel)
                result[f'x_{i}'].append(x)
                result[f'phi_{i}'].append(phi)
            if integrator.status == 'finished':
                break

        return result
