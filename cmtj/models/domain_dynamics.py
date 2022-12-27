import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable

from numba import njit
from scipy.integrate import RK45

from ..utils import bohr_magneton, echarge, gyromagnetic_ratio, hbar, mu0
from .sb import VectorObj

gyro = gyromagnetic_ratio
pi2 = math.pi / 2.
"""
Change DMI sign for UP/DOWN domain wall (same chirality)
"""


class DW:
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
def compute_dynamics(X, phi, hx, hy, hz, alpha, dw, bj, beta, Hshe, je, Hdmi,
                     Hk, Ms, V0_pin, pinning, Ly, Lz):
    reduced_alpha = (1. + alpha**2)
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

    dXdt = dXdt / reduced_alpha
    dPhidt = dPhidt / reduced_alpha

    return dXdt, dPhidt


@dataclass
class DomainWallDynamics:
    alpha: float
    domain_width: float
    H: VectorObj
    Hk: float
    Ms: float
    thickness: float
    SHE_angle: float
    D: float
    beta: float = 1
    p: float = 1
    V0_pin: float = 1.65e-20
    pinning: float = 30e-9
    Lz: float = 3e-9
    Ly: float = 120e-9

    def __post_init__(self):
        self.dw = self.domain_width
        denom = (2 * self.Ms * mu0 * echarge * self.thickness)
        self.Hshe = hbar * self.SHE_angle / denom
        self.Hdmi = self.D / (mu0 * self.Ms * self.dw)

        self.bj = bohr_magneton * self.p / (echarge * self.Ms)
        self.hx, self.hy, self.hz = self.H.get_cartesian()
        self.je_driver = lambda t: 0

    def set_current_function(self, driver: Callable):
        """
        :param driver: A function of time that returns the current density
        """
        self.je_driver = driver

    def llg(self, t, vec):
        X, phi = vec
        je_at_t = self.je_driver(t=t)
        return compute_dynamics(X,
                                phi,
                                hx=self.hx,
                                hy=self.hy,
                                hz=self.hz,
                                alpha=self.alpha,
                                dw=self.dw,
                                bj=self.bj,
                                beta=self.beta,
                                Hshe=self.Hshe,
                                je=je_at_t,
                                Hdmi=self.Hdmi,
                                Hk=self.Hk,
                                Ms=self.Ms,
                                V0_pin=self.V0_pin,
                                pinning=self.pinning,
                                Lz=self.Lz,
                                Ly=self.Ly)

    def run(self,
            sim_time: float,
            x0: float = 0,
            phi0: float = 0,
            max_step=1e-10):
        integrator = RK45(fun=self.llg,
                          t0=0.,
                          first_step=1e-16,
                          max_step=max_step,
                          y0=[x0, phi0],
                          rtol=1e-12,
                          t_bound=sim_time)
        result = defaultdict(list)
        while True:
            integrator.step()
            x, phi = integrator.y
            vel = (x - integrator.y_old[0]) / integrator.step_size
            result['t'].append(integrator.t)
            result['v'].append(vel)
            result['x'].append(x)
            result['phi'].append(phi)
            if integrator.status == 'failed':
                print("Failed to converge")
                break
            if integrator.status == 'finished':
                break

        return result
