import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Literal

from numba import njit
from scipy.integrate import RK45

from ..utils import bohr_magneton, echarge, gyromagnetic_ratio, hbar, mu0
from ..utils.general import VectorObj

gyro = gyromagnetic_ratio
pi2 = math.pi / 2.0


class DW:
    """Initial conditions for the phi of DW equation."""

    NEEL_RIGHT = 0
    NEEL_LEFT = math.pi
    BLOCH_UP = math.pi / 2.0
    BLOCH_DOWN = 3.0 * math.pi / 2.0


class DWRelax:
    NO_RELAX = 0
    STATIC = 1
    DYNAMIC = 2


@njit
def get_pinning_field(X, Ms, pinning, Ly, Lz, V0_pin):
    arg = X * math.pi / pinning
    dVdx = 2 * math.pi * V0_pin * math.sin(arg) * math.cos(arg)
    denom = 2 * mu0 * Ms * Lz * Ly
    return -(1.0 / denom) * dVdx


@njit
def get_edge_field(X, Lx, V0_edge):
    c = Lx / 2
    arg = (X - (Lx / 2)) / c
    p = 6
    return -p * V0_edge * (math.sinh(arg) * math.cosh(arg) ** (p - 1) / c)


@njit
def get_field_contribution(X, phi, hx, hy, hz, alpha, dw, Ms, V0_pin, pinning, Ly, Lz):
    pinning = get_pinning_field(X, Ms=Ms, pinning=pinning, Ly=Ly, Lz=Lz, V0_pin=V0_pin)
    dxdt = alpha * gyro * dw * (hz + pinning) + gyro * dw * pi2 * (-hy * math.cos(phi) + hx * math.sin(phi))
    dphidt = gyro * (hz + pinning) + alpha * gyro * pi2 * (hy * math.cos(phi) - hx * math.sin(phi))
    return dxdt, dphidt


@njit
def compute_gamma_a(X, phi, Q, dw, hk, hx, hy, hdmi, bj, IECterm):
    pi2 = math.pi / 2.0
    fact_gamma = (
        -0.5 * hk * math.sin(2 * phi)
        - pi2 * hy * math.cos(phi)
        + pi2 * hx * math.sin(phi)
        + Q * pi2 * hdmi * math.sin(phi)
        + IECterm
    )
    fact_stt = bj / dw
    return gyro * fact_gamma + fact_stt


@njit
def compute_gamma_b(X, phi, Q, dw, hshe, hz, hr, beta, bj, Ms, Lx, Ly, Lz, V0_pin, V0_edge, pinning):
    pi2 = math.pi / 2
    hp = get_pinning_field(X, Ms=Ms, pinning=pinning, Ly=Ly, Lz=Lz, V0_pin=V0_pin)
    he = get_edge_field(X, Lx, V0_edge)
    fact_gamma = Q * (he + hz + hp + pi2 * hshe * math.cos(phi)) - beta * pi2 * hr * math.cos(phi)
    fact_stt = beta * bj / dw
    return gyro * fact_gamma + fact_stt


@njit
def compute_dynamics(
    X,
    phi,
    delta,
    alpha,
    Q,
    hx,
    hy,
    hz,
    hk,
    hdmi,
    hr,
    hshe,
    beta,
    bj,
    Ms,
    Lx,
    Ly,
    Lz,
    V0_pin,
    V0_edge,
    pinning,
    IECterm,
    thickness,
    A,
    Ku,
    Kp,
):
    gamma_a = compute_gamma_a(X, phi, Q, delta, hk, hx, hy, hdmi, bj, IECterm)
    gamma_b = compute_gamma_b(X, phi, Q, delta, hshe, hz, hr, beta, bj, Ms, Lx, Ly, Lz, V0_pin, V0_edge, pinning)
    dXdt = delta * (gamma_a + alpha * gamma_b)
    dPhidt = -alpha * gamma_a + gamma_b
    pref = gyro / (alpha * mu0 * Ms * thickness)
    # domain width relaxation from Thiaville
    dDeltadt = pref * (A / delta - delta * (Ku + Kp * math.sin(phi) ** 2))
    # dDeltadt  = 0
    return dXdt, dPhidt, dDeltadt


@dataclass
class DomainWallDynamics:
    """Domain Wall dynamics class.
    :param H: applied magnetic field vector.
    :param alpha: Gilbert damping.
    :param Ms: magnetisation saturation [A/m].
    :param thickness: thickness of the FM material.
    :param SHE_angle: Spin Hall Effect angle.
    :param D: DMI constant.
    :param Ku: perpendicular anisotropy constant.
    :param Kp: inplane anisotropy constant.
    :param A: exchange constant.
    :param beta: STT beta parameter.
    :param p: STT polarisation efficiency.
    :param V0_pin: pinning voltage constant.
    :param V0_edge: edge voltage constant.
    :param pinning: the pinning period.
    :param Lx: z-dimension of the FM block.
    :param Ly: y-dimension of the FM block.
    :param Lz: z-dimension of the FM block.
    :param Q: up-down or down-up wall parameter (either 1 or -1).
    :param Hr: Rashba field [A/m].
    :param moving_field: whether the anisotropy field is perpendicular or parallel
    :param relax_dw: whether to relax the domain width. See DWRelax class.
    For classical formulation see:
    Current-driven dynamics of chiral ferromagnetic domain walls, Emori et al, 2013
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
    moving_field: Literal["perpendicular", "inplane"] = "perpendicular"
    relax_dw: DWRelax = DWRelax.STATIC
    dw0: float = field(init=False)

    def __post_init__(self):
        # in post init we already have p
        self.bj = bohr_magneton * self.p / (echarge * self.Ms)
        self.je_driver = lambda t: 0
        denom = 2 * self.Ms * mu0 * echarge * self.thickness
        self.Hshe = hbar * self.SHE_angle / denom
        self.hx, self.hy, self.hz = self.H.get_cartesian()
        self.dw0 = self.get_unrelaxed_domain_width()

        if self.moving_field == "perpendicular":
            self.Hk = self.get_perpendicular_anisotropy_field()
        elif self.moving_field == "inplane":
            self.Hk = self.get_inplane_anisotropy_field()

    def get_unrelaxed_domain_width(self, effective=False):
        """Domain width is based off the effective perpendicular anisotropy.
        We reduce the perpendicular anisotropy by demagnetising field"""
        # Keff = self.Ku - 0.5*mu0*(self.Ms)**2
        Keff = self.Ku - (0.5 * mu0) * (self.Ms**2) if effective else self.Ku
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
    layers: list[DomainWallDynamics]
    J: float = 0
    vector_size: int = 3  # 3 for X, phi, delta

    def __post_init__(self):
        if len(self.layers) > 2:
            raise ValueError("Wall dynamics supported up to 2 layers")

    def multilayer_dw_llg(self, t, vec):
        """Solve the Thiaville llg equation for LLG.
        :param t: current simulation time.
        :param vec: contains [X, phi, delta], current DW position, its angle and domain width.
        :returns (dXdt, dPhidt, dDeltad): velocity and change of angle and domain width.
        """
        # vector is X1, phi1, Delta1, X2, phi2, Delta2...
        layer: DomainWallDynamics
        new_vec = []
        for i, layer in enumerate(self.layers):
            je_at_t = layer.je_driver(t=t)
            reduced_alpha = 1.0 + layer.alpha**2
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
            dXdt, dPhidt, dDeltadt = compute_dynamics(
                X=lx,
                phi=lphi,
                delta=ldomain_width,
                Q=layer.Q,
                hx=layer.hx,
                hy=layer.hy,
                hz=layer.hz,
                alpha=layer.alpha,
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
                Lz=layer.Lz,
                A=layer.A,
                Ku=layer.Ku,
                Kp=layer.Kp,
                thickness=layer.thickness,
            )
            dXdt = dXdt / reduced_alpha
            dPhidt = dPhidt / reduced_alpha
            if layer.relax_dw != DWRelax.DYNAMIC:
                dDeltadt = 0  # no relaxation in the ODE
            new_vec.extend([dXdt, dPhidt, dDeltadt])
        return new_vec

    def run(self, sim_time: float, starting_conditions: list[float], max_step: float = 1e-10):
        """Run simulation of DW dynamics.
        :param sim_time: total simulation time (simulation units).
        :param starting_conditions: starting position and angle of the DW.
        :param max_step: maximum allowed step of the RK45 method.
        """
        integrator = RK45(
            fun=self.multilayer_dw_llg,
            t0=0.0,
            first_step=1e-16,
            max_step=max_step,
            y0=starting_conditions,
            rtol=1e-12,
            t_bound=sim_time,
        )
        result = defaultdict(list)
        while True:
            integrator.step()
            if integrator.status == "failed":
                print("Failed to converge")
                break
            layer_vecs = integrator.y
            result["t"].append(integrator.t)
            for i, layer in enumerate(self.layers):
                x, phi, dw = (
                    layer_vecs[self.vector_size * i],
                    layer_vecs[self.vector_size * i + 1],
                    layer_vecs[self.vector_size * i + 2],
                )
                # static relaxation Thiaville
                if layer.relax_dw == DWRelax.STATIC:
                    ratio = layer.Kp / layer.Ku
                    dw = layer.dw0 / math.sqrt(1 + ratio * math.sin(phi) ** 2)
                vel = (x - integrator.y_old[2 * i]) / integrator.step_size
                result[f"dw_{i}"].append(dw)
                result[f"v_{i}"].append(vel)
                result[f"x_{i}"].append(x)
                result[f"phi_{i}"].append(phi)
                result[f"je_{i}"].append(layer.je_driver(t=integrator.t))
            if integrator.status == "finished":
                break

        return result
