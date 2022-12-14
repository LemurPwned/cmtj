import math
from collections import defaultdict
from dataclasses import dataclass

from scipy.integrate import RK45

from cmtj.models.sb import VectorObj
from cmtj.utils import bohr_magneton, echarge, gamma, hbar, mu0


@dataclass
class DomainWallDynamics:
    alpha: float
    domain_width: float
    H: VectorObj
    Ku: float
    Ms: float
    je: float
    thickness: float
    SHE_angle: float
    D: float
    beta: float = 1
    p: float = 1

    def __post_init__(self):
        self.dw = self.domain_width
        denom = (2 * self.Ms * mu0 * echarge * self.thickness)
        self.Hshe = hbar * self.SHE_angle * self.je / denom
        self.Hk = 2 * self.Ku / (mu0 * self.Ms)
        self.Hdmi = self.D / (mu0 * self.Ms * self.dw)

        self.bj = bohr_magneton * self.p * self.je / (echarge * self.thickness)

    def llg(self, t, vec):
        _, phi = vec
        reduced_alpha = (1 + self.alpha**2)
        hx, hy, hz = self.H.get_cartesian()
        ag = self.alpha * gamma
        dg = self.dw * gamma * math.pi / 2
        pi2 = math.pi / 2
        dXdt = ag * self.dw * hz - dg * (self.Hk / 2) * math.sin(2 * phi) + (
            1 + self.alpha * self.beta) * self.bj + ag * dg * math.cos(
                phi) * self.Hshe + dg * self.Hdmi * math.sin(
                    phi) - dg * hy * math.cos(phi) + dg * hx * math.sin(phi)

        dPhidt = gamma * hz + ag * (self.Hk / 2) * math.sin(
            2 * phi) + (self.beta - self.alpha) * (
                self.bj / self.dw) + gamma * pi2 * self.Hshe * math.cos(
                    phi) - ag * pi2 * self.Hdmi * math.sin(
                        phi) + ag * pi2 * hy * math.cos(
                            phi) - ag * pi2 * hx * math.sin(phi)
        return [dXdt / reduced_alpha, dPhidt / reduced_alpha]

    def run(self,
            sim_time: float,
            x0: float = 0,
            phi0: float = 0,
            max_steps: int = 1000):
        max_steps = int(max_steps)
        integrator = RK45(fun=self.llg,
                          t0=0,
                          first_step=1e-16,
                          y0=[x0, phi0],
                          rtol=1e-8,
                          t_bound=sim_time)
        result = defaultdict(list)
        for step in range(max_steps):
            integrator.step()
            result['t'].append(integrator.t)
            result['x'].append(integrator.y[0])
            result['phi'].append(integrator.y[1])

            if integrator.status == 'finished':
                break

        return result
