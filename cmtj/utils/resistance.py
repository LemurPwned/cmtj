from functools import lru_cache
from typing import Callable, Union

import numpy as np
import sympy as sym

from .filters import Filters

EPS = np.finfo("float64").resolution


def compute_sd(dynamic_r: np.ndarray, dynamic_i: np.ndarray, integration_step: float) -> np.ndarray:
    """Computes the SD voltage.

    :param dynamic_r: magnetoresistance from log
    :param dynamic_i: excitation current
    :param integration_step: integration paramemter from run_simulation
    """
    SD = -dynamic_i * dynamic_r
    fs = 1.0 / integration_step
    SD_dc = Filters.butter_lowpass_filter(SD, cutoff=10e6, fs=fs, order=3)
    return np.mean(SD_dc)


def compute_resistance(
    Rx0: list[float],
    Ry0: list[float],
    AMR: list[float],
    AHE: list[float],
    SMR: list[float],
    m: Union[list[float], np.ndarray],
    l: list[float],
    w: list[float],
) -> tuple[list[float], list[float]]:
    """Computes the resistance of the system.

    If you want to compute the resistance for an entire time series, pass m as a 3D array
    with shape [number_of_layers, 3, T], where T is the time component.
    [number_of_layers, 3, T] where T is the time component.

    :param Rx0: resistance offset in longitudinal direction
    :param Ry0: resistance offset in transverse direction
    :param AMR: anisotropic magnetoresistance
    :param AHE: anomalous Hall effect
    :param SMR: spin Hall magnetoresistance
    :param m: magnetisation of the layers. Shape [number_of_layers, 3, T]
    :param l: length of the layers
    :param w: width of the layers
    """
    number_of_layers = len(Rx0)
    if not isinstance(m, np.ndarray):
        m = np.asarray(m)
    if m.ndim == 2:
        SxAll = np.zeros((number_of_layers,))
        SyAll = np.zeros((number_of_layers,))

    elif m.ndim == 3:
        SxAll = np.zeros((number_of_layers, m.shape[2]))
        SyAll = np.zeros((number_of_layers, m.shape[2]))

    for i in range(number_of_layers):
        w_l = w[i] / l[i]
        SxAll[i] = Rx0[i] + (AMR[i] * m[i, 0] ** 2 + SMR[i] * m[i, 1] ** 2)
        SyAll[i] = Ry0[i] + 0.5 * AHE[i] * m[i, 2] + (w_l) * (SMR[i] + AMR[i]) * m[i, 0] * m[i, 1]
    return SxAll, SyAll


def compute_gmr(Rp: float, Rap: float, m1: np.ndarray, m2: np.ndarray) -> np.ndarray:
    """Computes the GMR using parallel and antiparallel resistance.
    :param Rp: parallel resistance
    :param Rap: antiparallel resistance
    :param m1: magnetisation of layer 1
    :param m2: magnetisation of layer 2"""
    return Rp + 0.5 * (Rap - Rp) * (1 - np.sum(m1 * m2, axis=0))


def calculate_magnetoresistance(Rp: float, Rap: float, m: np.ndarray) -> np.ndarray:
    """Computes the magnetoresistance using parallel and antiparallel resistance.
    :param Rp: parallel resistance
    :param Rap: antiparallel resistance
    :param m: magnetisation, 2 layers of shape [2, 3, T] where T is the time component
    """
    if not isinstance(m, np.ndarray):
        m = np.asarray(m)
    if m.shape[0] != 2:
        raise ValueError("The magnetoresistance can only be computed for 2 layers" f". Current shape {m.shape}")
    return Rp + 0.5 * (Rap - Rp) * np.sum(m[0] * m[1], axis=0)


def calculate_resistance_parallel(
    Rx0: list[float],
    Ry0: list[float],
    AMR: list[float],
    AHE: list[float],
    SMR: list[float],
    m: list[float],
    l: list[float],
    w: list[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Calculates the resistance of the system in parallel.
    If you want to compute the resistance for an entire time series, pass m as a 3D array.
    [number_of_layers, 3, T] where T is the time component.
    Uses Kim's formula from the paper:
    https://link.aps.org/doi/10.1103/PhysRevLett.116.097201

    :param Rx0: resistance offset in longitudinal direction
    :param Ry0: resistance offset in transverse direction
    :param AMR: anisotropic magnetoresistance
    :param AHE: anomalous Hall effect
    :param SMR: spin Hall magnetoresistance
    :param m: magnetisation of the layers. Shape [number_of_layers, 3, T]
    :param l: length of the layers
    :param w: width of the layers
    """
    SxAll, SyAll = compute_resistance(Rx0, Ry0, AMR, AHE, SMR, m, l, w)
    Rx = 1.0 / np.sum(1.0 / SxAll, axis=0)
    Ry = 1.0 / np.sum(1.0 / SyAll, axis=0)
    return Rx, Ry


def calculate_resistance_series(
    Rx0: list[float],
    Ry0: list[float],
    AMR: list[float],
    AHE: list[float],
    SMR: list[float],
    m: list[float],
    l: list[float],
    w: list[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Calculates the resistance of the system in series.
    If you want to compute the resistance for an entire time series, pass m as a 3D array.
    [number_of_layers, 3, T] where T is the time component.
    Uses Kim's formula from the paper:
    https://link.aps.org/doi/10.1103/PhysRevLett.116.097201

    :param Rx0: resistance offset in longitudinal direction
    :param Ry0: resistance offset in transverse direction
    :param AMR: anisotropic magnetoresistance
    :param AHE: anomalous Hall effect
    :param SMR: spin Hall magnetoresistance
    :param m: magnetisation of the layers. Shape [number_of_layers, 3, T]
    :param l: length of the layers
    :param w: width of the layers
    """
    SxAll, SyAll = compute_resistance(Rx0, Ry0, AMR, AHE, SMR, m, l, w)
    Rx = np.sum(SxAll, axis=0)
    Ry = np.sum(SyAll, axis=0)
    return Rx, Ry


def angular_calculate_resistance_gmr(
    Rp: float,
    Rap: float,
    theta_1: np.ndarray,
    phi_1: np.ndarray,
    theta_2: np.ndarray,
    phi_2: np.ndarray,
) -> np.ndarray:
    """Computes the GMR using parallel and antiparallel resistance.
    :param Rp: parallel resistance
    :param Rap: antiparallel resistance
    :param theta_1: angle of layer 1
    :param phi_1: angle of layer 1
    :param theta_2: angle of layer 2
    :param phi_2: angle of layer 2
    """
    m1 = np.array(
        [
            np.cos(theta_1) * np.cos(phi_1),
            np.cos(theta_1) * np.sin(phi_1),
            np.sin(theta_1),
        ]
    )
    m2 = np.array(
        [
            np.cos(theta_2) * np.cos(phi_2),
            np.cos(theta_2) * np.sin(phi_2),
            np.sin(theta_2),
        ]
    )
    return compute_gmr(Rp, Rap, m1, m2)


@lru_cache(maxsize=5)
def Rxx_symbolic(id: int, AMR: float, SMR: float) -> tuple[float, sym.Symbol, sym.Symbol, sym.Matrix]:
    """Compute the Rxx resistance for a given layer.
    :param id: layer id
    :param AMR: anisotropic magnetoresistance
    :param SMR: spin Hall magnetoresistance
    """
    theta1 = sym.Symbol(r"\theta_" + str(id))
    phi1 = sym.Symbol(r"\phi_" + str(id))
    m = sym.Matrix(
        [
            sym.sin(theta1) * sym.cos(phi1),
            sym.sin(theta1) * sym.sin(phi1),
            sym.cos(theta1),
        ]
    )
    return AMR * m[0] ** 2 + SMR * m[1] ** 2, theta1, phi1, m


@lru_cache(maxsize=5)
def Rxy_symbolic(
    id: int, AMR: float, SMR: float, AHE: float, w_l: float
) -> tuple[float, sym.Symbol, sym.Symbol, sym.Matrix]:
    """Compute the Rxy resistance for a given layer.
    :param id: layer id
    :param AMR: anisotropic magnetoresistance
    :param SMR: spin Hall magnetoresistance
    :param AHE: anomalous Hall effect
    :param w_l: width to length ratio
    """
    theta1 = sym.Symbol(r"\theta_" + str(id))
    phi1 = sym.Symbol(r"\phi_" + str(id))
    m = sym.Matrix(
        [
            sym.sin(theta1) * sym.cos(phi1),
            sym.sin(theta1) * sym.sin(phi1),
            sym.cos(theta1),
        ]
    )
    return (0.5 * AHE * m[-1]) + w_l * (SMR - AMR) * m[0] * m[1], theta1, phi1, m


def calculate_linearised_resistance(
    GMR: float,
    AMR: list[float],
    SMR: list[float],
) -> tuple[float, float, float, sym.Symbol, sym.Symbol, sym.Symbol, sym.Symbol]:
    """
    Compute the resistance of the two FM bilayer system from the linearised angles.

    :param GMR: GMR
    :param AMR: AMR
    :param SMR: SMR
    """

    Rxx1, theta1, phi1, m1 = Rxx_symbolic(1, AMR[0], SMR[0])
    Rxx2, theta2, phi2, m2 = Rxx_symbolic(2, AMR[1], SMR[1])
    GMR_resistance = GMR * (1 - (m1.dot(m2))) / 2.0
    return Rxx1, Rxx2, GMR_resistance, theta1, phi1, theta2, phi2


def Rxx_parallel_bilayer_expr() -> (
    tuple[
        Callable[
            [
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
            ],
            float,
        ],
        Callable[
            [
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
            ],
            float,
        ],
    ]
):
    """Get the symbolic expressions for the parallel and linearised resistance of a bilayer system.

    Signals:
    - GMR: GMR
    - AMR1: AMR of layer 1
    - SMR1: SMR of layer 1
    - AMR2: AMR of layer 2
    - SMR2: SMR of layer 2
    - stationary angles: [t1, p1, t2, p2]
    - linearised angles: [dt1, dp1, dt2, dp2]

    Function signatures
    - Rlin_func: linearised resistance function
        f(GMR, AMR1, SMR1, AMR2, SMR2, [t1, p1, t2, p2], [dt1, dp1, dt2, dp2])
    - R_func: series resistance function
        f(GMR, AMR1, SMR1, AMR2, SMR2, [t1, p1, t2, p2])

    :returns: linearised and parallel resistance functions
    """
    AMR_1 = sym.Symbol(r"\mathrm{AMR}_1")
    SMR_1 = sym.Symbol(r"\mathrm{SMR}_1")
    AMR_2 = sym.Symbol(r"\mathrm{AMR}_2")
    SMR_2 = sym.Symbol(r"\mathrm{SMR}_2")
    GMR_s = sym.Symbol(r"\mathrm{GMR}")
    R_1, t1, p1, m1 = Rxx_symbolic(1, AMR_1, SMR_1)
    R_2, t2, p2, m2 = Rxx_symbolic(2, AMR_2, SMR_2)
    gmr_term = GMR_s * (1 - m1.dot(m2)) / 2

    Rparallel = gmr_term + (R_1 * R_2) / (R_1 + R_2 + EPS)
    linearised_terms = sym.symbols(r"\partial\theta_1, \partial\phi_1, \partial\theta_2, \partial\phi_2")
    dRdtheta1 = sym.diff(Rparallel, t1) * linearised_terms[0]
    dRdphi1 = sym.diff(Rparallel, p1) * linearised_terms[1]
    dRdtheta2 = sym.diff(Rparallel, t2) * linearised_terms[2]
    dRdphi2 = sym.diff(Rparallel, p2) * linearised_terms[3]

    linearised_R = dRdtheta1 + dRdtheta2 + dRdphi1 + dRdphi2

    Rlin_func = sym.lambdify(
        [GMR_s, AMR_1, SMR_1, AMR_2, SMR_2, [t1, p1, t2, p2], linearised_terms],
        linearised_R,
    )
    R_func = sym.lambdify([GMR_s, AMR_1, SMR_1, AMR_2, SMR_2, [t1, p1, t2, p2]], Rparallel)

    return Rlin_func, R_func


def GMR_expr() -> Callable[[float, float, float, float], float]:
    """Get the symbolic expression for the GMR.

    :returns: GMR function: Callable[[theta1, phi1, theta2, phi2], float]
    """
    GMR_s = sym.Symbol(r"\mathrm{GMR}")
    theta1 = sym.Symbol(r"\theta_1")
    phi1 = sym.Symbol(r"\phi_1")
    m1 = sym.Matrix(
        [
            sym.sin(theta1) * sym.cos(phi1),
            sym.sin(theta1) * sym.sin(phi1),
            sym.cos(theta1),
        ]
    )
    theta2 = sym.Symbol(r"\theta_2")
    phi2 = sym.Symbol(r"\phi_2")
    m2 = sym.Matrix(
        [
            sym.sin(theta2) * sym.cos(phi2),
            sym.sin(theta2) * sym.sin(phi2),
            sym.cos(theta2),
        ]
    )
    Rf = GMR_s * (1 - m1.dot(m2)) / 2
    dRdt1 = sym.diff(Rf, theta1)
    dRdp1 = sym.diff(Rf, phi1)
    dRdt2 = sym.diff(Rf, theta2)
    dRdp2 = sym.diff(Rf, phi2)
    linearised_terms = sym.symbols(r"\partial\theta_1, \partial\phi_1, \partial\theta_2, \partial\phi_2")
    dRf = (
        dRdt1 * linearised_terms[0]
        + dRdp1 * linearised_terms[1]
        + dRdt2 * linearised_terms[2]
        + dRdp2 * linearised_terms[3]
    )

    Rf_func = sym.lambdify([GMR_s, [theta1, phi1, theta2, phi2]], Rf)
    dRf_func = sym.lambdify([GMR_s, [theta1, phi1, theta2, phi2], linearised_terms], dRf)
    return Rf_func, dRf_func


def Rxx_series_bilayer_expr() -> (
    tuple[
        Callable[
            [
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
            ],
            float,
        ],
        Callable[
            [
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
                float,
            ],
            float,
        ],
    ]
):
    """Get the symbolic expressions for the series and linearised resistance of a bilayer system.

    Signals:
    - GMR: GMR
    - AMR1: AMR of layer 1
    - SMR1: SMR of layer 1
    - AMR2: AMR of layer 2
    - SMR2: SMR of layer 2
    - stationary angles: [t1, p1, t2, p2]
    - linearised angles: [dt1, dp1, dt2, dp2]

    Function signatures
    - Rlin_func: linearised resistance function
        f(GMR, AMR1, SMR1, AMR2, SMR2, [t1, p1, t2, p2], [dt1, dp1, dt2, dp2])
    - R_func: series resistance function
        f(GMR, AMR1, SMR1, AMR2, SMR2, [t1, p1, t2, p2])

    :returns: linearised and series resistance functions
    """
    AMR_1 = sym.Symbol(r"\mathrm{AMR}_1")
    SMR_1 = sym.Symbol(r"\mathrm{SMR}_1")
    AMR_2 = sym.Symbol(r"\mathrm{AMR}_2")
    SMR_2 = sym.Symbol(r"\mathrm{SMR}_2")
    GMR_s = sym.Symbol(r"\mathrm{GMR}")
    R_1, t1, p1, m1 = Rxx_symbolic(1, AMR_1, SMR_1)
    R_2, t2, p2, m2 = Rxx_symbolic(2, AMR_2, SMR_2)
    gmr_term = GMR_s * (1 - m1.dot(m2)) / 2

    Rseries = gmr_term + R_1 + R_2
    linearised_terms = sym.symbols(r"\partial\theta_1, \partial\phi_1, \partial\theta_2, \partial\phi_2")
    dRdtheta1 = sym.diff(Rseries, t1) * linearised_terms[0]
    dRdphi1 = sym.diff(Rseries, p1) * linearised_terms[1]
    dRdtheta2 = sym.diff(Rseries, t2) * linearised_terms[2]
    dRdphi2 = sym.diff(Rseries, p2) * linearised_terms[3]

    linearised_R = dRdtheta1 + dRdtheta2 + dRdphi1 + dRdphi2

    Rlin_func = sym.lambdify(
        [GMR_s, AMR_1, SMR_1, AMR_2, SMR_2, [t1, p1, t2, p2], linearised_terms],
        linearised_R,
    )
    R_func = sym.lambdify([GMR_s, AMR_1, SMR_1, AMR_2, SMR_2, [t1, p1, t2, p2]], Rseries)

    return Rlin_func, R_func


def calculate_linearised_resistance_parallel(
    GMR: float,
    AMR: list[float],
    SMR: list[float],
    stationary_angles: list[float],
    linearised_angles: list[float],
):
    """
    Compute the parallel resistance of the two FM bilayer system from the linearised angles.

    :param GMR: GMR
    :param AMR: AMR
    :param SMR: SMR
    :param stationary_angles: stationary angles [t1, p1, t2, p2]
    :param linearised_angles: linearised angles [dt1, dp1, dt2, dp2]
    """
    t01, p01 = stationary_angles[:2]
    t02, p02 = stationary_angles[2:]
    dt1, dp1 = linearised_angles[:2]
    dt2, dp2 = linearised_angles[2:]
    Rxx1, Rxx2, GMR_resistance, theta1, phi1, theta2, phi2 = calculate_linearised_resistance(GMR, AMR, SMR)
    Rparallel = GMR_resistance
    if any(AMR) or any(SMR):
        Rparallel += (Rxx1 * Rxx2) / (Rxx1 + Rxx2 + EPS)
    elif GMR == 0:
        return 0, 0
    dRparallel = (
        sym.diff(Rparallel, theta1) * dt1
        + sym.diff(Rparallel, phi1) * dp1
        + sym.diff(Rparallel, theta2) * dt2
        + sym.diff(Rparallel, phi2) * dp2
    )
    if isinstance(dRparallel, (list, np.ndarray)):
        dRparallel = dRparallel[0]
    dRparallel = dRparallel.subs(
        {
            theta1: t01,
            phi1: p01,
            theta2: t02,
            phi2: p02,
        }
    ).evalf()
    Rparallel = Rparallel.subs(
        {
            theta1: t01,
            phi1: p01,
            theta2: t02,
            phi2: p02,
        }
    ).evalf()
    return dRparallel, Rparallel


def calculate_linearised_resistance_series(
    GMR: float,
    AMR: list[float],
    SMR: list[float],
    stationary_angles: list[float],
    linearised_angles: list[float],
) -> tuple[float, float]:
    """
    Compute the resistance of the two FM bilayer system from the linearised angles.
    :param GMR: GMR
    :param AMR: AMR
    :param SMR: SMR
    :param stationary_angles: stationary angles [t1, p1, t2, p2]
    :param linearised_angles: linearised angles [dt1, dp1, dt2, dp2]
    """
    t01, p01 = stationary_angles[:2]
    t02, p02 = stationary_angles[2:]
    dt1, dp1 = linearised_angles[:2]
    dt2, dp2 = linearised_angles[2:]
    Rxx1, Rxx2, GMR_resistance, theta1, phi1, theta2, phi2 = calculate_linearised_resistance(GMR, AMR, SMR)
    Rseries = GMR_resistance + Rxx1 + Rxx2
    dRseries = (
        sym.diff(Rseries, theta1) * dt1
        + sym.diff(Rseries, phi1) * dp1
        + sym.diff(Rseries, theta2) * dt2
        + sym.diff(Rseries, phi2) * dp2
    )

    dRseries = dRseries.subs(
        {
            theta1: t01,
            phi1: p01,
            theta2: t02,
            phi2: p02,
        }
    ).evalf()
    Rseries = Rseries.subs(
        {
            theta1: t01,
            phi1: p01,
            theta2: t02,
            phi2: p02,
        }
    ).evalf()
    return dRseries, Rseries
