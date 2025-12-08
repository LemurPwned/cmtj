"""
This piece of code was inspired by: https://github.com/micromagnetics/70LinesOfNumpy
Reference: https://arxiv.org/abs/1411.7188.

Hence, this piece of code is a modified version of the original code and still follows LGPL license.
Note on dipole and demag tensor:

    - The dipole tensor excludes "itself" interaction.
    - The demag tensor includes "itself" interaction.

    Ad. dipole
    For macrospin mask the other object in the M field, compute the field and return value of the
    field in the masked region -- that's the magnitude of the dipole coupling field.
"""

from math import asinh, atan, log, pi, sqrt

import numpy as np
from numba import jit

EPS = 1e-18


def _aharoni_demag_factor_z(a: float, b: float, c: float) -> float:
    """
    Demagnetising factor for a rectangular prism magnetised along the z-axis.

    The implementation follows Aharoni, J. Appl. Phys. 83 (1998) 3432.
    All edge lengths must be positive and expressed in the same units.
    """
    a, b, c = map(float, (a, b, c))
    if min(a, b, c) <= 0.0:
        raise ValueError("Edge lengths must be positive.")

    r = sqrt(a * a + b * b + c * c)
    r_ab = sqrt(a * a + b * b)
    r_ac = sqrt(a * a + c * c)
    r_bc = sqrt(b * b + c * c)

    def _log_ratio(numerator: float, denominator: float) -> float:
        ratio = numerator / denominator
        if ratio <= 0.0:
            raise ValueError("Logarithm argument must be positive.")
        return log(ratio)

    pi_dz = 0.0
    pi_dz += (b * b - c * c) / (2.0 * b * c) * _log_ratio(r - a, r + a)
    pi_dz += (a * a - c * c) / (2.0 * a * c) * _log_ratio(r - b, r + b)
    pi_dz += (b / (2.0 * c)) * _log_ratio(r_ab + a, r_ab - a)
    pi_dz += (a / (2.0 * c)) * _log_ratio(r_ab + b, r_ab - b)
    pi_dz += (c / (2.0 * a)) * _log_ratio(r_bc - b, r_bc + b)
    pi_dz += (c / (2.0 * b)) * _log_ratio(r_ac - a, r_ac + a)
    pi_dz += 2.0 * atan(a * b / (c * r))
    pi_dz += (a**3 + b**3 - 2.0 * c**3) / (3.0 * a * b * c)
    pi_dz += (a * a + b * b - 2.0 * c * c) / (3.0 * a * b * c) * r
    pi_dz += (c / (a * b)) * (r_ac + r_bc)
    pi_dz -= (r_ab**3 + r_bc**3 + r_ac**3) / (3.0 * a * b * c)
    return pi_dz / pi


def rectangular_prism_demag_factors(width: float, length: float, height: float):
    """
    Compute (Dx, Dy, Dz) for a uniformly magnetised rectangular prism.

    Args:
        width: Prism extent along the x-axis.
        length: Prism extent along the y-axis.
        height: Prism extent along the z-axis.

    Returns:
        Tuple of demagnetising factors (Dx, Dy, Dz).
    """
    dz = _aharoni_demag_factor_z(width, length, height)
    dx = _aharoni_demag_factor_z(length, height, width)
    dy = _aharoni_demag_factor_z(height, width, length)
    return dx, dy, dz


@jit
def f(p):
    x, y, z = abs(p[0]), abs(p[1]), abs(p[2])
    return (
        +y / 2.0 * (z**2 - x**2) * asinh(y / (sqrt(x**2 + z**2) + EPS))
        + z / 2.0 * (y**2 - x**2) * asinh(z / (sqrt(x**2 + y**2) + EPS))
        - x * y * z * atan(y * z / (x * sqrt(x**2 + y**2 + z**2) + EPS))
        + 1.0 / 6.0 * (2 * x**2 - y**2 - z**2) * sqrt(x**2 + y**2 + z**2)
    )


# newell g
@jit
def g(p):
    x, y, z = p[0], p[1], abs(p[2])
    return (
        +x * y * z * asinh(z / (sqrt(x**2 + y**2) + EPS))
        + y / 6.0 * (3.0 * z**2 - y**2) * asinh(x / (sqrt(y**2 + z**2) + EPS))
        + x / 6.0 * (3.0 * z**2 - x**2) * asinh(y / (sqrt(x**2 + z**2) + EPS))
        - z**3 / 6.0 * atan(x * y / (z * sqrt(x**2 + y**2 + z**2) + EPS))
        - z * y**2 / 2.0 * atan(x * z / (y * sqrt(x**2 + y**2 + z**2) + EPS))
        - z * x**2 / 2.0 * atan(y * z / (x * sqrt(x**2 + y**2 + z**2) + EPS))
        - x * y * sqrt(x**2 + y**2 + z**2) / 3.0
    )


def set_n_demag(n_demag, n, dx, c, permute, func):
    it = np.nditer(n_demag[:, :, :, c], flags=["multi_index"], op_flags=["writeonly"])
    drrr = np.prod(dx)
    while not it.finished:
        value = 0.0
        for i in np.rollaxis(np.indices((2,) * 6), 0, 7).reshape(64, 6):
            idx = list(
                map(
                    lambda k: (it.multi_index[k] + n[k] - 1) % (2 * n[k] - 1) - n[k] + 1,
                    range(3),
                )
            )
            value += (-1) ** sum(i) * func(list(map(lambda j: (idx[j] + i[j] - i[j + 3]) * dx[j], permute)))
        it[0] = -value / (4 * pi * drrr)
        it.iternext()


def get_full_demag_tensor(n, dx):
    """
    Get the full demag tensor for a given cells n and their sizes dx.
    """
    n_demag = np.zeros([2 * i - 1 for i in n] + [6])
    for i, t in enumerate(
        (
            (f, 0, 1, 2),
            (g, 0, 1, 2),
            (g, 0, 2, 1),
            (f, 1, 2, 0),
            (g, 1, 2, 0),
            (f, 2, 0, 1),
        )
    ):
        set_n_demag(n_demag, n, dx, i, t[1:], t[0])
    n_demag = n_demag[: n[0], : n[1], : n[2], :]
    # N11, N12, N13
    # N21, N22, N23
    # N31, N32, N33
    # ===
    # 0  1  2
    # 1  3  4
    # 2  4  5
    tensor = np.zeros_like((*n, 3, 3), dtype=n_demag.dtype)
    tensor[..., 0, 0] = n_demag[..., 0]
    tensor[..., 0, 1] = n_demag[..., 1]
    tensor[..., 0, 2] = n_demag[..., 2]
    tensor[..., 1, 0] = n_demag[..., 1]
    tensor[..., 1, 1] = n_demag[..., 3]  # N22
    tensor[..., 1, 2] = n_demag[..., 4]  # N23
    tensor[..., 2, 0] = n_demag[..., 2]  # N31 = N13
    tensor[..., 2, 1] = n_demag[..., 4]  # N32 = N23
    tensor[..., 2, 2] = n_demag[..., 5]  # N33
    return tensor
