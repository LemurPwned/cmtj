import contextlib
from typing import Union

import numpy as np
import sympy as sym
from numba import njit

from ..utils import VectorObj
from ..utils.solvers import RootFinder

EPS = np.finfo("float64").resolution
OMEGA = sym.Symbol(r"\omega", complex=True)


def real_deocrator(fn):
    """Using numpy real cast is way faster than sympy."""

    def wrap_fn(*args):
        return np.real(fn(*args))

    return wrap_fn


def _lu_decomposition_det(matrix: sym.Matrix):
    _, U, perm = matrix.LUdecomposition()
    # sgn = Permutation(perm).signature()
    sgn = 1
    return sgn * U.det()


def _berkowitz_det(matrix: sym.Matrix):
    return matrix.det(method="berkowitz")


def _root_solver_numerical(
    root_expr,
    n_layers: int = None,
    *args,
    normalise_roots_by_2pi: bool = False,
    ftol: float = 0.01e9,
    max_freq: float = 80e9,
):
    xtol = 1e-4
    if n_layers is None or n_layers <= 3:
        # makes it faster for small systems, otherise jit cost too high
        y = real_deocrator(njit(sym.lambdify(OMEGA, root_expr, "math")))
    else:
        y = real_deocrator(sym.lambdify(OMEGA, root_expr, "math"))
    mfreq = max_freq * 2 * np.pi if normalise_roots_by_2pi else max_freq
    r = RootFinder(xtol, mfreq, step=ftol, xtol=xtol, root_dtype="float16")
    roots = r.find(y)
    # convert to GHz
    # reduce unique solutions to 2 decimal places
    roots = np.unique(np.around(roots / 1e9, 2))
    return roots / (2 * np.pi) if normalise_roots_by_2pi else roots


# we add n_layers and args to make the signature compatible with numerical solver
def _root_solver_analytical(
    root_expr,
    n_layers: int = None,
    normalise_roots_by_2pi: bool = False,
    solve_direct: bool = False,
    **kwargs,
):
    # Try multiple approaches to find roots
    all_solutions = []

    # # Approach 1: Direct solving
    if solve_direct:
        with contextlib.suppress(Exception):
            direct_solutions = sym.solve(root_expr, OMEGA)
            all_solutions.extend(direct_solutions)

    # Approach 2: Factorized solving
    with contextlib.suppress(Exception):
        factorised = sym.factor(root_expr)
        factored_solutions = sym.solve(factorised, OMEGA)
        all_solutions.extend(factored_solutions)

    # Remove duplicates
    all_solutions = list(set(all_solutions))

    # More robust real check - evaluate numerically and check if imaginary part is negligible
    real_solutions = []
    for sol in all_solutions:
        try:
            numeric_val = complex(sol.evalf())
            if abs(numeric_val.imag) < 1e-10:  # More lenient than is_real
                real_solutions.append(numeric_val.real)
        except Exception:
            continue

    # Convert to numpy and filter
    if not real_solutions:
        return np.array([])

    all_roots = np.asarray(real_solutions) / 1e9

    # Filter: positive and within frequency range
    all_roots = all_roots[(all_roots > 0)]

    # Remove duplicates with tolerance
    if len(all_roots) > 1:
        all_roots = np.unique(np.around(all_roots, 6))  # Higher precision than numerical

    return all_roots / (2 * np.pi) if normalise_roots_by_2pi else all_roots


def _default_matrix_conversion(
    vector_or_matrix: Union[VectorObj, list[float], sym.Matrix],
) -> sym.ImmutableMatrix:
    if isinstance(vector_or_matrix, VectorObj):
        vector_or_matrix = vector_or_matrix.get_cartesian()
    elif isinstance(vector_or_matrix, list):
        vector_or_matrix = sym.ImmutableMatrix(vector_or_matrix)
    return sym.ImmutableMatrix(vector_or_matrix)
