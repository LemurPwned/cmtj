"""
Fast Smit-Beljers Hysteresis Loops with relax_field_path

GOAL:
This example computes in-plane hysteresis loops of an antiferromagnetically coupled Co/Ru/Co-like
bilayer at several field angles with the static Smit-Beljers (SB) model, in two ways: the classic
point-by-point sweep (`Solver.adam_gradient_descent` at every field, started from the previous
equilibrium, as in docs/tutorials/SBModel.ipynb), and `Solver.relax_field_path`, which follows the
whole sweep for all field angles at once. It prints the run time of both and how far apart the
equilibria are.

CRITICAL SIMULATION MECHANISMS:
1. SB statics: equilibria are energy minima of the solver's full energy (Zeeman, uniaxial in-plane
   anisotropy, demagnetisation, bilinear + biquadratic interlayer exchange). SB models take Ms in A/m.
2. Hysteresis: each field step starts from the previous equilibrium, so metastable states persist until
   they lose stability (spin-flop / switching fields depend on the sweep direction).
3. relax_field_path: per-direction damped Newton with line search, compiled once per layout; all
   field angles are relaxed together, and new material parameters need no recompile.
4. Comparison: the loops (net moment along the field) agree; per layer the methods may settle in
   mirror-image scissor states about the field, within ~1e-2 mJ/m^2 in energy
   (`Solver.create_energy`), either one lower. Near zero field the loops can switch one step apart.

The figure shows the net moment along the field for both methods at each field angle.
"""

import contextlib
import time

import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

from cmtj.models import LayerSB, Solver, VectorObj
from cmtj.utils import mu0

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Co/Ru/Co-like bilayer (SB model: Ms in A/m, Kv in J/m^3, J1/J2 in J/m^2)
Ms = 1.65 / mu0  # 1.65 T
thicknesses = (4e-9, 3e-9)  # unequal Co layers break the layer-swap symmetry [m]
Ku = 2e3  # weak in-plane uniaxial anisotropy along x [J/m^3]
J1 = -1.2e-4  # antiferromagnetic bilinear IEC [J/m^2]
J2 = -5e-5  # biquadratic IEC [J/m^2]
demag = VectorObj.from_cartesian(0, 0, 1)  # thin film

layers = [LayerSB(_id=i, thickness=thicknesses[i], Kv=VectorObj(0, 0, Ku), Ks=0, Ms=Ms, Ndemag=demag) for i in range(2)]
solver = Solver(layers=layers, J1=[J1], J2=[J2])

# field sweep: +H -> -H -> +H, in-plane at several angles to the easy axis
Hmax, n_half = 80e3, 40  # [A/m]
down = np.linspace(Hmax, -Hmax, n_half)
Hpath = np.concatenate([down, down[::-1]])
angles_deg = np.array([0.0, 30.0, 60.0])
angles = np.deg2rad(angles_deg)


def m_along(coords, phi_h):
    """Net moment along the field, thickness-weighted over the layers."""
    t, p = coords[..., 0::2], coords[..., 1::2]
    return np.average(np.sin(t) * np.cos(p - phi_h), axis=-1, weights=thicknesses)


def cartesian(coords):
    t, p = coords[..., 0::2], coords[..., 1::2]
    return np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)], -1)


# (a) classic: one gradient descent per field point, per angle
t0 = time.perf_counter()
eq_adam = np.empty((len(Hpath), len(angles), 4))
for j, phi_h in enumerate(angles):
    # start along the first field, layer phi split by +-1e-5 rad (same start as relax_field_path)
    position = [np.pi / 2, phi_h + 1e-5, np.pi / 2, phi_h - 1e-5]
    for k, H in enumerate(Hpath):
        solver.set_H(VectorObj(np.pi / 2, phi_h, H))
        position = solver.adam_gradient_descent(position, max_steps=int(1e6), perturbation=0)
        eq_adam[k, j] = position
t_adam = time.perf_counter() - t0

# (b) relax_field_path: the whole sweep for all angles at once (first call includes compilation)
t0 = time.perf_counter()
eq_path = solver.relax_field_path(Hpath, angles)
t_path_first = time.perf_counter() - t0
t0 = time.perf_counter()
eq_path = solver.relax_field_path(Hpath, angles)
t_path = time.perf_counter() - t0

dm = np.linalg.norm(cartesian(eq_adam) - cartesian(eq_path), axis=-1)  # (steps, angles, layers)
print(f"Field points: {len(Hpath)} x {len(angles)} angles")
print(f"adam_gradient_descent per point: {t_adam:.2f} s")
print(f"relax_field_path: {t_path_first:.2f} s (first call, with compilation), {t_path:.3f} s (cached)")
dm_par = np.abs(m_along(eq_adam, angles[:, None]) - m_along(eq_path, angles[:, None]))
print(f"loop |delta m_par| between methods: max {dm_par.max():.2e}, median {np.median(dm_par):.2e}")

# per layer the two methods can pick mirror-image scissor states (about the field): compare their energies
Hsym = sym.symbols("H_x H_y H_z")
coords = [c for layer in solver.layers for c in layer.get_coord_sym()]
energy = sym.lambdify([*coords, *Hsym], solver.create_energy(H=sym.ImmutableMatrix(Hsym)))
Hxyz = np.moveaxis(Hpath[:, None, None] * np.stack([np.cos(angles), np.sin(angles), np.zeros_like(angles)], -1), -1, 0)
E_adam, E_path = energy(*np.moveaxis(eq_adam, -1, 0), *Hxyz), energy(*np.moveaxis(eq_path, -1, 0), *Hxyz)
print(f"per-layer |delta m|: max {dm.max():.2e}, points > 1e-3: {(dm.max(-1) > 1e-3).sum()} of {dm[..., 0].size}")
print(
    f"energy E_adam - E_path [mJ/m^2]: min {1e3 * (E_adam - E_path).min():.2e}, max {1e3 * (E_adam - E_path).max():.2e}"
)

with plt.style.context(["science", "no-latex"]):
    fig, axs = plt.subplots(1, len(angles), figsize=(7, 2.4), dpi=300, sharey=True)
    for j, ax in enumerate(axs):
        ax.plot(Hpath / 1e3, m_along(eq_adam[:, j], angles[j]), "o", ms=2.5, color="gray", label="adam per point")
        ax.plot(Hpath / 1e3, m_along(eq_path[:, j], angles[j]), "-", lw=1, color="crimson", label="relax_field_path")
        ax.set_title(rf"$\phi_H$ = {angles_deg[j]:.0f}$^\circ$")
        ax.set_xlabel("H (kA/m)")
    axs[0].set_ylabel(r"$m_\parallel$")
    axs[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig("./curated-examples/figures/sb-hysteresis-fast.png", dpi=300, bbox_inches="tight")
