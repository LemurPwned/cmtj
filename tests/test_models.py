import numpy as np
from cmtj.models import Solver, LayerDynamic, VectorObj, LayerSB
from typing import Tuple


def test_sb_dynamic(two_layer_symbolic_dyn: Tuple[LayerDynamic]):
    J1 = -1e-3
    J2 = 1e-4

    for Hmag in np.linspace(-300e3, 300e3, 8):
        H = VectorObj(np.deg2rad(88), np.deg2rad(0.1), Hmag)
        solver_dyn = Solver(layers=two_layer_symbolic_dyn, J1=[J1], J2=[J2], H=H)
        pos = np.asarray(
            [np.deg2rad(88), np.deg2rad(0.1), np.deg2rad(88), np.deg2rad(0.1)]
        )
        # set perturbation to 0 to avoid numerical errors
        eq_sb, f_sb = solver_dyn.solve(
            init_position=pos, perturbation=0, max_steps=1e6, force_sb=True
        )
        eq_dyn, f_dyn, _ = solver_dyn.solve(
            init_position=pos, max_steps=1e6, perturbation=0
        )
        f_sb.sort()
        f_dyn.sort()
        assert np.allclose(eq_sb, eq_dyn)
        assert np.allclose(f_sb, f_dyn, atol=0.2)
        pos = eq_dyn


def test_sb_classic_dipole(two_layer_symbolic_classic: Tuple[LayerSB]):
    J1 = -1e-3
    J2 = 1e-4
    dipole = [
        VectorObj.from_cartesian(0, 1e-4, 0),
        VectorObj.from_cartesian(0, -1e-4, 0),
        VectorObj.from_cartesian(0, 0, 1e-6),
    ]
    for Hmag in np.linspace(-300e3, 300e3, 8):
        H = VectorObj(np.deg2rad(88), np.deg2rad(0.1), Hmag)
        solver_dyn = Solver(
            layers=two_layer_symbolic_classic, J1=[J1], J2=[J2], H=H, Ndipole=[dipole]
        )
        pos = np.asarray(
            [np.deg2rad(88), np.deg2rad(0.1), np.deg2rad(88), np.deg2rad(0.1)]
        )
        # set perturbation to 0 to avoid numerical errors
        eq_sb, f_sb = solver_dyn.solve(
            init_position=pos, perturbation=0, max_steps=1e6, force_sb=True
        )
        f_sb.sort()
        pos = eq_sb
