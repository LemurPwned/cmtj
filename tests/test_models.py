import numpy as np
import pytest
from cmtj.models import Solver, LayerDynamic, VectorObj, LayerSB, spin_diode_mr_channels
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


def _diode_solver(hoe_direction=None, hoe_dl_ratio=0):
    layers = [
        LayerDynamic(
            i,
            thickness=4e-9,
            Kv=VectorObj(np.deg2rad(90), np.deg2rad(20 + 40 * i), 5e3),
            Ks=0,
            Ms=1.6 / (4e-7 * np.pi),
            alpha=0.02,
            Ndemag=VectorObj.from_cartesian(0, 0, 1),
            hoe_direction=hoe_direction,
            hoe_dl_ratio=hoe_dl_ratio,
        )
        for i in range(2)
    ]
    H = VectorObj.from_cartesian(40e3, 0, 0)
    return Solver(layers=layers, J1=[-2e-4], J2=[0.0], H=H), H


def _diode(solver, H, pos):
    return solver.linearised_N_spin_diode(
        H=H,
        frequency=8e9,
        Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol(),
        Vdc_ex_value=1.0,
        zero_pos=pos,
        return_complex=True,
    ).ravel()


def test_spin_diode_lab_frame_drive():
    """A fixed rf field direction reproduces the default (phi-hat) drive when it is
    along phi-hat of both layers, and scales with its projection on phi-hat."""
    pos = np.array([np.pi / 2, 0.0, np.pi / 2, 0.0])  # both layers along x
    solver, H = _diode_solver()
    ref = _diode(solver, H, pos)
    solver, H = _diode_solver(VectorObj.from_cartesian(0, 1, 0))
    assert np.allclose(_diode(solver, H, pos), ref, rtol=1e-6)
    solver, H = _diode_solver(VectorObj.from_cartesian(1, 1, 0))
    assert np.allclose(_diode(solver, H, pos), ref / np.sqrt(2), rtol=1e-6)
    # complex output agrees with the real one
    real = solver.linearised_N_spin_diode(
        H=H,
        frequency=8e9,
        Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol(),
        Vdc_ex_value=1.0,
        zero_pos=pos,
    ).ravel()
    assert np.allclose(real, _diode(solver, H, pos).real)


def test_spin_diode_damping_like_ratio():
    """hoe_dl_ratio adds the drive of the effective field r (m x d): for m = x, d = y that is
    r z, i.e. -r along theta-hat; the response is linear in the (complex) ratio."""
    pos = np.array([np.pi / 2, 0.0, np.pi / 2, 0.0])
    d = VectorObj.from_cartesian(0, 1, 0)
    base = _diode(*_diode_solver(d), pos)
    r1 = _diode(*_diode_solver(d, 1.0), pos) - base
    r2 = _diode(*_diode_solver(d, 0.5 + 2j), pos) - base
    assert np.abs(r1).max() > 1e-3 * np.abs(base).max()
    assert np.allclose(r2, (0.5 + 2j) * r1, rtol=1e-6, atol=1e-12 * np.abs(r1).max())


def test_spin_diode_linear_in_excitation():
    """The response is linear in Vdc_ex_value also away from the plane (the Jacobian
    is taken at zero excitation)."""
    pos = np.array([np.pi / 3, 0.2, np.pi / 3, 0.4])
    solver, H = _diode_solver()
    kw = dict(H=H, frequency=8e9, Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol(), zero_pos=pos, return_complex=True)
    one = solver.linearised_N_spin_diode(Vdc_ex_value=1.0, **kw)
    big = solver.linearised_N_spin_diode(Vdc_ex_value=1e4, **kw)
    assert np.allclose(big, 1e4 * one, rtol=1e-6)


def test_spin_diode_zero_direction_rejected():
    with pytest.raises(ValueError):
        _diode_solver(VectorObj.from_cartesian(0, 0, 0))


def test_spin_diode_scan_matches_pointwise():
    """The vectorised scan equals linearised_N_spin_diode point by point (lab-frame drive, complex DL)."""
    solver, _ = _diode_solver(VectorObj.from_cartesian(-0.6, 0.8, 0), 1.5 - 0.7j)
    rng = np.random.default_rng(1)
    pos = np.column_stack([rng.uniform(1.2, 1.9, 4), rng.uniform(-1, 1, 4), rng.uniform(1.2, 1.9, 4), rng.uniform(-1, 1, 4)])
    H = np.column_stack([np.linspace(-50e3, 50e3, 4), np.full(4, 3e3), np.zeros(4)])
    kw = dict(frequency=6e9, Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol(), Vdc_ex_value=2.0)
    scan = solver.linearised_spin_diode_scan(H=H, zero_pos=pos, **kw)
    ref = [solver.linearised_N_spin_diode(H=VectorObj.from_cartesian(*h), zero_pos=x, return_complex=True, **kw).ravel()
           for h, x in zip(H, pos)]
    assert np.allclose(scan, ref, rtol=1e-8)


def test_spin_diode_mr_channels_first_order():
    """Channels equal the first-order change of AMR sum (m.j)^2 and GMR (1 - m1.m2)/2 for a small real step."""
    pos = np.array([[1.3, 0.4, 1.7, 2.5]])
    step = np.array([[2e-6, -1e-6, 3e-6, 1e-6]])
    j = np.array([np.cos(0.7), np.sin(0.7), 0.0])

    def R(x):
        t, p = x[0::2], x[1::2]
        m = np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)], -1)
        return np.array([np.sum((m @ j) ** 2), (1 - m[0] @ m[1]) / 2])

    expected = R(pos[0] + step[0]) - R(pos[0])
    assert np.allclose(spin_diode_mr_channels(pos, step, j)[0].real, expected, rtol=1e-4)


def test_spin_diode_scan_reuses_layout_compile():
    """Solvers with the same layout share one compiled system; each still gets its own parameters
    (including the dipole coupling and a complex DL ratio)."""
    from cmtj.models.general_sb import _DIODE_FN_CACHE

    rng = np.random.default_rng(2)
    pos = np.column_stack([rng.uniform(1.2, 1.9, 3), rng.uniform(-1, 1, 3), rng.uniform(1.2, 1.9, 3), rng.uniform(-1, 1, 3)])
    H = np.column_stack([np.linspace(-30e3, 30e3, 3), np.zeros(3), np.full(3, 1e3)])
    kw = dict(frequency=5e9, Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol())
    dipole = [VectorObj.from_cartesian(0, 1e-4, 0), VectorObj.from_cartesian(0, -1e-4, 0), VectorObj.from_cartesian(0, 0, 1e-6)]
    for k, (J1, ratio, Ku) in enumerate([(-2e-4, 0.3 + 1j, 5e3), (6e-4, -2.0 - 0.5j, 9e3)]):
        layers = [
            LayerDynamic(i, thickness=3e-9 + 1e-9 * i, Kv=VectorObj(np.deg2rad(90), 0.3 * i, Ku * (1 + i)), Ks=1e3,
                         Ms=1.5 / (4e-7 * np.pi), alpha=0.01 + 0.01 * i, Hdmi=[1e3 * i, -5e2, 0],
                         hoe_direction=[0.2, 1.0, 0.1], hoe_dl_ratio=ratio)
            for i in range(2)
        ]
        solver = Solver(layers=layers, J1=[J1], J2=[1e-5], ilD=[VectorObj.from_cartesian(1e-5, 0, 2e-5)], Ndipole=[dipole])
        n_cached = len(_DIODE_FN_CACHE)
        scan = solver.linearised_spin_diode_scan(H=H, zero_pos=pos, **kw)
        assert k == 0 or len(_DIODE_FN_CACHE) == n_cached  # second solver: no new compile
        ref = [solver.linearised_N_spin_diode(H=VectorObj.from_cartesian(*h), zero_pos=x, Vdc_ex_value=1.0,
                                              return_complex=True, **kw).ravel() for h, x in zip(H, pos)]
        assert np.allclose(scan, ref, rtol=1e-8)
