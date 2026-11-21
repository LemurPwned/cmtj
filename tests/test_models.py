import numpy as np
import pytest
from cmtj.models import Solver, LayerDynamic, VectorObj, LayerSB, spin_diode_mr_channels
from pathlib import Path
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


def test_hoe_direction_needs_three_components():
    with pytest.raises(ValueError, match="3 components"):
        LayerDynamic(0, thickness=1e-9, Kv=VectorObj(0, 0, 1e3), Ks=0, Ms=1e6, hoe_direction=[0, 1])


def test_spin_diode_scan_numpy_complex_ratio():
    """numpy complex scalars (not `complex` subclasses) are passed as complex parameters."""
    pos = np.array([[1.5, 0.3]])
    H = np.array([[1e4, 0, 1e3]])
    kw = dict(frequency=5e9, Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol())
    scans = []
    for ratio in (np.complex64(0.5 + 1j), 0.5 + 1j):
        layer = LayerDynamic(0, thickness=3e-9, Kv=VectorObj(np.deg2rad(90), 0, 5e3), Ks=1e3, Ms=1.2e6,
                             hoe_direction=[0, 1, 0], hoe_dl_ratio=ratio)
        scans.append(Solver(layers=[layer], J1=[], J2=[]).linearised_spin_diode_scan(H=H, zero_pos=pos, **kw))
    assert np.allclose(scans[0], scans[1], rtol=1e-6)


def test_diode_cache_is_bounded_and_clearable(monkeypatch):
    from cmtj.models import general_sb as gsb

    monkeypatch.setattr(gsb, "_DIODE_FN_CACHE_SIZE", 1)
    gsb.clear_diode_cache()
    pos, H = np.array([[1.5, 0.3]]), np.array([[1e4, 0, 1e3]])
    for hoe in (None, [0, 1, 0]):  # two layouts
        layer = LayerDynamic(0, thickness=3e-9, Kv=VectorObj(np.deg2rad(90), 0, 5e3), Ks=1e3, Ms=1.2e6,
                             hoe_direction=hoe)
        Solver(layers=[layer], J1=[], J2=[]).linearised_spin_diode_scan(
            H=H, zero_pos=pos, frequency=5e9, Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol())
    assert len(gsb._DIODE_FN_CACHE) == 1
    gsb.clear_diode_cache()
    assert len(gsb._DIODE_FN_CACHE) == 0


def _one_layer_scan(**layer_kw):
    layer = LayerDynamic(0, thickness=3e-9, Kv=VectorObj(np.deg2rad(90), 0, 5e3), Ks=1e3, Ms=1.2e6, **layer_kw)
    solver = Solver(layers=[layer], J1=[], J2=[])
    kw = dict(H=np.array([[1e4, 0, 1e3]]), zero_pos=np.array([[1.5, 0.3]]), frequency=5e9,
              Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol())
    return solver, kw


def test_hoe_direction_read_at_use_time():
    """Reassigning hoe_direction after construction is used by the next scan (no stale snapshot)."""
    solver, kw = _one_layer_scan(hoe_direction=[1, 0, 0])
    solver.layers[0].hoe_direction = [0, 1, 0.3]
    fresh, _ = _one_layer_scan(hoe_direction=[0, 1, 0.3])
    assert np.allclose(solver.linearised_spin_diode_scan(**kw), fresh.linearised_spin_diode_scan(**kw))
    solver.layers[0].hoe_direction = [0, 1]
    with pytest.raises(ValueError, match="3 components"):
        solver.linearised_spin_diode_scan(**kw)


def test_spin_diode_scan_follows_gamma():
    """Constants are part of the compile key: changing gamma after a compile is not served stale."""
    from cmtj.utils import Constants

    gamma = Constants.gamma_rad()
    solver, kw = _one_layer_scan(hoe_direction=[0, 1, 0])
    base = solver.linearised_spin_diode_scan(**kw)
    try:
        Constants.set_gamma_rad(1.5 * gamma)
        scan = solver.linearised_spin_diode_scan(**kw)
        ref = solver.linearised_N_spin_diode(H=VectorObj.from_cartesian(*kw["H"][0]), zero_pos=kw["zero_pos"][0],
                                             frequency=kw["frequency"], Vdc_ex_variable=kw["Vdc_ex_variable"],
                                             Vdc_ex_value=1.0, return_complex=True).ravel()
    finally:
        Constants.set_gamma_rad(gamma)
    assert not np.allclose(scan, base)
    assert np.allclose(scan.ravel(), ref, rtol=1e-6)


def test_spin_diode_hoe_amplitude():
    """hoe_amplitude scales a layer's drive (field-like and damping-like); response is linear in it per layer."""
    pos = np.array([np.pi / 3, 0.2, np.pi / 2.5, 0.9])
    d = VectorObj.from_cartesian(-0.6, 0.8, 0)

    def solver(a0, a1):
        s, H = _diode_solver(d, 0.7 - 0.3j)
        for layer, a in zip(s.layers, (a0, a1)):
            layer.hoe_amplitude = a
        return _diode(s, H, pos)

    base = solver(1, 1)
    assert np.allclose(solver(1, 1), _diode(*_diode_solver(d, 0.7 - 0.3j), pos))
    assert np.allclose(solver(2.5, 2.5), 2.5 * base, rtol=1e-6)
    # superposition: drive on layer 1 only + layer 0 only = both
    assert np.allclose(solver(1, 0) + solver(0, 1.0 + 2j), solver(1, 1.0 + 2j), rtol=1e-6, atol=1e-12)


def test_spin_diode_mr_channels_per_layer_amr_smr():
    """Per-layer AMR/SMR channel equals the first-order change of sum_i AMR_i (m_i.j)^2 + SMR_i (m_i.y)^2, y = z x j,
    i.e. compute_resistance's R_xx with x along j; defaults reproduce the equal-weight AMR channel."""
    pos = np.array([[1.3, 0.4, 1.7, 2.5]])
    step = np.array([[2e-6, -1e-6, 3e-6, 1e-6]])
    c = 0.7
    j = np.array([np.cos(c), np.sin(c), 0.0])
    y = np.array([-np.sin(c), np.cos(c), 0.0])
    amr, smr = [0.8, 1.3], [0.0, -0.4]

    def R(x):
        t, p = x[0::2], x[1::2]
        m = np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)], -1)
        return np.sum(np.array(amr) * (m @ j) ** 2 + np.array(smr) * (m @ y) ** 2)

    got = spin_diode_mr_channels(pos, step, j, AMR=amr, SMR=smr)[0, 0].real
    assert np.isclose(got, R(pos[0] + step[0]) - R(pos[0]), rtol=1e-4)
    # defaults: the equal-weight AMR channel sum_i 2 (m_i.j)(dm_i.j), computed independently
    t, p = pos[0, 0::2], pos[0, 1::2]
    m = np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)], -1)
    dt, dp = step[0, 0::2], step[0, 1::2]
    theta_hat = np.stack([np.cos(t) * np.cos(p), np.cos(t) * np.sin(p), -np.sin(t)], -1)
    phi_hat_sin = np.stack([-np.sin(t) * np.sin(p), np.sin(t) * np.cos(p), 0 * t], -1)  # dm/dphi
    dm = dt[:, None] * theta_hat + dp[:, None] * phi_hat_sin
    assert np.isclose(spin_diode_mr_channels(pos, step, j)[0, 0], np.sum(2 * (m @ j) * (dm @ j)), rtol=1e-12)


def test_spin_diode_scan_hoe_amplitude_no_recompile(monkeypatch):
    """hoe_amplitude is an argument of the compiled scan: changing it on the same solver rescales the
    response without a new compile, and per-layer amplitudes match the direct (symbolic) solution."""
    from cmtj.models import general_sb as gsb

    compiles = []
    lambdify = gsb.sym.lambdify
    monkeypatch.setattr(gsb.sym, "lambdify", lambda *a, **k: compiles.append(1) or lambdify(*a, **k))
    pos = np.array([[np.pi / 3, 0.2, np.pi / 2.5, 0.9]])
    s, H = _diode_solver(VectorObj.from_cartesian(-0.6, 0.8, 0), 0.7 - 0.3j)
    kw = dict(H=np.array([H.get_cartesian()]), zero_pos=pos, frequency=8e9, Vdc_ex_variable=LayerDynamic.get_hoe_ex_symbol())
    base = s.linearised_spin_diode_scan(**kw)
    n_compiles = len(compiles)
    for layer, a in zip(s.layers, (2 - 1j, 0.5)):
        layer.hoe_amplitude = a
    scan = s.linearised_spin_diode_scan(**kw)
    assert not np.allclose(scan, base)
    assert np.allclose(scan.ravel(), _diode(s, H, pos[0]), rtol=1e-6)
    for layer in s.layers:
        layer.hoe_amplitude = 3.0
    assert np.allclose(s.linearised_spin_diode_scan(**kw), 3.0 * base, rtol=1e-6)
    assert len(compiles) == n_compiles  # amplitude changes never recompile


def test_spin_diode_mr_channels_validation():
    pos = np.array([[1.3, 0.4, 1.7, 2.5]])
    step = np.array([[2e-6, -1e-6, 3e-6, 1e-6]])
    with pytest.raises(ValueError, match="one value per layer"):
        spin_diode_mr_channels(pos, step, [1, 0, 0], AMR=[1.0])
    with pytest.raises(ValueError, match="one value per layer"):
        spin_diode_mr_channels(pos, step, [1, 0, 0], SMR=[0.1, 0.1, 0.1])
    with pytest.raises(ValueError, match="along z"):
        spin_diode_mr_channels(pos, step, [0, 0, 1], SMR=[0.0, -0.4])
    # current along z without SMR: AMR channel only, finite
    assert np.all(np.isfinite(spin_diode_mr_channels(pos, step, [0, 0, 1])))


@pytest.mark.parametrize("field", ["hoe_amplitude", "hoe_dl_ratio", "hoe_direction"])
def test_drive_parameters_reject_strings(field):
    """Strings would go through sympy's parser (which evaluates them): rejected up front and at use time."""
    kw = dict(hoe_direction=[0, 1, 0])
    kw[field] = "__import__('os')" if field != "hoe_direction" else ["1", "0", "0"]
    with pytest.raises(TypeError, match="must be a number"):
        LayerDynamic(0, thickness=1e-9, Kv=VectorObj(0, 0, 1e3), Ks=0, Ms=1e6, **kw)
    s, H = _diode_solver(VectorObj.from_cartesian(0, 1, 0))
    s.layers[0].hoe_amplitude = "1"
    with pytest.raises(TypeError, match="must be a number"):
        _diode(s, H, np.array([np.pi / 3, 0.2, np.pi / 2.5, 0.9]))


def _coruco_solver(params, ms=1.65 / (4e-7 * np.pi), thickness=4e-9, ndemag=(2.9e-05, 1.33e-4, 0.99984)):
    """Co/Ru/Co-like bilayer: in-plane uniaxial axes, Hdmi, J1/J2 and interlayer DMI."""
    Ku, Kphi, asym, dphi, J1, J2, Dx, Dy, Dz, h1x, h1y, h2x, h2y = params
    layers = [
        LayerSB(_id=i, thickness=thickness, Kv=VectorObj(0, Kphi + (0.5 - i) * dphi, Ku * (1 + (1 - 2 * i) * asym)),
                Ks=0, Ms=ms, Ndemag=VectorObj.from_cartesian(*ndemag), Hdmi=[(h1x, h2x)[i], (h1y, h2y)[i], 0])
        for i in range(2)
    ]
    return Solver(layers=layers, J1=[J1], J2=[J2], ilD=[VectorObj.from_cartesian(Dx, Dy, Dz)])


def _m(coords):
    t, p = coords[..., 0::2], coords[..., 1::2]
    return np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)], -1)


def test_relax_field_path_matches_per_point_solver():
    """Each equilibrium of the path is a minimum the per-point Adam solver agrees with (same history)."""
    params = (2e3, 0.3, -0.2, 0.2, -1.2e-4, -5e-5, 0, 0, 3e-5, 2e3, 0, -2e3, 0)
    solver = _coruco_solver(params)
    path = np.linspace(-60e3, 60e3, 25)
    eq = solver.relax_field_path(path, np.deg2rad([20.0]))
    for k in (0, 8, 12, 16, 24):
        solver.set_H(VectorObj(np.pi / 2, np.deg2rad(20.0), path[k]))
        adam = solver.adam_gradient_descent(eq[k, 0], max_steps=200000, tol=1e-12, learning_rate=1e-4, perturbation=0)
        assert np.allclose(_m(adam), _m(eq[k, 0]), atol=1e-4)


def test_relax_field_path_matches_reference_hysteresis():
    """Hysteresis loops of the vectorised Newton relaxation used in the Co/Ru/Co polar-coercivity fits."""
    ref = np.load(Path(__file__).parent / "assets" / "sb_relax_path_reference.npz")
    solver = _coruco_solver(ref["params"], float(ref["ms"]), float(ref["thickness"]), ref["ndemag"])
    eq = solver.relax_field_path(ref["path"], np.deg2rad(ref["angles_deg"]))
    assert np.allclose(_m(eq), ref["m"], atol=1e-5)


def test_relax_field_path_new_parameters_no_recompile(monkeypatch):
    from cmtj.models import general_sb as gsb

    compiles = []
    kernel = gsb._row_kernel
    monkeypatch.setattr(gsb, "_row_kernel", lambda *a, **k: compiles.append(1) or kernel(*a, **k))
    path, angles = np.linspace(-50e3, 50e3, 11), np.deg2rad([0.0, 90.0])
    a = _coruco_solver((2e3, 0.3, 0, 0, -1e-4, 0, 0, 0, 0, 0, 0, 0, 0)).relax_field_path(path, angles)
    n_compiles = len(compiles)
    b, hess = _coruco_solver((5e3, 1.0, 0.1, 0, -3e-4, -1e-4, 0, 0, 1e-5, 1e3, 0, 0, 0)).relax_field_path(
        path, angles, return_hessian=True)
    assert len(compiles) == n_compiles
    assert not np.allclose(a, b)
    assert hess.shape == (11, 2, 4, 4) and np.allclose(hess, np.swapaxes(hess, -1, -2))


def test_relax_field_path_rejects_zero_direction():
    with pytest.raises(ValueError, match="non-zero"):
        _coruco_solver((2e3, 0.3, 0, 0, -1e-4, 0, 0, 0, 0, 0, 0, 0, 0)).relax_field_path([1e3], [[1, 0, 0], [0, 0, 0]])


def test_relax_field_path_warns_when_unconverged():
    solver = _coruco_solver((2e3, 0.3, 0, 0, -1e-4, 0, 0, 0, 0, 0, 0, 0, 0))
    with pytest.warns(RuntimeWarning, match="did not converge"):
        solver.relax_field_path(np.linspace(-50e3, 50e3, 5), np.deg2rad([10.0]), max_iter=1)


def test_relax_field_path_follows_mu0_change():
    from cmtj.models import clear_diode_cache
    from cmtj.utils import Constants

    params = (2e3, 0.3, 0, 0, -1e-4, 0, 0, 0, 0, 0, 0, 0, 0)
    path, angles = np.linspace(-50e3, 50e3, 11), np.deg2rad([10.0, 80.0])
    mu0 = Constants.mu0()
    try:
        _coruco_solver(params).relax_field_path(path, angles)  # compile and cache at the default mu0
        Constants.set_mu0(2 * mu0)
        changed = _coruco_solver(params).relax_field_path(path, angles)
        clear_diode_cache()
        Solver._create_energy.cache_clear()
        fresh = _coruco_solver(params).relax_field_path(path, angles)
    finally:
        Constants.set_mu0(mu0)
    assert np.allclose(changed, fresh)


@pytest.mark.parametrize(
    "path, directions, match",
    [([1e3], np.ones((2, 2)), "directions must be"), ([1e3], np.ones((1, 2, 3)), "directions must be"),
     ([], [0.0], "Hpath must be"), ([[1e3]], [0.0], "Hpath must be")],
)
def test_relax_field_path_rejects_bad_shapes(path, directions, match):
    with pytest.raises(ValueError, match=match):
        _coruco_solver((2e3, 0.3, 0, 0, -1e-4, 0, 0, 0, 0, 0, 0, 0, 0)).relax_field_path(path, directions)
