from collections import defaultdict
from itertools import groupby

from new_sb import LayerDynamic
import numpy as np
import streamlit as st

from cmtj import Layer, CVector, constantDriver, Junction, Axis, AxialDriver
from cmtj.models import LayerSB, Solver
from cmtj.utils import FieldScan, VectorObj, mu0
from cmtj.utils.procedures import PIMM_procedure, ResistanceParameters, VSD_procedure


def create_single_domain(id_: str) -> Layer:
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    # Kdir1 = get_axis_cvector(st.session_state[f"anisotropy_axis{id_}"])
    Kdir = FieldScan.angle2vector(theta=st.session_state[f"theta_K{id_}"], phi=st.session_state[f"phi_K{id_}"])
    layer = Layer(
        id=f"domain_{id_}",
        mag=Kdir,
        anis=Kdir,
        Ms=st.session_state["Ms_shared"],
        thickness=st.session_state["thickness_shared"] * 1e-9,
        cellSurface=1e-16,
        demagTensor=demag,
        damping=st.session_state["alpha_shared"],
    )
    layer.setAnisotropyDriver(constantDriver(st.session_state[f"K{id_}"] * 1e3))
    return layer


def create_single_layer(id_: str) -> tuple:
    """Do not forget to rescale the units!"""
    nxx = st.session_state[f"Nxx{id_}"]
    nyy = st.session_state[f"Nyy{id_}"]
    nzz = st.session_state[f"Nzz{id_}"]
    demag = [
        CVector(nxx, 0, 0),
        CVector(0, nyy, 0),
        CVector(0, 0, nzz),
    ]
    demag_sum = nxx + nyy + nzz
    if abs(demag_sum - 1.0) > 1e-5:
        st.warning(f"Warning: Demagnetization tensor components should sum to 1.0 (Layer {id_})")
    Kdir = FieldScan.angle2vector(theta=st.session_state[f"theta_K{id_}"], phi=st.session_state[f"phi_K{id_}"])
    layer = Layer(
        id=f"layer_{id_}",
        mag=Kdir,
        anis=Kdir,
        Ms=st.session_state[f"Ms{id_}"],
        thickness=st.session_state[f"thickness{id_}"] * 1e-9,
        cellSurface=1e-16,
        demagTensor=demag,
        damping=st.session_state[f"alpha{id_}"],
    )
    layer.setAnisotropyDriver(constantDriver(st.session_state[f"K{id_}"] * 1e3))
    rp = ResistanceParameters(
        Rxx0=100,
        Rxy0=1,
        Rsmr=-0.46,
        Rahe=-2.7,
        Ramr=-0.24,
        l=st.session_state[f"length{id_}"] * 1e-6,
        w=st.session_state[f"width{id_}"] * 1e-6,
    )
    return layer, rp


def create_sb_layer(id_: str) -> tuple[LayerDynamic, list[float]]:
    nxx = st.session_state[f"Nxx{id_}"]
    nyy = st.session_state[f"Nyy{id_}"]
    nzz = st.session_state[f"Nzz{id_}"]
    Ks = st.session_state[f"Ks{id_}"]
    Kv = st.session_state[f"Kv{id_}"]
    kphi = st.session_state[f"phi_K{id_}"]
    demag = VectorObj.from_cartesian(nxx, nyy, nzz)
    layer = LayerDynamic(
        _id=int(id_),
        thickness=st.session_state[f"thickness{id_}"] * 1e-9,
        Kv=VectorObj(np.deg2rad(0.0), np.deg2rad(kphi), Kv),
        Ks=Ks,
        Ms=st.session_state[f"Ms{id_}"] / mu0,
        demagTensor=demag,
        damping=st.session_state[f"alpha{id_}"],
    )
    ktheta = 90 if Kv > Ks else 0  # if Kv is smaller than Ks, assume in plane
    return layer, [np.deg2rad(ktheta), np.deg2rad(kphi)]


def get_axis_cvector(axis: str):
    if axis == "x":
        return CVector(1, 0, 0)
    elif axis == "y":
        return CVector(0, 1, 0)
    elif axis == "z":
        return CVector(0, 0, 1)
    elif axis == "xy":
        return CVector(1, 1, 0)
    elif axis == "xz":
        return CVector(1, 0, 1)
    elif axis == "yz":
        return CVector(0, 1, 1)
    else:
        raise ValueError(f"Invalid axis {axis}")


def get_axis(axis: str) -> Axis:
    if axis == "x":
        return Axis.xaxis
    elif axis == "y":
        return Axis.yaxis
    elif axis == "z":
        return Axis.zaxis
    else:
        raise ValueError(f"Invalid axis {axis}")


def get_axis_angles(axis: str):
    """Returns (theta, phi)"""
    if axis == "x":
        return 90, 0
    elif axis == "y":
        return 90, 90
    elif axis == "z":
        return 0, 0
    elif axis == "xy":
        return 90, 45
    elif axis == "xz":
        return 45, 0
    elif axis == "yz":
        return 45, 90
    else:
        raise ValueError(f"Invalid axis {axis}")


def prepare_simulation():
    layers = []
    rparams = []
    N = st.session_state["N"]
    for i in range(N):
        layer, rp = create_single_layer(i)
        layers.append(layer)
        rparams.append(rp)
    j = Junction(layers=layers)
    for jvals in range(N - 1):
        J = st.session_state[f"J{jvals}"] * 1e-6  # rescale GUI units
        J2 = st.session_state[f"J2{jvals}"] * 1e-6  # rescale GUI units
        ilD = st.session_state[f"ilD{jvals}"] * 1e-6  # rescale GUI units
        l1_name = layers[jvals].id
        l2_name = layers[jvals + 1].id
        j.setIECDriver(l1_name, l2_name, constantDriver(J))
        j.setQuadIECDriver(l2_name, l1_name, constantDriver(J2))
        j.setIDMIDriver(l1_name, l2_name, AxialDriver(0, 0, ilD))
    return j, rparams


def prepare_sb_simulation() -> tuple:
    layers = []
    init_pos = []
    N = st.session_state["N"]
    for i in range(N):
        layer, init_pos_i = create_sb_layer(i)
        layers.append(layer)
        init_pos.append(init_pos_i)
    Js = []
    Js2 = []
    for jvals in range(N - 1):
        J = st.session_state[f"J{jvals}"] * 1e-6  # rescale GUI units
        J2 = st.session_state[f"J2{jvals}"] * 1e-6  # rescale GUI units
        Js.append(J)
        Js2.append(J2)
    return layers, init_pos, Js, Js2


def get_spectrum_sb_data(H_axis, Hmin, Hmax, Hsteps, run_vsd: bool = False):
    layers, init_pos, Js, Js2 = prepare_sb_simulation()

    htheta, hphi = get_axis_angles(H_axis)
    hmin, hmax = min([Hmin, Hmax]), max([Hmin, Hmax])  # fix user input
    _, Hvecs = FieldScan.amplitude_scan(hmin, hmax, Hsteps, htheta, hphi)
    force_single_layer = not any(Js) and not any(Js2)
    all_sb_data = defaultdict(list)
    for H in Hvecs:
        Hvec = VectorObj.from_cartesian(*H)
        solver = Solver(layers=layers, J1=Js, J2=Js2, H=Hvec)
        eq, frequencies = solver.solve(
            init_position=init_pos,
            perturbation=1e-4,
            force_single_layer=force_single_layer,
        )
        for freq in frequencies:
            all_sb_data["Hmag"].append(Hvec.mag / 1e3)
            all_sb_data["frequency"].append(freq)

        if run_vsd:
            res = solver.linearised_N_spin_diode(
                H=Hvec,
                frequency=freq * 1e9,
                Vdc_ex_variable=LayerDynamic.get_Vp_symbol(),
                Vdc_ex_value=st.session_state["Hoe_mag"] * 1e3,
                zero_pos=eq.tolist(),
                phase_shift=0,
                cache_var="H",
            )

        all_sb_data["pos"].append(eq.tolist())
        init_pos = eq.tolist()
    return all_sb_data


# @st.cache_data
def get_pimm_data(
    H_axis,
    Hmin,
    Hmax,
    Hsteps,
    int_step=1e-12,
    sim_time=16e-9,
):
    htheta, hphi = get_axis_angles(H_axis)
    hmin, hmax = min([Hmin, Hmax]), max([Hmin, Hmax])  # fix user input
    Hscan, Hvecs = FieldScan.amplitude_scan(hmin, hmax, Hsteps, htheta, hphi)
    if st.session_state["Hreturn"]:
        Hscan = np.concatenate((Hscan, Hscan[::-1]))
        Hvecs = np.concatenate((Hvecs, Hvecs[::-1]))
    j, rparams = prepare_simulation()
    # avoid wait if the user sim's too short
    wtime = 4e-9 if sim_time >= 6e-9 else 0.0
    spec, freqs, out = PIMM_procedure(
        j,
        Hvecs=Hvecs,
        Hoe_direction=get_axis(st.session_state["Hoeaxis"]),
        Hoe_excitation=st.session_state["Hoe_mag"] * 1e3,
        int_step=int_step,
        resistance_params=rparams,
        max_frequency=st.session_state["max_freq"] * 1e9,
        simulation_duration=sim_time,
        disturbance=1e-6,
        wait_time=wtime,
    )
    return spec, freqs, out, Hscan


def get_vsd_data(
    H_axis,
    Hmin,
    Hmax,
    Hsteps,
    Hoex,
    fmin=0,
    fmax=30e9,
    fstep=0.5e9,
    int_step=1e-12,
    sim_time=16e-9,
    Rtype="Ry",
    Hoex_mag=500,
):
    htheta, hphi = get_axis_angles(H_axis)
    hmin, hmax = min([Hmin, Hmax]), max([Hmin, Hmax])  # fix user input
    Hscan, Hvecs = FieldScan.amplitude_scan(hmin, hmax, Hsteps, htheta, hphi)
    j, rparams = prepare_simulation()
    frequencies = np.arange(fmin, fmax, step=fstep)
    spec = VSD_procedure(
        j,
        Hvecs=Hvecs,
        int_step=int_step,
        frequencies=frequencies,
        resistance_params=rparams,
        simulation_duration=sim_time,
        Rtype=Rtype,
        Hoe_excitation=Hoex_mag,
        Hoe_direction=get_axis(Hoex),
        disturbance=1e-6,
    )
    return spec, frequencies, Hscan


def compute_sb_mse(target, data):
    # if there are multiple layers, and target is a single list
    # we need to take the closest resonance
    if len(data["frequency"]) < len(target):
        return float("inf")

    mse = 0
    for i, (_, f) in enumerate(groupby(zip(data["Hmag"], data["frequency"]), key=lambda x: x[0])):
        # f is a list of tuples (Hmag, frequency)
        contrib = min(((target[i] / 1e9 - f_[1]) ** 2).sum() for f_ in f)
        mse += contrib

    return mse


def simulate_sb(hvals: list[float]):
    layers = []
    N = st.session_state.N
    init_pos = []
    for i in range(N):
        ktheta, kphi = get_axis_angles(st.session_state[f"anisotropy_axis{i}"])
        Kval = st.session_state[f"K{i}"] * 1e3  # rescale GUI units
        if ktheta == 0:
            Ks = Kval
            Kv = 10
        else:
            Ks = 10
            Kv = Kval
        layer = LayerSB(
            _id=i,
            thickness=st.session_state[f"thickness{i}"] * 1e-9,  # rescale GUI units
            Kv=VectorObj(np.deg2rad(0.0), np.deg2rad(kphi), Kv),
            Ks=Ks,
            Ms=st.session_state[f"Ms{i}"] / mu0,
            # alpha=st.session_state[f"alpha{i}"],
        )
        init_pos.extend([np.deg2rad(ktheta), np.deg2rad(kphi)])
        layers.append(layer)
    Js = [st.session_state[f"J{i}"] * 1e-6 for i in range(N - 1)]
    result_dictionary = defaultdict(list)
    # we perform a sweep over the field magnitude
    htheta, hphi = get_axis_angles(st.session_state.H_axis)
    for i in range(len(hvals)):
        solver = Solver(
            layers=layers,
            J1=Js,
            J2=[0 for _ in Js],
            H=VectorObj(theta=htheta, phi=hphi, mag=hvals[i]),
        )
        eq, frequencies = solver.solve(init_position=init_pos, perturbation=1e-4)
        for freq in frequencies:
            result_dictionary["Hmag"].append(hvals[i] / 1e3)
            result_dictionary["frequency"].append(freq)
        init_pos = eq

    return result_dictionary


def kwargs_to_list(kwargs: dict, N: int):
    return {
        "J": [kwargs[f"J{i}"] for i in range(N - 1)],
        "K": [kwargs[f"K{i}"] for i in range(N)],
        "Ms": [kwargs[f"Ms{i}"] for i in range(N)],
    }


def simulate_sb_wrapper(
    hvals: list[float],
    N: int,
    thickness: list[float],
    anisotropy_axis: list[str],
    H_axis: str,
    **kwargs,
):
    J, K, Ms = kwargs_to_list(kwargs, N).values()
    layers = []
    init_pos = []
    for i in range(N):
        ktheta, kphi = get_axis_angles(anisotropy_axis[i])
        Kval = K[i] * 1e3
        if ktheta == 0:
            Ks = Kval
            Kv = 1
        else:
            Ks = 1
            Kv = Kval
        layer = LayerSB(
            _id=i,
            thickness=thickness[i] * 1e-9,  # rescale GUI units
            Kv=VectorObj(np.deg2rad(0.0), np.deg2rad(kphi), Kv),
            Ks=Ks,
            Ms=Ms[i] / mu0,
        )
        init_pos.extend([np.deg2rad(ktheta), np.deg2rad(kphi)])
        layers.append(layer)
    Js = [J[i] * 1e-6 for i in range(N - 1)]
    result_dictionary = defaultdict(list)
    # we perform a sweep over the field magnitude
    htheta, hphi = get_axis_angles(H_axis)
    for i in range(len(hvals)):
        solver = Solver(
            layers=layers,
            J1=Js,
            J2=[0 for _ in Js],
            H=VectorObj(theta=htheta, phi=hphi, mag=hvals[i]),
        )
        eq, frequencies = solver.solve(init_position=init_pos, perturbation=1e-4)
        for freq in frequencies:
            result_dictionary["Hmag"].append(hvals[i] / 1e3)
            result_dictionary["frequency"].append(freq)
        init_pos = eq

    return result_dictionary
