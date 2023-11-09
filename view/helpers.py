import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from cmtj import *
from cmtj.utils import FieldScan
from cmtj.utils.procedures import (PIMM_procedure, ResistanceParameters,
                                   VSD_procedure)


def create_single_layer(id_: str) -> tuple:
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    Kdir1 = get_axis_cvector(st.session_state[f"anisotropy_axis{id_}"])
    layer = Layer(
        id=f"layer_{id_}",
        mag=Kdir1,
        anis=Kdir1,
        Ms=st.session_state[f"Ms{id_}"],
        thickness=st.session_state[f"thickness{id_}"],
        cellSurface=1e-16,
        demagTensor=demag,
        damping=st.session_state[f"alpha{id_}"],
    )
    layer.setAnisotropyDriver(
        ScalarDriver.getConstantDriver(st.session_state[f"K{id_}"])
    )
    rp = ResistanceParameters(
        Rxx0=100,
        Rxy0=1,
        Rsmr=-0.46,
        Rahe=-2.7,
        Ramr=-0.24,
        l=st.session_state[f"length{id_}"],
        w=st.session_state[f"width{id_}"],
    )
    return layer, rp


def get_axis_cvector(axis: str):
    if axis == "x":
        return CVector(1, 0, 0)
    elif axis == "y":
        return CVector(0, 1, 0)
    elif axis == "z":
        return CVector(0, 0, 1)
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
    if axis == "x":
        return 0, 0
    elif axis == "y":
        return 0, 90
    elif axis == "z":
        return 90, 0
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
        J = st.session_state[f"J{jvals}"]
        l1_name = layers[jvals].id
        l2_name = layers[jvals + 1].id
        j.setIECDriver(l1_name, l2_name, ScalarDriver.getConstantDriver(J))
    return j, rparams


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
    Hscan, Hvecs = FieldScan.amplitude_scan(Hmin, Hmax, Hsteps, htheta, hphi)
    j, rparams = prepare_simulation()
    spec, freqs, out = PIMM_procedure(
        j,
        Hvecs=Hvecs,
        int_step=int_step,
        resistance_params=rparams,
        max_frequency=60e9,
        simulation_duration=sim_time,
        disturbance=1e-6,
        wait_time=4e-9,
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
    nf=50,
    int_step=1e-12,
    sim_time=16e-9,
):
    htheta, hphi = get_axis_angles(H_axis)
    Hscan, Hvecs = FieldScan.amplitude_scan(Hmin, Hmax, Hsteps, htheta, hphi)
    j, rparams = prepare_simulation()
    frequencies = np.linspace(fmin, fmax, nf)
    spec = VSD_procedure(
        j,
        Hvecs=Hvecs,
        int_step=int_step,
        frequencies=frequencies,
        resistance_params=rparams,
        simulation_duration=sim_time,
        Rtype="Rx",
        Hoe_excitation=500,
        Hoe_direction=get_axis(Hoex),
        disturbance=1e-6,
    )
    return spec, frequencies, Hscan


def read_data():
    filedata = st.session_state.upload.read().decode("utf-8")
    lines = filedata.split("\n")
    fields, freqs = [], []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        fields.append(float(line.split()[0]))
        freqs.append(float(line.split()[1]))
    return np.asarray(fields), np.asarray(freqs)


def plot_data(Hscan, freqs, spec, title="Resonance spectrum"):
    with plt.style.context(["dark_background"]):
        fig, ax = plt.subplots(dpi=300)
        ax.pcolormesh(
            Hscan / 1e3,
            freqs / 1e9,
            10 * np.log10(np.abs(spec.T)),
            shading="auto",
            cmap="inferno",
            rasterized=True,
        )
        ax.set_xlabel("H (kA/m)")
        ax.set_ylabel("Frequency (GHz)")
        ax.set_title(title)

        try:
            fields, freqs = read_data()
            ax.plot(fields / 1e3, freqs / 1e9, "o", color="white", label="user data")
        except (ValueError, AttributeError):
            ...
        st.pyplot(fig)


def simulate_vsd():
    with st.spinner("Simulating VSD..."):
        spec, freqs, Hscan = get_vsd_data(
            int_step=st.session_state.int_step,
            fmin=st.session_state.fmin * 1e9,
            fmax=st.session_state.fmax * 1e9,
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin * 1e3,
            Hmax=st.session_state.Hmax * 1e3,
            Hsteps=st.session_state.Hsteps,
            Hoex=st.session_state.Hoex,
            nf=st.session_state.nf,
        )
    plot_data(Hscan, freqs, spec, title="VSD spectrum")


def simulate_pimm():
    st.write("### PIMM")
    with st.spinner("Simulating PIMM..."):
        spec, freqs, _, Hscan = get_pimm_data(
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin * 1e3,
            Hmax=st.session_state.Hmax * 1e3,
            Hsteps=st.session_state.Hsteps,
            int_step=st.session_state.int_step,
        )

    plot_data(Hscan, freqs, spec, title="PIMM spectrum")
