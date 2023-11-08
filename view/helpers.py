import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from cmtj import *
from cmtj.utils import FieldScan
from cmtj.utils.procedures import (PIMM_procedure, ResistanceParameters,
                                   VSD_procedure)


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


def prepare_simulation(
    Ms1,
    Ms2,
    K1,
    K2,
    alpha1,
    alpha2,
    thickness1,
    thickness2,
    width1,
    width2,
    length1,
    length2,
    anisotropy_axis1,
    anisotropy_axis2,
    J,
):
    demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1)]
    Kdir1 = get_axis_cvector(anisotropy_axis1)
    Kdir2 = get_axis_cvector(anisotropy_axis2)
    layer1 = Layer(
        "top",
        mag=Kdir1,
        anis=Kdir1,
        Ms=Ms1,
        thickness=thickness1,
        cellSurface=1e-16,
        demagTensor=demag,
        damping=alpha1,
    )
    layer2 = Layer(
        id="bottom",
        mag=Kdir2,
        anis=Kdir2,
        Ms=Ms2,
        damping=alpha2,
        demagTensor=demag,
        cellSurface=1e-16,
        thickness=thickness2,
    )
    j = Junction([layer1, layer2], 100, 200)
    j.setLayerAnisotropyDriver("top", ScalarDriver.getConstantDriver(K1))
    j.setLayerAnisotropyDriver("bottom", ScalarDriver.getConstantDriver(K2))
    j.setIECDriver("top", "bottom", ScalarDriver.getConstantDriver(J))

    rparams = [
        ResistanceParameters(
            Rxx0=100, Rxy0=1, Rsmr=-0.46, Rahe=-2.7, Ramr=-0.24, l=width1, w=length1
        ),
        ResistanceParameters(
            Rxx0=100, Rxy0=1, Rsmr=-0.46, Rahe=-2.7, Ramr=-0.24, l=width2, w=length2
        ),
    ]
    return j, rparams


@st.cache_data
def get_pimm_data(
    Ms1,
    Ms2,
    K1,
    K2,
    alpha1,
    alpha2,
    thickness1,
    thickness2,
    width1,
    width2,
    length1,
    length2,
    anisotropy_axis1,
    anisotropy_axis2,
    H_axis,
    Hmin,
    Hmax,
    Hsteps,
    J,
    int_step=1e-12,
    sim_time=16e-9,
):
    htheta, hphi = get_axis_angles(H_axis)
    Hscan, Hvecs = FieldScan.amplitude_scan(Hmin, Hmax, Hsteps, htheta, hphi)
    j, rparams = prepare_simulation(
        Ms1,
        Ms2,
        K1,
        K2,
        alpha1,
        alpha2,
        thickness1,
        thickness2,
        width1,
        width2,
        length1,
        length2,
        anisotropy_axis1,
        anisotropy_axis2,
        J,
    )
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


@st.cache_data
def get_vsd_data(
    Ms1,
    Ms2,
    K1,
    K2,
    alpha1,
    alpha2,
    thickness1,
    thickness2,
    width1,
    width2,
    length1,
    length2,
    anisotropy_axis1,
    anisotropy_axis2,
    H_axis,
    Hmin,
    Hmax,
    Hsteps,
    Hoex,
    J,
    fmin=0,
    fmax=30e9,
    nf=50,
    int_step=1e-12,
    sim_time=16e-9,
):
    htheta, hphi = get_axis_angles(H_axis)
    Hscan, Hvecs = FieldScan.amplitude_scan(Hmin, Hmax, Hsteps, htheta, hphi)
    j, rparams = prepare_simulation(
        Ms1,
        Ms2,
        K1,
        K2,
        alpha1,
        alpha2,
        thickness1,
        thickness2,
        width1,
        width2,
        length1,
        length2,
        anisotropy_axis1,
        anisotropy_axis2,
        J,
    )
    frequencies = np.linspace(fmin, fmax, nf)
    spec = VSD_procedure(
        j,
        Hvecs=Hvecs,
        int_step=int_step,
        frequencies=frequencies,
        resistance_params=rparams,
        simulation_duration=sim_time,
        Rtype="Rz",
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


def plot_data(Hscan, freqs, spec):
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
        ax.set_title("Resonance spectrum")

        try:
            fields, freqs = read_data()
            ax.plot(fields / 1e3, freqs / 1e9, "o", color="white", label="user data")
        except (ValueError, AttributeError):
            ...
        st.pyplot(fig)


def simulate_vsd():
    st.write("### VSD")
    with st.spinner("Simulating VSD..."):
        spec, freqs, Hscan = get_vsd_data(
            Ms1=st.session_state.Ms0,
            Ms2=st.session_state.Ms1,
            K1=st.session_state.K0 * 1e3,
            K2=st.session_state.K1 * 1e3,
            alpha1=st.session_state.alpha0,
            alpha2=st.session_state.alpha1,
            thickness1=st.session_state.thickness0 * 1e-9,
            thickness2=st.session_state.thickness1 * 1e-9,
            width1=st.session_state.width0 * 1e-6,
            width2=st.session_state.width1 * 1e-6,
            length1=st.session_state.length0 * 1e-6,
            length2=st.session_state.length1 * 1e-6,
            anisotropy_axis1=st.session_state.anisotropy_axis0,
            anisotropy_axis2=st.session_state.anisotropy_axis1,
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin * 1e3,
            Hmax=st.session_state.Hmax * 1e3,
            Hsteps=st.session_state.Hsteps,
            J=st.session_state.J * 1e-3,
            int_step=st.session_state.int_step,
            fmin=st.session_state.fmin * 1e9,
            fmax=st.session_state.fmax * 1e9,
            Hoex=st.session_state.Hoex,
            nf=st.session_state.nf,
        )
    plot_data(Hscan, freqs, spec)


def simulate_pimm():
    st.write("### PIMM")
    with st.spinner("Simulating PIMM..."):
        spec, freqs, _, Hscan = get_pimm_data(
            Ms1=st.session_state.Ms0,
            Ms2=st.session_state.Ms1,
            K1=st.session_state.K0 * 1e3,
            K2=st.session_state.K1 * 1e3,
            alpha1=st.session_state.alpha0,
            alpha2=st.session_state.alpha1,
            thickness1=st.session_state.thickness0 * 1e-9,
            thickness2=st.session_state.thickness1 * 1e-9,
            width1=st.session_state.width0 * 1e-6,
            width2=st.session_state.width1 * 1e-6,
            length1=st.session_state.length0 * 1e-6,
            length2=st.session_state.length1 * 1e-6,
            anisotropy_axis1=st.session_state.anisotropy_axis0,
            anisotropy_axis2=st.session_state.anisotropy_axis1,
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin * 1e3,
            Hmax=st.session_state.Hmax * 1e3,
            Hsteps=st.session_state.Hsteps,
            J=st.session_state.J * 1e-3,
            int_step=st.session_state.int_step,
        )

    plot_data(Hscan, freqs, spec)
