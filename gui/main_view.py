import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from cmtj import *
from cmtj.utils import FieldScan
from cmtj.utils.procedures import PIMM_procedure, ResistanceParameters

apptitle = "CMTJ simulator"
st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")
st.title(apptitle)


def get_axis_cvector(axis: str):
    if axis == "x":
        return CVector(1, 0, 0)
    elif axis == "y":
        return CVector(0, 1, 0)
    elif axis == "z":
        return CVector(0, 0, 1)
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
    J,
    int_step=1e-12,
    sim_time=16e-9,
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
    j = Junction([layer1, layer2])
    j.setLayerAnisotropyDriver("top", ScalarDriver.getConstantDriver(K1))
    j.setLayerAnisotropyDriver("bottom", ScalarDriver.getConstantDriver(K2))
    j.setIECDriver("top", "bottom", ScalarDriver.getConstantDriver(J))
    htheta, hphi = get_axis_angles(H_axis)
    Hscan, Hvecs = FieldScan.amplitude_scan(Hmin, Hmax, 100, htheta, hphi)
    rparams = [
        ResistanceParameters(
            Rxx0=100, Rxy0=1, Rsmr=-0.46, Rahe=-2.7, Ramr=-0.24, l=width1, w=length1
        ),
        ResistanceParameters(
            Rxx0=100, Rxy0=1, Rsmr=-0.46, Rahe=-2.7, Ramr=-0.24, l=width2, w=length2
        ),
    ]

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


N = 2
st.markdown(
    """
    This app simulates the resonance characteristics of a MTJ device.
    The device is composed of two layers, each with its own magnetic properties.
    """
)


with st.sidebar:
    st.markdown("## Simulation parameters")
    st.markdown("### Layer parameters")
    for i in range(N):
        st.slider(
            f"Ms ({i})", min_value=0.2, max_value=2.0, value=1.0, step=0.1, key=f"Ms{i}"
        )
        st.number_input(
            f"K ({i})", min_value=0.1, max_value=100e6, value=50e3, key=f"K{i}"
        )
        st.number_input(
            f"alpha ({i})", min_value=1e-5, max_value=0.1, value=0.001, key=f"alpha{i}"
        )
        st.number_input(
            f"thickness ({i})",
            min_value=1e-9,
            max_value=10e-9,
            value=1e-9,
            key=f"thickness{i}",
        )
        st.number_input(
            f"width ({i})",
            min_value=1e-6,
            max_value=500e-6,
            value=10e-6,
            key=f"width{i}",
        )
        st.number_input(
            f"length ({i})",
            min_value=1e-6,
            max_value=500e-6,
            value=10e-6,
            key=f"length{i}",
        )
        st.radio(
            f"anisotropy axis ({i})", options=["x", "y", "z"], key=f"anisotropy_axis{i}"
        )
    st.number_input("J", min_value=-1e3, max_value=1e3, value=0.0, key="J")

    st.markdown("### External field")
    st.radio("H axis", options=["x", "y", "z"], key="H_axis")
    st.number_input(
        "Hmin", min_value=-1000e3, max_value=1000e3, value=-400e3, key="Hmin"
    )
    st.number_input("Hmax", min_value=0.0, max_value=1000e3, value=400e3, key="Hmax")

global placeholder


def simulate():
    with st.spinner("Simulating..."):
        spec, freqs, _, Hscan = get_pimm_data(
            Ms1=st.session_state.Ms0,
            Ms2=st.session_state.Ms1,
            K1=st.session_state.K0,
            K2=st.session_state.K1,
            alpha1=st.session_state.alpha0,
            alpha2=st.session_state.alpha1,
            thickness1=st.session_state.thickness0,
            thickness2=st.session_state.thickness1,
            width1=st.session_state.width0,
            width2=st.session_state.width1,
            length1=st.session_state.length0,
            length2=st.session_state.length1,
            anisotropy_axis1=st.session_state.anisotropy_axis0,
            anisotropy_axis2=st.session_state.anisotropy_axis1,
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin,
            Hmax=st.session_state.Hmax,
            J=st.session_state.J,
        )
    with plt.style.context(["science", "nature"]):
        fig, ax = plt.subplots(dpi=300)
        ax.pcolormesh(
            Hscan / 1e3,
            freqs / 1e9,
            10 * np.log10(spec.T),
            shading="auto",
            cmap="inferno",
            rasterized=True,
        )
        ax.set_xlabel("H (kA/m)")
        ax.set_ylabel("Frequency (GHz)")
        ax.set_title("Resonance spectrum")

        st.pyplot(fig)


st.button("Simulate", on_click=simulate)
