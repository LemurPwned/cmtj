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


def create_junction(
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
    Kdir1 = get_axis_cvector(anisotropy_axis1)
    Kdir2 = get_axis_cvector(anisotropy_axis2)
    layer1 = Layer(
        Ms=Ms1,
        Ku=K1,
        alpha=alpha1,
        thickness=thickness1,
        anis=Kdir1,
    )
    layer2 = Layer(
        Ms=Ms2,
        Ku=K2,
        alpha=alpha2,
        thickness=thickness2,
        anis=Kdir2,
    )
    j = Junction([layer1, layer2])
    j.setLayerAnisotropyDriver("top", K1)
    j.setLayerAnisotropyDriver("bottom", K2)
    j.setIECDriver("top", "bottom", J)
    htheta, hphi = get_axis_cvector(H_axis)
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

    """
)


with st.sidebar:
    st.markdown("## Simulation parameters")
    st.markdown("### Layer parameters")
    for i in range(N):
        st.slider(f"Ms {i}", min_value=0.2, max_value=2.0, value=1.0, step=0.1)
        st.number_input(f"K {i}", min_value=0.1, max_value=100e6, value=50e3)
        st.number_input(f"alpha {i}", min_value=1e-5, max_value=0.1, value=0.001)
        st.number_input(f"thickness {i}", min_value=1e-9, max_value=10e-9, value=1e-9)
        st.number_input(f"width {i}", min_value=1e-6, max_value=500e-6, value=10e-6)
        st.number_input(f"length {i}", min_value=1e-6, max_value=500e-6, value=10e-6)
        st.radio(f"anisotropy axis {i}", options=["x", "y", "z"])

    st.markdown("### External field")
    st.radio(f"H axis {i}", options=["x", "y", "z"])
    st.number_input(f"Hmin {i}", min_value=-1000e3, max_value=1000e3, value=-400e3)
    st.number_input(f"Hmax {i}", min_value=0.0, max_value=1000e3, value=400e3)


def simulate():
    st.write("Simulating...")
    st.write(st.session_state)

st.button("Simulate", on_click=simulate)
