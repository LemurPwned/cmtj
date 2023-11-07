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
    ## Simulation info
    This app simulates the resonance characteristics of a MTJ device.
    The device is composed of two layers, each with its own magnetic properties.
    The number in bracket indicates the layer number.

    ## Data Upload
    If you want to upload data, to overlay it on the plot, please upload
    a file with two columns: H and f. Put H in (A/m) and f in (Hz).
    They will be rescaled to (kA/m) and (GHz) automatically.
    """
)

display_format = "%.3f"

with st.sidebar:
    st.markdown("## Simulation parameters")
    st.markdown("### Layer parameters")
    for i in range(N):
        st.slider(
            f"Ms ({i}) (T)",
            min_value=0.2,
            max_value=2.0,
            value=1.0,
            step=0.1,
            key=f"Ms{i}",
        )
        st.number_input(
            f"K ({i}) (kJ/m^3)",
            min_value=0.1,
            max_value=10e3,
            value=50.0,
            step=10.0,
            key=f"K{i}",
        )
        st.number_input(
            f"alpha ({i})",
            min_value=1e-3,
            max_value=0.1,
            value=1e-3,
            key=f"alpha{i}",
            format="%.3f",
        )
        st.number_input(
            f"thickness ({i}) (nm)",
            min_value=1.0,
            max_value=10.0,
            value=1.0,
            key=f"thickness{i}",
        )
        st.number_input(
            f"width ({i}) (um)",
            min_value=1.0,
            max_value=500.0,
            value=10.0,
            key=f"width{i}",
        )
        st.number_input(
            f"length ({i}) (um)",
            min_value=1.0,
            max_value=500.0,
            value=10.0,
            key=f"length{i}",
        )
        st.radio(
            f"anisotropy axis ({i})", options=["x", "y", "z"], key=f"anisotropy_axis{i}"
        )
    st.number_input(
        "J (mJ/m^2)", min_value=-1.0, max_value=1.0, value=0.0, key="J", format="%.2f"
    )

    st.markdown("### External field")
    st.radio("H axis", options=["x", "y", "z"], key="H_axis")
    st.number_input(
        "Hmin (kA/m)", min_value=-1000.0, max_value=1000.0, value=-400.0, key="Hmin"
    )
    st.number_input(
        "Hmax (kA/m)", min_value=0.0, max_value=1000.0, value=400.0, key="Hmax"
    )
    st.number_input(
        "int_step",
        min_value=1e-14,
        max_value=1e-12,
        value=1e-12,
        key="int_step",
        format="%.1e",
    )

global fig, ax
fig = None
ax = None


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


def simulate():
    with st.spinner("Simulating..."):
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
            J=st.session_state.J * 1e-3,
            int_step=st.session_state.int_step,
        )
    with plt.style.context(["dark_background"]):
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

        try:
            fields, freqs = read_data()
            ax.plot(fields / 1e3, freqs / 1e9, "o", color="white", label="user data")
        except (ValueError, AttributeError):
            ...
        st.pyplot(fig)


def overlay_fig():
    try:
        fields, freqs = read_data()
        if fig is not None:
            ax.plot(fields, freqs, "o", color="white")
    except (ValueError, AttributeError):
        st.error(
            "Invalid file format. Must be `\t` separated values and have H and f headers."
        )


st.button("Simulate", on_click=simulate)
st.file_uploader(
    "Upload your data here",
    help="Upload your data here. Must be `\t` separated values and have H and f headers.",
    type=["txt", "dat"],
    accept_multiple_files=False,
    key="upload",
    on_change=overlay_fig,
)
