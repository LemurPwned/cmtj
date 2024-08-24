from threading import RLock

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from helpers import read_mh_data
from simulation_fns import create_single_domain, get_axis_angles
from utils import GENERIC_BOUNDS, GENERIC_UNITS

from cmtj import AxialDriver, Junction
from cmtj.utils import FieldScan

_lock = RLock()

with st.container():
    st.write("# Domain fitting")

    st.write(
        "Fit M(H) to multidomain model. "
        "Upload your file with experimental data: with columns H, mx, my, mz.\n"
        "First column is always H in A/m, the rest are the components of the magnetisation in range (-1, 1)."
        "If you upload just two columns, the script will assume that the data is in the form H, mx."
    )

    progress_bar = st.progress(0)


def simulate_domains():
    htheta, hphi = get_axis_angles(st.session_state.H_axis)
    Hscan, Hvecs = FieldScan.amplitude_scan(
        st.session_state.Hmin * 1e3,
        st.session_state.Hmax * 1e3,
        st.session_state.Hsteps,
        htheta,
        hphi,
        back=st.session_state.back_to_back_scan,
    )

    domains = []
    for i in range(st.session_state.N):
        dom = create_single_domain(i)
        domains.append(dom)
    j = Junction(layers=domains)
    maxlen = len(Hvecs)
    # plot data
    weights = [st.session_state[f"area{i}"] * 1e-18 for i in range(st.session_state.N)]
    wsum = sum(weights)
    weights = [w / wsum for w in weights]
    Mmixed = np.zeros((len(Hscan), 3))
    for i, H in enumerate(Hvecs):
        j.clearLog()
        j.setLayerExternalFieldDriver("all", AxialDriver(*H))
        j.runSimulation(
            st.session_state.sim_time * 1e-9,
            st.session_state.int_step,
            st.session_state.int_step,
        )
        log = j.getLog()
        for l in range(st.session_state.N):
            mx = np.mean(log[f"domain_{l}_mx"][-st.session_state.last_N :])
            my = np.mean(log[f"domain_{l}_my"][-st.session_state.last_N :])
            mz = np.mean(log[f"domain_{l}_mz"][-st.session_state.last_N :])
            Mmixed[i] += np.array([mx, my, mz]) * weights[l]
        # for each layer take last N values
        progress_bar.progress(int((i / maxlen) * 100))
    st.session_state["Mmixed"] = Mmixed
    st.session_state["Hscan"] = Hscan


with st.sidebar:
    experimental_file = st.file_uploader(
        "Upload your experimental data here",
        help="Upload your data here. "
        "Must be `\t` separated values and have H and m headers. "
        "Express H in A/m and m in units",
        type=["txt", "dat", "tsv"],
        accept_multiple_files=False,
        key="upload",
    )
    with st.expander("Simulation parameters", expanded=True):
        st.selectbox(
            "H axis", options=["x", "y", "z", "xy", "xz", "yz"], key="H_axis", index=0
        )
        st.number_input(
            "Hmin (kA/m)", min_value=-1000.0, max_value=1000.0, value=-400.0, key="Hmin"
        )
        st.number_input(
            "Hmax (kA/m)", min_value=0.0, max_value=1000.0, value=400.0, key="Hmax"
        )
        st.number_input(
            "H steps", min_value=1, max_value=1000, value=50, key="Hsteps", format="%d"
        )
        # boolean if to back to back scan
        st.checkbox("Back to back scan", key="back_to_back_scan", value=True)
        st.number_input(
            "int_step",
            min_value=1e-14,
            max_value=1e-12,
            value=1e-12,
            key="int_step",
            format="%.1e",
        )
        st.number_input(
            "sim_time (ns)",
            min_value=1,
            max_value=500,
            value=16,
            key="sim_time",
            format="%d",
        )
        st.number_input(
            "last N values to consider",
            min_value=1,
            max_value=1000,
            value=100,
            key="last_N",
        )

    shared_params = st.expander("Shared parameters", expanded=True)
    with shared_params:
        st.slider(
            f"Ms ({GENERIC_UNITS['Ms']})",
            min_value=GENERIC_BOUNDS["Ms"][0],
            max_value=GENERIC_BOUNDS["Ms"][1],
            value=1.8,
            step=0.01,
            key="Ms_shared",
        )
        st.number_input(
            "thickness (nm)",
            min_value=1.0,
            max_value=10.0,
            value=1.0,
            key="thickness_shared",
        )
        st.number_input(
            "alpha",
            min_value=1e-3,
            max_value=0.1,
            value=1e-3,
            key="alpha_shared",
            format="%.3f",
        )
    domain_params = st.expander("Domain parameters", expanded=True)
    with domain_params:
        N = st.number_input(
            "Number of domains", min_value=1, max_value=10, value=1, key="N"
        )
        for i in range(N):
            st.markdown(f"#### Domain {i+1}")
            unit_k = GENERIC_UNITS["K"]
            st.number_input(
                f"K ({i+1}) " r" ($\mathrm{" f"{unit_k}" "}$)",
                min_value=GENERIC_BOUNDS["K"][0],
                max_value=GENERIC_BOUNDS["K"][1],
                value=1.2e3,
                step=10.0,
                key=f"K{i}",
            )

            st.number_input(
                f"area ({i+1}) " r"($\mathrm{nm^2}$)",
                min_value=1.0,
                max_value=500.0,
                value=10.0,
                key=f"area{i}",
            )
            st.radio(
                f"anisotropy axis ({i+1})",
                options=["x", "y", "z"],
                key=f"anisotropy_axis{i}",
                index=2,
            )
            st.markdown("-----\n")


def render(Hscan, Mmixed):
    if len(Hscan) <= 0 or len(Mmixed) <= 0:
        return
    with _lock:
        with plt.style.context(["science", "nature"]):
            w, h = plt.figaspect(1 / 3)
            fig, ax = plt.subplots(
                1, 3, dpi=400, figsize=(w, h), sharex=True, sharey=True
            )
            if experimental_file is not None:
                fields, mh = read_mh_data()
                render_from_exp(ax, fields=fields, mh=mh)
            render_from_sim(ax, Hscan, Mmixed)

            st.pyplot(fig)


def render_from_exp(ax, fields, mh):
    if len(fields) <= 0 or len(mh) <= 0:
        st.toast("No data to render")

    m = np.asarray(mh)
    for i, l in zip(range(m.shape[1]), "xyz"):
        ax[i].plot(
            fields,
            m[:, i],
            label=f"$m_{l}$",
            color="black",
            marker="x",
            alpha=0.5,
            linestyle="--",
        )


def render_from_sim(ax, Hscan, Mmixed):
    ax[0].plot(Hscan / 1e3, Mmixed[:, 0], label="$m_x$", color="crimson")
    ax[1].plot(Hscan / 1e3, Mmixed[:, 1], label="$m_y$", color="forestgreen")
    ax[2].plot(Hscan / 1e3, Mmixed[:, 2], label="$m_z$", color="royalblue")

    ax[0].set_title(r"$m_x$")
    ax[1].set_title(r"$m_y$")
    ax[2].set_title(r"$m_z$")
    ax[0].set_xlabel("H (kA/m)")
    ax[1].set_xlabel("H (kA/m)")
    ax[2].set_xlabel("H (kA/m)")


simulate_btn = st.button(
    "Simulate",
    key="simulate",
    on_click=simulate_domains,
    type="primary",
)

render(st.session_state.get("Hscan", []), st.session_state.get("Mmixed", []))
