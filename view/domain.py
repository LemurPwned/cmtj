from threading import RLock

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from helpers import read_mh_data
from simulation_fns import create_single_domain, get_axis_angles
from utils import GENERIC_BOUNDS, GENERIC_UNITS, get_init_kval

from cmtj import AxialDriver, Junction
from cmtj.utils import FieldScan

_lock = RLock()


with st.container():
    st.write("# Domain fitting")
    with st.expander("Explanation"):
        st.write(
            "#### Fit M(H) to multidomain model. \n"
            "##### Experimental data\n"
            "Upload your file with experimental data: with columns H, mx, my, mz.\n"
            "First column is always H in A/m, the rest are the components of the magnetisation in range (-1, 1)."
            "If you upload just two columns, the script will assume that the data is in the form H, mx.\n"
            "##### Simulation\n"
            "Add new domains in the sidebar panel. Area of each domain is a weight. The resulting $m_\mathrm{mixed}(h)$ plot is produces with\n\n"
            r"$m_\mathrm{mixed} = (1/{\sum_i a_i})\sum_i a_i m_i$"
            "\n\n where $a_i$ is the area of $i$th domain and $m_i$ is $m(h)$ for that domain.\n"
            "##### Coordinate system\n"
            r"$\theta$ is the polar angle and $\phi$ is the azimuthal angle."
            r" $\theta=90^\circ$ denotes fully in-plane magnetisation, "
            r"$\theta=0^\circ$ denotes out-of-plane magnetisation (positive z-direction)"
            r"$\phi=0^\circ$ denotes magnetisation pointing in the positive x direction."
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
    Mmixed = np.zeros((len(Hscan), 3), dtype=np.float16)
    Mdomains = np.zeros((st.session_state.N, len(Hscan), 3), np.float16)
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
            Mdomains[l, i] = np.array([mx, my, mz])
            Mmixed[i] += Mdomains[l, i] * weights[l]
        # for each layer take last N values
        progress_bar.progress(int((i / maxlen) * 100))
    st.session_state["Mmixed"] = Mmixed
    st.session_state["Hscan"] = Hscan
    st.session_state["Mdomains"] = Mdomains


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
    st.checkbox("Show domains in polar angles", key="domain_use_angle", value=False)
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
                help="Uniaxial anisotropy constant of the domain. Does not include shape anisotropy.",
            )

            st.number_input(
                f"area ({i+1}) " r"($\mathrm{nm^2}$)",
                min_value=1.0,
                max_value=500.0,
                value=10.0,
                key=f"area{i}",
                help="Area of the domain. In fact, this is a weight of the domain.",
            )

            st.number_input(
                r"$\theta_\mathrm{K}$ (deg.)",
                min_value=0.0,
                max_value=180.0,
                value=get_init_kval(i),
                key=f"theta_K{i}",
                help="Polar angle of the anisotropy axis",
            )
            st.number_input(
                r"$\phi_\mathrm{K}$ (deg.)",
                min_value=0.0,
                max_value=180.0,
                value=get_init_kval(i),
                key=f"phi_K{i}",
                help="Azimuthal angle of the anisotropy axis",
            )
            st.markdown("-----\n")

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
            help="Saturation magnetisation of the system. Affects the shape anisotropy",
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


def render(Hscan, Mmixed, Mdomains=None):
    if len(Hscan) <= 0 or len(Mmixed) <= 0:
        return
    with _lock:
        with plt.style.context(["dark_background"]):
            w, h = plt.figaspect(1 / 3)
            fig, ax = plt.subplots(
                1, 3, dpi=400, figsize=(w, h), sharex=True, sharey=True
            )
            if experimental_file is not None:
                fields, mh = read_mh_data()
                render_from_exp(ax, fields=fields, mh=mh)
            render_from_sim(ax, Hscan, Mmixed)
            fig.suptitle("Mixed Domain model")
            st.pyplot(fig)
            n = 3
            if st.session_state["domain_use_angle"]:
                n = 2
            w, h = plt.figaspect(st.session_state.N / n)
            fig, ax = plt.subplots(
                st.session_state.N, n, dpi=400, figsize=(w, h), sharex=True, sharey=True
            )
            if st.session_state.N <= 1:
                ax = np.array([ax])
            try:
                render_individual_domains(fig, ax, Hscan, Mdomains)
            except IndexError:
                pass
            st.pyplot(fig)


def render_from_exp(ax, fields, mh):
    if len(fields) <= 0 or len(mh) <= 0:
        st.toast("No data to render")
    fields = np.asarray(fields)
    m = np.asarray(mh)
    for i, l in zip(range(m.shape[1]), "xyz"):
        ax[i].plot(
            fields / 1e3,  # A/m --> kA/m
            m[:, i],
            label=f"$m_{l}$",
            color="yellow",
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


def render_individual_domains(fig, ax, Hscan, M: list[list[float]]):
    n = 3
    if st.session_state["domain_use_angle"]:
        n = 2
        ax[0, 0].set_title(r"$\theta$")
        ax[0, 1].set_title(r"$\phi$")
        for i, domain_m in enumerate(M):
            # convert to polar
            theta_m, phi_m, _ = FieldScan.vector2angle(
                domain_m[:, 0], domain_m[:, 1], domain_m[:, 2]
            )
            ax[i, 0].plot(
                Hscan / 1e3, theta_m, label=rf"$\theta_{i+1}$", color="crimson"
            )
            ax[i, 1].plot(
                Hscan / 1e3, phi_m, label=rf"$\phi_{i+1}$", color="forestgreen"
            )
            ax[i, 0].set_ylabel(rf"$\theta_{i+1}$ (deg.)")
            ax[i, 1].set_ylabel(rf"$\phi_{i+1}$ (deg.)")
    else:
        n = 3
        ax[0, 0].set_title(r"$m_\mathrm{x}$")
        ax[0, 1].set_title(r"$m_\mathrm{y}$")
        ax[0, 2].set_title(r"$m_\mathrm{z}$")
        for i, domain_m in enumerate(M):
            ax[i, 0].plot(
                Hscan / 1e3, domain_m[:, 0], label=rf"$m_{i+1}$", color="crimson"
            )
            ax[i, 1].plot(
                Hscan / 1e3, domain_m[:, 1], label=rf"$m_{i+1}$", color="forestgreen"
            )
            ax[i, 2].plot(
                Hscan / 1e3, domain_m[:, 2], label=rf"$m_{i+1}$", color="royalblue"
            )
            ax[i, 0].set_ylabel(rf"$m_{i+1}$")
    for j in range(n):
        ax[-1, j].set_xlabel("H (kA/m)")
    fig.suptitle("Individual Domains")  # fig.legend()


simulate_btn = st.button(
    "Simulate",
    key="simulate",
    on_click=simulate_domains,
    type="primary",
)

render(
    st.session_state.get("Hscan", []),
    st.session_state.get("Mmixed", []),
    Mdomains=st.session_state.get("Mdomains", []),
)
