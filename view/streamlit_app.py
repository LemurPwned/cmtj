# Authors: LemurPwned
from functools import partial

import streamlit as st
from autofit import autofit
from helpers import simulate_pimm, simulate_vsd
from utils import GENERIC_BOUNDS, GENERIC_UNITS
import json


def export_session_state():
    export_dict = {}
    opts = ["_btn", "_file", "_state", "low_", "up_", "check_", "upload"]
    for k, v in st.session_state.items():
        skip = any(forb_opts in k for forb_opts in opts)
        # Skip pimm_history as it contains bytes that aren't JSON serializable
        if k == "pimm_history":
            skip = True
        if not skip:
            export_dict[k] = v

    return json.dumps(export_dict)


def import_session_state(file):
    try:
        data = json.load(file)
        for k, v in data.items():
            st.session_state[k] = v
    except json.JSONDecodeError:
        st.error("Error: Invalid JSON file format. Please upload a valid JSON file.")


with st.expander("# Read me"):
    st.write(
        """
        ### Data Upload
        If you want to upload data, to overlay it on the plot, please upload
        a file with \t separated columns: H, f1, f2, ... Put H in (A/m) and f in (Hz).
        They will be rescaled to (kA/m) and (GHz) automatically.
        """
    )

with st.sidebar:
    with st.expander("Export/Import"):
        st.download_button(
            label="Export session state",
            data=export_session_state(),
            file_name="session_state.json",
            mime="application/json",
            type="primary",
            help="Export the current session state to a JSON file. "
            "You can use this to save your current settings and load them later or "
            "share them with others.",
        )

        st.file_uploader(
            "Upload session state",
            help="Upload your data here. Must be `\t` separated values and have H and f headers.",
            type=["json"],
            accept_multiple_files=False,
            key="import_file",
        )
        if st.session_state.import_file:
            import_session_state(st.session_state.import_file)

        st.file_uploader(
            "Upload your data here",
            help="Upload your data here. Must be `\t` separated values and have H and f headers.",
            type=["txt", "dat"],
            accept_multiple_files=False,
            key="upload",
        )
    N = st.number_input(
        "Number of layers", min_value=1, max_value=10, value=1, key="N", format="%d"
    )
    with st.expander("Layer parameters", expanded=True):
        for i in range(N):
            st.markdown(f"#### Layer {i+1}")
            st.slider(
                f"Ms ({i+1}) ({GENERIC_UNITS['Ms']})",
                min_value=GENERIC_BOUNDS["Ms"][0],
                max_value=GENERIC_BOUNDS["Ms"][1],
                value=0.52,
                step=0.01,
                key=f"Ms{i}",
                help="Saturation magnetisation (T). Affectes shape aniostropy.",
            )
            st.number_input(
                f"K ({i+1}) ({GENERIC_UNITS['K']})",
                min_value=GENERIC_BOUNDS["K"][0],
                max_value=GENERIC_BOUNDS["K"][1],
                value=150.0,
                step=10.0,
                key=f"K{i}",
            )
            st.number_input(
                f"alpha ({i+1})",
                min_value=1e-3,
                max_value=0.1,
                value=1e-3,
                key=f"alpha{i}",
                format="%.3f",
                help="Gilbert damping parameter",
            )
            st.number_input(
                f"thickness ({i+1}) (nm)",
                min_value=1.0,
                max_value=10.0,
                value=1.0,
                key=f"thickness{i}",
            )
            st.number_input(
                f"width ({i+1}) (um)",
                min_value=1.0,
                max_value=500.0,
                value=10.0,
                key=f"width{i}",
                help="Width of the sample. For resistance calculations",
            )
            st.number_input(
                f"length ({i+1}) (um)",
                min_value=1.0,
                max_value=500.0,
                value=10.0,
                key=f"length{i}",
                help="Length of the sample. For resistance calculations",
            )
            st.number_input(
                r"$\theta_\mathrm{K}$ (deg.)",
                min_value=0.0,
                max_value=180.0,
                value=0.0,
                key=f"theta_K{i}",
                help="Polar angle of the anisotropy axis",
            )
            st.number_input(
                r"$\phi_\mathrm{K}$ (deg.)",
                min_value=0.0,
                max_value=180.0,
                value=0.0,
                key=f"phi_K{i}",
                help="Azimuthal angle of the anisotropy axis",
            )
            st.write("Demagnetization field")
            st.number_input(
                f"Nxx ({i+1})",
                value=0.0,
                key=f"Nxx{i}",
                format="%0.5f",
                help="Demagnetization field component diagonal(xx)",
            )
            st.number_input(
                f"Nyy ({i+1})",
                value=0.0,
                key=f"Nyy{i}",
                format="%0.5f",
                help="Demagnetization field component diagonal(yy)",
            )
            st.number_input(
                f"Nzz ({i+1})",
                value=1.0,
                key=f"Nzz{i}",
                format="%0.5f",
                help="Demagnetization field component diagonal(zz)",
            )
            st.markdown("-----\n")

    with st.expander("Interlayer parameters"):
        for j in range(N - 1):
            st.number_input(
                f"J ({j+1}<-->{j+2}) (uJ/m^2)",
                min_value=GENERIC_BOUNDS["J"][0],
                max_value=GENERIC_BOUNDS["J"][1],
                value=0.0,
                key=f"J{j}",
                format="%.3f",
                help="Interlayer exchange coupling constant",
            )

            st.number_input(
                f"$J_2$ ({j+1}<-->{j+2}) (uJ/m^2)",
                min_value=GENERIC_BOUNDS["J"][0],
                max_value=GENERIC_BOUNDS["J"][1],
                value=0.0,
                key=f"J2{j}",
                format="%.4f",
                help="Biquadratic interlayer exchange coupling constant",
            )

            st.number_input(
                f"ilD ({j+1}<-->{j+2}) (uJ/m^2)",
                min_value=GENERIC_BOUNDS["ilD"][0],
                max_value=GENERIC_BOUNDS["ilD"][1],
                value=0.0,
                key=f"ilD{j}",
                format="%.4f",
                help="Interlayer DMI constant (D vector along z axis)",
            )

    with st.expander("Simulation & control parameters"):
        st.selectbox(
            "H axis", options=["x", "y", "z", "xy", "xz", "yz"], key="H_axis", index=0
        )
        st.number_input(
            "Hmin (kA/m)", min_value=-1e6, max_value=1e6, value=-400.0, key="Hmin"
        )
        st.number_input(
            "Hmax (kA/m)", min_value=-1e6, max_value=1e6, value=400.0, key="Hmax"
        )
        st.number_input(
            "H steps", min_value=1, max_value=1000, value=50, key="Hsteps", format="%d"
        )
        st.checkbox("Back-and-forth-field (PIMM only)", value=False, key="Hreturn")
        st.number_input(
            "int_step",
            min_value=1e-14,
            max_value=1e-12,
            value=1e-12,
            key="int_step",
            format="%.1e",
            help="Integration step for the simulation",
        )
        st.number_input(
            "sim_time (ns)",
            min_value=1,
            max_value=500,
            value=16,
            key="sim_time",
            format="%d",
            help="Simulation time in nanoseconds",
        )
        st.number_input(
            "max_freq (GHz)",
            min_value=1,
            max_value=100,
            value=50,
            key="max_freq",
            format="%d",
            help="Maximum frequency (cutoff) visible in plot",
        )

pimm_tab, vsd_tab, opt_tab = st.tabs(["PIMM", "VSD", "Optimization"])
with opt_tab:
    with st.expander("Optimization parameters"):
        st.markdown(
            """
        ### Optimization

        Run Bayesian optimisation -- fitting data is source from file upload.
        Select the number of iterations.
        The only optimised values are: Ms, and K and J.
        All other parameters are treated as constants.
        Narrow the bounds to improve the optimisation.

        #### Fixed parameters
        Parameters marked as fixed will not be optimised.
        In the panel below, __the lower bound value will be fixed__

        #### Values to set on the left hand side panel
        - H axis
        - K axis for each layer
        - alpha for each layer
        - thickness for each layer
        """
        )
    st.number_input(
        "Number of iterations",
        min_value=1,
        max_value=100,
        value=10,
        key="n_iters",
        format="%d",
    )
    st.number_input(
        "Number of suggestions per iteration",
        min_value=1,
        max_value=8,
        value=4,
        key="n_suggestions",
    )
    placeholder = st.empty()
    st.button(
        "Run optimization",
        on_click=partial(autofit, placeholder=placeholder),
        key="opt_btn",
    )

    st.markdown("#### Optimisation bounds")
    grid = st.columns([2, 3, 3, 2])
    grid[0].text("Param")
    grid[1].text("Min")
    grid[2].text("Max")
    grid[3].text("Fix")
    for i in range(st.session_state.N):
        for param_name in ("Ms", "K", "J"):
            if param_name == "J" and i >= (st.session_state.N - 1):
                # we don't have a J
                continue
            grid = st.columns([2, 3, 3, 2])

            grid[0].write(f"{param_name} {i+1} ({GENERIC_UNITS[param_name]})")
            grid[1].text_input(
                f"{param_name} ({i})",
                label_visibility="collapsed",
                value=GENERIC_BOUNDS[param_name][0],
                key=f"low_{param_name}{i}",
            )
            grid[2].text_input(
                f"{param_name} ({i})",
                label_visibility="collapsed",
                value=GENERIC_BOUNDS[param_name][1],
                key=f"up_{param_name}{i}",
            )
            grid[3].toggle(f"fix {param_name} ({i+1})", key=f"check_{param_name}{i}")

with vsd_tab:
    with st.expander("VSD simulation parameters"):
        st.markdown(
            """### Simulation info

        This is Voltage Spin Diode experiment.
        The frequency is the sinusoidal excitation frequency fed into the system.
        """
        )
    st.number_input(
        "Frequency min (GHz)", min_value=0.0, max_value=50.0, value=0.0, key="fmin"
    )
    st.number_input(
        "Frequency max (GHz)", min_value=0.0, max_value=50.0, value=20.0, key="fmax"
    )
    st.number_input(
        "Frequency step (GHz)",
        min_value=0.1,
        max_value=10.0,
        value=0.5,
        key="fstep",
        step=0.1,
    )
    st.number_input(
        "Excitation (A/m)", min_value=1e-5, max_value=1e5, value=500.0, key="Hoex_mag"
    )
    st.selectbox(
        "Resistance type",
        options=["Rx", "Ry"],
        key="res_type",
        index=1,
        help="Resistance type for the VSD simulation. PRL, Kim et al. 2016",
    )
    st.selectbox("excitation axis", options=["x", "y", "z"], key="Hoex", index=1)

    st.button("Simulate VSD", on_click=simulate_vsd, key="VSD_btn")

with pimm_tab:
    fn = simulate_pimm
    with st.expander("PIMM simulation parameters"):
        st.write(
            """### Simulation info\n
        This app simulates PIMM experiment wherin we apply a step impulse.
        Here, the impulse is field impulse (in the experiment that could be Oersted field impulse).
        This impulse is called $\mathbf{H}_{oe}$.
        """
        )
    st.button("Simulate PIMM", on_click=simulate_pimm, key="PIMM_btn")
    st.number_input(
        "Hoe (kA/m)",
        min_value=-500.0,
        max_value=500.0,
        value=50.0,
        key="Hoe_mag",
        help="Magnitude of the Oersted field impulse (PIMM excitation)",
    )
    st.radio(
        "Hoe axis",
        options=["x", "y", "z"],
        key="Hoeaxis",
        index=1,
        help="Direction of the Oersted field impulse",
    )

    # History settings
    st.markdown("---")
    st.markdown("### Simulation History")
    st.number_input(
        "Max history size",
        min_value=1,
        max_value=20,
        value=5,
        key="pimm_max_history",
        format="%d",
        help="Maximum number of past simulations to keep in history",
    )

    # Display history gallery
    if "pimm_history" in st.session_state and len(st.session_state.pimm_history) > 0:
        st.markdown("#### Past Simulations")
        history = st.session_state.pimm_history

        # Display in a grid layout (3 columns for smaller images)
        n_cols = min(4, len(history))
        cols = st.columns(n_cols)

        for idx, entry in enumerate(history):
            col_idx = idx % n_cols
            with cols[col_idx]:
                st.image(entry["image_bytes"], width=600)
                st.caption(f"Sim {idx + 1} - {entry['timestamp']}")

        # Clear history button
        if st.button("Clear History", key="clear_history_btn"):
            st.session_state.pimm_history = []
            st.rerun()
    else:
        st.info("No simulation history yet. Run a PIMM simulation to see it here.")
