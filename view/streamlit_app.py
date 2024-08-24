# Authors: LemurPwned
from functools import partial

import streamlit as st
from autofit import autofit
from helpers import simulate_pimm, simulate_vsd
from utils import GENERIC_BOUNDS, GENERIC_UNITS

apptitle = "CMTJ simulator"

st.title(apptitle)

container = st.container()
N = container.number_input(
    "Number of layers", min_value=1, max_value=10, value=1, key="N", format="%d"
)
container.markdown(
    """
    ## Data Upload
    If you want to upload data, to overlay it on the plot, please upload
    a file with \t separated columns: H, f1, f2, ... Put H in (A/m) and f in (Hz).
    They will be rescaled to (kA/m) and (GHz) automatically.
    """
)
container.file_uploader(
    "Upload your data here",
    help="Upload your data here. Must be `\t` separated values and have H and f headers.",
    type=["txt", "dat"],
    accept_multiple_files=False,
    key="upload",
)

with st.sidebar:
    st.markdown("## Simulation parameters")
    st.markdown("### Layer parameters")
    for i in range(N):
        st.markdown(f"#### Layer {i+1}")
        st.slider(
            f"Ms ({i+1}) ({GENERIC_UNITS['Ms']})",
            min_value=GENERIC_BOUNDS["Ms"][0],
            max_value=GENERIC_BOUNDS["Ms"][1],
            value=0.52,
            step=0.01,
            key=f"Ms{i}",
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
        )
        st.number_input(
            f"length ({i+1}) (um)",
            min_value=1.0,
            max_value=500.0,
            value=10.0,
            key=f"length{i}",
        )
        st.radio(
            f"anisotropy axis ({i+1})",
            options=["x", "y", "z"],
            key=f"anisotropy_axis{i}",
            index=2,
        )
        st.markdown("-----\n")

    st.markdown("### Interlayer parameters")
    for j in range(N - 1):
        st.number_input(
            f"J ({j+1}<-->{j+2}) (uJ/m^2)",
            min_value=GENERIC_BOUNDS["J"][0],
            max_value=GENERIC_BOUNDS["J"][1],
            value=0.0,
            key=f"J{j}",
            format="%.3f",
        )
    st.markdown("-----\n")
    st.markdown("## Control parameters")
    st.markdown("### External field")
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
        "max_freq (GHz)",
        min_value=1,
        max_value=100,
        value=50,
        key="max_freq",
        format="%d",
    )


pimm_tab, vsd_tab, opt_tab = st.tabs(["PIMM", "VSD", "Optimization"])
with opt_tab:
    st.markdown(
        """
    ## Optimization

    Run Bayesian optimisation -- fitting data is source from file upload.
    Select the number of iterations.
    The only optimised values are: Ms, and K and J.
    All other parameters are treated as constants.
    Narrow the bounds to improve the optimisation.

    ### Fixed parameters
    Parameters marked as fixed will not be optimised.
    In the panel below, __the lower bound value will be fixed__

    ### Values to set on the left hand side panel
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
    )
    st.selectbox("excitation axis", options=["x", "y", "z"], key="Hoex", index=1)
    st.markdown(
        """### Simulation info

    This is Voltage Spin Diode experiment (WIP).
    """
    )

    st.button("Simulate VSD", on_click=simulate_vsd, key="VSD_btn")

with pimm_tab:
    fn = simulate_pimm
    st.markdown(
        """### Simulation info

    This app simulates PIMM experiment.
    """
    )
    st.button("Simulate PIMM", on_click=simulate_pimm, key="PIMM_btn")
    st.number_input(
        "Hoe (kA/m)", min_value=0.05, max_value=50.0, value=0.05, key="Hoe_mag"
    )
    st.radio(
        "Hoe axis",
        options=["x", "y", "z"],
        key="Hoeaxis",
        index=1,
    )
