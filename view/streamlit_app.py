# Authors: LemurPwned
from functools import partial

import streamlit as st
from autofit import autofit
from helpers import simulate_pimm, simulate_vsd

apptitle = "CMTJ simulator"

st.set_page_config(page_title=apptitle, page_icon=":rocket:")
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
            f"Ms ({i+1}) (T)",
            min_value=0.2,
            max_value=2.5,
            value=0.52,
            step=0.01,
            key=f"Ms{i}",
        )
        st.number_input(
            f"K ({i+1}) (kJ/m^3)",
            min_value=0.1,
            max_value=300e3,
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
            min_value=-5000.0,
            max_value=5000.0,
            value=0.0,
            key=f"J{j}",
            format="%.3f",
        )
    st.markdown("-----\n")
    st.markdown("## Control parameters")
    st.markdown("### External field")
    st.selectbox("H axis", options=["x", "y", "z"], key="H_axis", index=0)
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
    The only optimised values are: Ms, and K.
    All other parameters are treated as constants.
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
    placeholder = st.empty()
    st.button(
        "Run optimization",
        on_click=partial(autofit, placeholder=placeholder),
        key="opt_btn",
    )

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
        "Excitation (kA/m)", min_value=0.5, max_value=50.0, value=5.0, key="Hoex_mag"
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
