# Authors: LemurPwned
import streamlit as st
from helpers import simulate_pimm, simulate_vsd

N = 2
apptitle = "CMTJ simulator"

st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")
st.title(apptitle)
container = st.container()
container.markdown(
    """
    ## Data Upload
    If you want to upload data, to overlay it on the plot, please upload
    a file with two columns: H and f. Put H in (A/m) and f in (Hz).
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

display_format = "%.3f"

with st.sidebar:
    st.markdown("## Simulation parameters")
    st.markdown("### Layer parameters")
    for i in range(N):
        st.slider(
            f"Ms ({i}) (T)",
            min_value=0.2,
            max_value=2.0,
            value=1.2,
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
            f"anisotropy axis ({i})",
            options=["x", "y", "z"],
            key=f"anisotropy_axis{i}",
        )
    st.number_input(
        "J (mJ/m^2)",
        min_value=-1.0,
        max_value=1.0,
        value=0.0,
        key="J",
        format="%.2f",
    )

    st.markdown("### External field")
    st.radio("H axis", options=["x", "y", "z"], key="H_axis", index=2)
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


pimm_tab, vsd_tab = st.tabs(["PIMM", "VSD"])
with vsd_tab:
    st.number_input(
        "Frequency min (GHz)", min_value=0.0, max_value=50.0, value=0.0, key="fmin"
    )
    st.number_input(
        "Frequency max (GHz)", min_value=0.0, max_value=50.0, value=30.0, key="fmax"
    )
    st.number_input(
        "Number of frequencies",
        min_value=1,
        max_value=100,
        value=30,
        key="nf",
        format="%d",
    )

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
