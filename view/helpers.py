import io
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from datetime import datetime
from simulation_fns import get_pimm_data, get_vsd_data

from cmtj.utils import Filters


def read_mh_data():
    filedata = st.session_state.upload.read().decode("utf-8")
    lines = filedata.split("\n")
    fields, mag = [], []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        vals = line.split()
        if len(vals):
            fields.append(float(vals[0]))
            mag.append([float(val) for val in vals[1:]])
    return fields, mag


def read_data():
    filedata = st.session_state.upload.read().decode("utf-8")
    lines = filedata.split("\n")
    fields, freqs = [], []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        vals = line.split()
        if len(vals):
            fields.append(float(vals[0]))
            freqs.append([float(val) for val in vals[1:]])
    return fields, freqs


def plot_optim(hvals, target, simulation, title="Optimisation"):
    hvals = np.asarray(hvals)
    with plt.style.context(["dark_background"]):
        fig, ax = plt.subplots(dpi=300)
        ax.plot(
            hvals / 1e3, np.asarray(target) / 1e9, "o", color="crimson", label="Target"
        )
        ax.plot(
            simulation["Hmag"], simulation["frequency"], "x", color="white", label="Sim"
        )
        ax.set_xlabel("H (kA/m)")
        ax.set_ylabel("Frequency (GHz)")
        ax.legend()
        fig.suptitle(title)
        st.pyplot(fig)


def plot_data(
    Hscan, freqs, spec, mag=None, title: str = "Resonance spectrum", mode: str = "none", return_fig: bool = False
):
    Hscan_spec = Hscan
    if st.session_state["Hreturn"] and mode.lower() == "pimm":
        print("Adding spectrum")
        Hlen = len(Hscan) // 2
        Hscan_spec = Hscan[:Hlen]
        spec_left_right = spec[:Hlen, :]
        spec_right_left = spec[Hlen:, :]
        # invert the right side
        spec_right_left = spec_right_left[::-1, :]
        spec = spec_left_right + spec_right_left

    with plt.style.context(["dark_background"]):
        if mag is not None:
            fig = plt.figure(dpi=300)
            gs = plt.GridSpec(3, 6)
            indx = 4
            ax1 = fig.add_subplot(gs[:, :indx])
            prev_ax = []
            for i, lab in enumerate("xyz"):
                _ax = fig.add_subplot(gs[i, indx:])
                prev_ax.append(_ax)
                _ax.plot(Hscan / 1e3, mag[:, i], color="white")
                # move ticks to the right side
                _ax.yaxis.tick_right()
                _ax.yaxis.set_label_position("right")
                _ax.set_ylim(-1.1, 1.1)
                if lab != "z":
                    # remove x ticks
                    _ax.set_xticks([])
                _ax.set_ylabel(rf"$m_{lab}$", rotation=270, labelpad=15)
            prev_ax[-1].set_xlabel("H (kA/m)")
            fig.subplots_adjust(wspace=0.05, hspace=0.05)
            fig.align_ylabels()
        else:
            fig, ax1 = plt.subplots(dpi=300)
        ax1.pcolormesh(
            Hscan_spec / 1e3,
            freqs / 1e9,
            10 * np.log10(np.abs(spec.T)),
            shading="auto",
            cmap="inferno",
            rasterized=True,
        )
        ax1.set_xlabel("H (kA/m)")
        ax1.set_ylabel("Frequency (GHz)")
        fig.suptitle(title)
        try:
            fields, freqs = read_data()
            freqs = np.asarray(freqs)
            fields = np.asarray(fields)
            cmap = plt.colormaps["Pastel2"]
            for f_indx in range(freqs.shape[1]):
                ax1.plot(
                    fields / 1e3,
                    freqs[:, f_indx] / 1e9,
                    "o",
                    markeredgecolor="white",
                    color=cmap(f_indx),
                    label=f"Data {f_indx}",
                )
        except (ValueError, AttributeError) as e:
            print(f"Error plotting data: {e}")
        st.pyplot(fig)
        if return_fig:
            return fig
        return None


def simulate_vsd():
    with st.spinner("Simulating VSD..."):
        spec, freqs, Hscan = get_vsd_data(
            int_step=st.session_state.int_step,
            fmin=st.session_state.fmin * 1e9,
            fmax=st.session_state.fmax * 1e9,
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin * 1e3,
            Hmax=st.session_state.Hmax * 1e3,
            Hsteps=st.session_state.Hsteps,
            Hoex=st.session_state.Hoex,
            fstep=st.session_state.fstep * 1e9,
            Rtype=st.session_state.res_type,
            sim_time=st.session_state.sim_time * 1e-9,
            Hoex_mag=st.session_state.Hoex_mag,
        )
        spec = Filters.detrend_axis(spec, axis=0)
    plot_data(Hscan, freqs, spec, title="VSD spectrum", mode="vsd")


def fig_to_bytes(fig):
    """Convert matplotlib figure to PNG bytes."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()


def add_to_history(fig, max_history: int):
    """Add figure to history in session state."""
    if "pimm_history" not in st.session_state:
        st.session_state.pimm_history = []

    # Convert figure to bytes
    fig_bytes = fig_to_bytes(fig)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Add to history (most recent first)
    history_entry = {
        "timestamp": timestamp,
        "image_bytes": fig_bytes,
    }
    st.session_state.pimm_history.insert(0, history_entry)

    # Keep only max_history entries
    if len(st.session_state.pimm_history) > max_history:
        st.session_state.pimm_history = st.session_state.pimm_history[:max_history]


def simulate_pimm():
    with st.spinner("Simulating PIMM..."):
        spec, freqs, output, Hscan = get_pimm_data(
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin * 1e3,
            Hmax=st.session_state.Hmax * 1e3,
            Hsteps=st.session_state.Hsteps,
            int_step=st.session_state.int_step,
            sim_time=st.session_state.sim_time * 1e-9,
        )
    mag = np.asarray(output["m_avg"])
    fig = plot_data(Hscan, freqs, spec, mag=mag, title="PIMM spectrum", mode="pimm", return_fig=True)

    # Save to history before closing
    if fig is not None:
        max_history = st.session_state.get("pimm_max_history", 5)
        add_to_history(fig, max_history)
        plt.close(fig)  # Close to free memory
