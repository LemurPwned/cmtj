import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from simulation_fns import get_pimm_data, get_vsd_data


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


def plot_data(Hscan, freqs, spec, mag=None, title="Resonance spectrum"):
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
            Hscan / 1e3,
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
            for f_indx in range(freqs.shape[1]):
                ax1.plot(
                    fields / 1e3,
                    freqs[:, f_indx] / 1e9,
                    "o",
                    color="white",
                    label="user data",
                )
        except (ValueError, AttributeError):
            ...
        st.pyplot(fig)


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
            nf=st.session_state.nf,
        )
    plot_data(Hscan, freqs, spec, title="VSD spectrum")


def simulate_pimm():
    with st.spinner("Simulating PIMM..."):
        spec, freqs, output, Hscan = get_pimm_data(
            H_axis=st.session_state.H_axis,
            Hmin=st.session_state.Hmin * 1e3,
            Hmax=st.session_state.Hmax * 1e3,
            Hsteps=st.session_state.Hsteps,
            int_step=st.session_state.int_step,
        )
    mag = np.asarray(output["m_avg"])
    plot_data(Hscan, freqs, spec, mag=mag, title="PIMM spectrum")
