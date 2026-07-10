"""
Antiferromagnetic (AFM) Layer: Exchange-Enhanced Resonance

GOAL:
This simulation demonstrates the defining physical signature of a two-sublattice
antiferromagnetic (AFM) layer: an exchange-enhanced resonance frequency that scales
with the square root of the intra-layer Neel exchange strength, sqrt(|J_afm|). This is
the AFM analogue of ferromagnetic resonance (FMR) -- antiferromagnetic resonance (AFMR)
-- and it is what pushes AFM resonances from the GHz range (typical FM) into the
hundreds-of-GHz to THz range.

CRITICAL SIMULATION MECHANISMS:
1. Single AFM Layer: A single `Layer` with two exchange-coupled sublattices (mag/mag2),
   created via `Layer.createAFMLayer`. Both sublattices share Ms/anisotropy/damping
   (same material), coupled by a strongly negative intra-layer exchange driver.
2. Impulsive Excitation: A brief, strong transverse field pulse cants the two
   sublattices away from their antiparallel (Neel) ground state. Because the intra-layer
   exchange field is enormous compared to typical interlayer/anisotropy fields, only a
   large transient kick produces measurable canting -- this is a numerical stand-in for
   the impulsive/THz excitation used in real AFMR experiments.
3. Free Ring-Down + FFT: After the pulse, the two sublattices precess coherently
   (antiparallel-locked, out of phase) and relax back to equilibrium. FFT of the
   transverse magnetization component extracts the resonance frequency.
4. Exchange Sweep: Repeating this for a range of |J_afm| values and plotting resonance
   frequency vs sqrt(|J_afm|) reproduces the expected AFMR scaling law -- a strong linear
   trend is the sanity check that the two-sublattice coupled dynamics are physically
   correct, not just numerically stable.

The simulation generates one figure with three panels: (1)-(2) the antiparallel-locked
precession trajectory of both sublattices after the pulse, at the weak- and
strong-exchange edges of the sweep (very different ring-down timescales), and (3)
resonance frequency vs sqrt(|J_afm|), confirming the exchange-enhanced AFMR scaling.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np

from cmtj import AxialDriver, CVector, Junction, Layer, ScalarDriver
from cmtj.utils.procedures import compute_spectrum_strip

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Shared AFM layer parameters
Ms = 0.5  # T
Ku = 5e4  # J/m^3, easy axis along x (the Neel axis)
Kdir = CVector(1, 0, 0)
thickness = 1e-9  # m
cell_surface = 1e-16  # m^2
damping = 0.005
demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]

dt = 1e-14  # s -- resolves oscillations up to a few THz (Nyquist ~5 THz)
sim_time = 5e-9  # s
pulse_amplitude = 5e6  # A/m -- large transient kick (see mechanism 2 above)
pulse_duration = 5 * dt


def make_afm_layer(j_afm: float) -> Layer:
    layer = Layer.createAFMLayer(
        "afm",
        CVector(1, 0, 0),
        CVector(-1, 0, 0),
        Kdir,
        Ms,
        thickness,
        cell_surface,
        demag,
        damping,
        ScalarDriver.getConstantDriver(j_afm),
    )
    layer.setAnisotropyDriver(ScalarDriver.getConstantDriver(Ku))
    layer.setExternalFieldDriver(
        AxialDriver(
            ScalarDriver.getConstantDriver(0),
            ScalarDriver.getStepDriver(0, pulse_amplitude, 0, pulse_duration),
            ScalarDriver.getConstantDriver(0),
        )
    )
    return layer


def resonance_frequency(j_afm: float) -> tuple[float, dict]:
    """Run the pulsed-excitation simulation and extract the dominant ring-down frequency."""
    junction = Junction([make_afm_layer(j_afm)])
    junction.runSimulation(sim_time, dt, dt)
    log = junction.getLog()
    mz = np.asarray(log["afm_mz"])
    spectrum, freqs = compute_spectrum_strip(mz - np.mean(mz), dt, max_frequency=5e12)
    peak_freq = float(freqs[np.argmax(spectrum)][0])
    return peak_freq, log


# --- Figure 2: resonance frequency vs sqrt(|J_afm|) exchange sweep ---
# Sign convention: NEGATIVE J_afm = antiferromagnetic (Neel) coupling, which is
# what this layer models. Positive J would be ferromagnetic intra-layer coupling,
# making the antiparallel initial state unstable and the resonance analysis
# meaningless -- keep the sweep strictly negative.
j_values = -np.linspace(0.02, 2.0, 100)

# --- Figure 1: antiparallel-locked precession trajectory at the two sweep edges ---
# Only the first slice is plotted: the ring-down decays (damping x resonance
# frequency) within a fraction of a ns, and the rest of the 3 ns window (needed
# for good FFT frequency resolution) is just flat equilibrium. Weak and strong
# exchange ring down on very different timescales (faster ring-down = higher
# |J_afm|), so each edge gets its own time window rather than sharing one axis.
j_edges = [j_values[0], j_values[-1]]
plot_windows_ns = [0.08, 0.008]  # weak exchange rings down slower than strong
edge_trajectories = []
for j_edge, window_ns in zip(j_edges, plot_windows_ns):
    _, edge_log = resonance_frequency(j_edge)
    time_ns = np.asarray(edge_log["time"]) * 1e9
    plot_mask = time_ns <= window_ns
    edge_trajectories.append(
        (
            time_ns[plot_mask],
            np.asarray(edge_log["afm_mz"])[plot_mask],
            np.asarray(edge_log["afm_m2z"])[plot_mask],
        )
    )
peak_freqs = np.array([resonance_frequency(j)[0] for j in j_values])
sqrt_j = np.sqrt(np.abs(j_values))

# linear fit through the origin: f = a * sqrt(|J_afm|)
slope = float(np.sum(sqrt_j * peak_freqs) / np.sum(sqrt_j**2))
correlation = float(np.corrcoef(sqrt_j, peak_freqs)[0, 1])

with plt.style.context(["science", "no-latex"]):
    fig, (ax1a, ax1b, ax2) = plt.subplots(1, 3, figsize=(14, 4), dpi=300)

    edge_axes = [ax1a, ax1b]
    edge_labels = ["Weak exchange", "Strong exchange"]
    for ax, (t_edge, m1z, m2z), j_edge, label in zip(edge_axes, edge_trajectories, j_edges, edge_labels):
        ax.plot(t_edge, m1z, color="crimson", linewidth=1.2, label="Sublattice 1 $m_z$")
        ax.plot(t_edge, m2z, color="navy", linewidth=1.2, label="Sublattice 2 $m_z$")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("$m_z$")
        ax.set_title(f"{label}\n$J_\\mathrm{{afm}}$ = {j_edge:.2f} J/m$^2$")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    ax2.plot(sqrt_j, peak_freqs / 1e9, "o", color="forestgreen", markersize=6, label="Simulation")
    fit_x = np.linspace(0, sqrt_j.max(), 50)
    ax2.plot(
        fit_x, slope * fit_x / 1e9, "--", color="gray", linewidth=1, label=f"$f = a\\sqrt{{|J|}}$, r={correlation:.3f}"
    )
    ax2.set_xlabel(r"$\sqrt{|J_\mathrm{afm}|}$ ($\sqrt{\mathrm{J/m^2}}$)")
    ax2.set_ylabel("Resonance frequency (GHz)")
    ax2.set_title("Exchange-enhanced AFMR scaling")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(
        "./curated-examples/figures/afm-resonance.png",
        dpi=300,
        bbox_inches="tight",
    )

print(f"J_afm sweep: {j_values} J/m^2")
print(f"Resonance frequencies: {np.round(peak_freqs / 1e9, 2)} GHz")
print(f"Linear fit slope (f = a*sqrt(|J|)): {slope / 1e9:.2f} GHz per sqrt(J/m^2)")
print(f"Correlation (f vs sqrt(|J_afm|)): {correlation:.4f}")
assert correlation > 0.95, "AFMR frequency should scale ~linearly with sqrt(|J_afm|)"
print("AFMR exchange-enhancement scaling confirmed.")
