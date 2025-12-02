"""
Ferromagnetic Resonance (FMR) Spectroscopy Simulation

GOAL:
This simulation performs ferromagnetic resonance (FMR) simulation on a multi-layer
magnetic system to characterize resonance modes and extract magnetic parameters.
The simulation maps the frequency-field response by applying pulsed excitation and
analyzing the resulting magnetization dynamics through Fourier analysis.

CRITICAL SIMULATION MECHANISMS:
1. Multi-layer System: Models three magnetic layers (free, perpendicular, buffer) with
   different magnetic properties (Ms, Ku) and coupling interactions
2. Pulsed Excitation: Applies short-duration step pulses in Oersted field to excite
   magnetization oscillations across a broad frequency spectrum
3. Field Sweep Protocol: Systematically varies external magnetic field to map
   field-dependent resonance conditions for each layer
4. Fourier Analysis: Performs FFT on time-domain magnetization signals to extract
   frequency-domain response and identify resonance peaks
5. Magnetization Dynamics: Solves coupled LLG equations for all layers including
   interlayer exchange coupling and demagnetization effects
6. Spectral Processing: Applies frequency filtering and signal processing to enhance
   resonance peak visibility and reduce noise

The simulation generates frequency-field spectrograms that reveal FMR modes,
enabling determination of layer-specific magnetic parameters (Ms, Ku, damping) and
characterization of interlayer coupling strength.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.ndimage import uniform_filter
from tqdm import tqdm

from cmtj import (
    AxialDriver,
    CVector,
    Junction,
    Layer,
    NullDriver,
    constantDriver,
    stepDriver,
)
from cmtj.utils import FieldScan

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

Ms1 = 1.03
Ms2 = 1.03
Ms3 = 0.12

Ku1 = 489e3
Ku2 = 514e3
Ku3 = 17e3
Kdir = CVector(0, 0, 1)
damping = 0.01
demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]
demag_buffer = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]

l1 = Layer(
    "free",
    mag=CVector(0, 0, 0.99),
    anis=Kdir,
    Ms=Ms1,
    thickness=1.0e-9,
    damping=damping,
    demagTensor=demag,
    cellSurface=1e-9,
)
l2 = Layer(
    "perpendicular",
    mag=CVector(0, 0, 0.99),
    anis=Kdir,
    Ms=Ms2,
    thickness=1.0e-9,
    damping=damping,
    cellSurface=1e-9,
    demagTensor=demag,
)
l3 = Layer(
    "buffer",
    mag=CVector(1, 0, 0),
    anis=CVector(1, 0, 0),
    Ms=Ms3,
    thickness=1.0e-9,
    damping=0.015,
    cellSurface=1e-9,
    demagTensor=demag_buffer,
)
l1.setAnisotropyDriver(constantDriver(Ku1))
l2.setAnisotropyDriver(constantDriver(Ku2))
l3.setAnisotropyDriver(constantDriver(Ku3))
j = Junction([l1, l2, l3])

Hmax = -700e3
Hscan, Hvecs = FieldScan.amplitude_scan(
    start=-Hmax,
    stop=Hmax,
    steps=100,
    theta=0,
    phi=0,
)
dt = 5e-13
sim_time = 60e-9
ampl = 500
dur = 2e-9
spectrum = []
wait_time = 4e-9
j.setLayerExternalFieldDriver("all", AxialDriver(*Hvecs[0]))
j.runSimulation(sim_time, dt, dt)
back = False


def compute_fft(mixed_signal: np.ndarray, dt: float) -> np.ndarray:
    fft_mixed = np.abs(fft(mixed_signal))
    ffreqs = fftfreq(len(mixed_signal), dt)
    max_freq = 30e9
    fft_mixed = fft_mixed[1 : len(fft_mixed) // 2]
    ffreqs = ffreqs[1 : len(ffreqs) // 2]
    fft_mixed = fft_mixed[np.where(ffreqs < max_freq)]
    ffreqs = ffreqs[np.where(ffreqs < max_freq)]
    return fft_mixed, ffreqs


def extract_fft(log: dict[str, list[float]], dt: float, wait_time: float) -> np.ndarray:
    tm = np.asarray(log["time"])
    tmindx = np.argwhere(tm > wait_time).ravel()
    my_free = np.asarray(log["free_mx"])[tmindx]
    my_perpendicular = np.asarray(log["perpendicular_mx"])[tmindx]
    buffer = np.asarray(log["buffer_mx"])[tmindx]
    mmixed = (1 / 3) * (my_free + my_perpendicular + buffer)
    fft_mixed, ffreqs = compute_fft(mmixed, dt)
    return fft_mixed, ffreqs


for Hv in tqdm(Hvecs):
    j.clearLog()
    j.setLayerExternalFieldDriver("all", AxialDriver(*Hv))
    j.setLayerOerstedFieldDriver("all", AxialDriver(NullDriver(), stepDriver(0, ampl, 0, dur), NullDriver()))
    j.runSimulation(sim_time, dt, dt)
    log = j.getLog()
    fft_mixed, ffreqs = extract_fft(log, dt, wait_time)

    spectrum.append(fft_mixed)

if back:
    for i, Hv in tqdm(enumerate(Hvecs[::-1])):
        j.clearLog()
        j.setLayerExternalFieldDriver("all", AxialDriver(*Hv))
        j.setLayerOerstedFieldDriver("all", AxialDriver(NullDriver(), stepDriver(0, ampl, 0, dur), NullDriver()))
        j.runSimulation(sim_time, dt, dt)
        log = j.getLog()
        fft_mixed, ffreqs = extract_fft(log, dt, wait_time)

        spectrum[len(spectrum) - i - 1] += fft_mixed

spectrum = np.array(spectrum) / 2.0


with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots(dpi=300)
    ax.pcolormesh(Hscan / 1e3, ffreqs / 1e9, 10 * np.log10(uniform_filter(spectrum.T, size=1)))
    ax.set_xlabel(r"$\mathrm{H}$ ($\mathrm{kA/m}$)")
    ax.set_ylabel(r"$f$ ($\mathrm{GHz}$)")
    ax.legend(loc="upper right", fontsize=6)
    fig.suptitle(
        rf"$\mu_0 M_\mathrm{{s,1}}$ = {Ms1:.2f} T, $\mu_0 M_\mathrm{{s,2}}$ = {Ms2:.2f} T, "
        rf"$\mu_0 M_\mathrm{{s,3}}$ = {Ms3:.2f} T"
        "\n"
        rf"$K_\mathrm{{u,1}}$ = {Ku1 / 1e3:.1f} kA/m, $K_\mathrm{{u,2}}$ = {Ku2 / 1e3:.1f} kA/m, "
        rf"$K_\mathrm{{u,3}}$ = {Ku3 / 1e3:.1f} kA/m",
        fontsize=6,
    )
    fig.savefig(
        "./curated-examples/figures/fmr-basic.png",
        dpi=350,
        bbox_inches="tight",
    )
