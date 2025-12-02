"""
Voltage Spin Diode (VSD) Frequency Scan Simulation

GOAL:
This simulation demonstrates Voltage Spin Diode (VSD) by performing frequency sweeps to
characterize the dynamic response of a free magnetic layer. The simulation computes the
DC voltage generated through the mixing of AC current and dynamic resistance changes due
to magnetization precession.

CRITICAL SIMULATION MECHANISMS:
1. Magnetization Dynamics: Uses Landau-Lifshitz-Gilbert (LLG) equation with RK4 solver
   to simulate magnetization precession in a perpendicular magnetic anisotropy (PMA) layer
2. AC Excitation: Applies sinusoidal Oersted field to drive magnetization oscillations
3. Resistance Calculation: Computes dynamic resistance using parallel resistance model
   incorporating AMR (Anisotropic Magnetoresistance), SMR (Spin Hall Resistance),
   and AHE (Anomalous Hall Effect) contributions
4. Voltage Mixing: Calculates DC voltage through multiplication of dynamic current and
   dynamic resistance, followed by low-pass filtering to extract DC component
5. Frequency Response: Sweeps excitation frequency to map resonance characteristics
6. Field Orientation: Compares 4-point (φ=0°) and 2-point (φ=45°) measurement geometries

The simulation outputs frequency-dependent voltage spectra that reveal ferromagnetic
resonance (FMR) characteristics and can be used to extract material parameters like
damping, anisotropy, and saturation magnetization.
"""

import contextlib
import math
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from tqdm import tqdm

import cmtj
from cmtj import AxialDriver, CVector, Junction, Layer, NullDriver, constantDriver, sineDriver
from cmtj.utils import Filters
from cmtj.utils.resistance import calculate_resistance_parallel

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

OeToAm = 79.57747


Rx0 = 304.306
Ry0 = 1.008
SMR = -0.464
AMR = -0.053
AHE = -5.71
w = 3e-5
l = 2e-5
INT_STEP = 1e-13
HMIN = 350e3
HMAX = 630e3
HSTEPS = 100


def compute_vsd2(dynamicR, integration_step, dynamicI):
    """Pymag"""
    SD = -dynamicI * dynamicR
    fs = 1.0 / integration_step
    SD_dc = Filters.butter_lowpass_filter(SD, cutoff=10e6, fs=fs, order=3)
    return np.mean(SD_dc)


def simulate_lorentz(Ms, Ku, frequency, orient, alpha=1e-4, Irf=0.5e-3):
    demagTensor = [
        CVector(0.0, 0.0, 0.0),
        CVector(0.0, 0.0, 0.0),
        CVector(0.0, 0.0, 1.0),
    ]
    thickness = 1.45e-9
    l1 = Layer(
        id="free",
        mag=CVector(0.0, 0.0, 1.0),
        anis=CVector(0, 0.0, 1.0),
        Ms=Ms,
        thickness=thickness,
        cellSurface=0,
        demagTensor=demagTensor,
        damping=alpha,
    )
    junction = Junction([l1])

    junction.setLayerAnisotropyDriver("free", constantDriver(Ku))

    HoeAmpl = 5000  # A/m
    Hspace = np.linspace(HMIN, HMAX, num=HSTEPS)
    theta = np.deg2rad(91)
    if orient == "4p":
        phideg = 0
    elif orient == "2p":
        phideg = 45
    else:
        raise ValueError("Unknown orient")
    phi = np.deg2rad(phideg)
    Hsweep = np.zeros(Hspace.shape[0])
    for i, H in enumerate(Hspace):
        junction.clearLog()
        HDriver = AxialDriver(
            constantDriver(H * math.sin(theta) * math.cos(phi)),
            constantDriver(H * math.sin(theta) * math.sin(phi)),
            constantDriver(H * math.cos(theta)),
        )

        HoeDriver = AxialDriver(
            NullDriver(),
            NullDriver(),
            sineDriver(0, -HoeAmpl, frequency, 0),
        )
        junction.setLayerExternalFieldDriver("all", HDriver)
        junction.setLayerOerstedFieldDriver("all", HoeDriver)
        junction.runSimulation(40e-9, INT_STEP, INT_STEP, solverMode=cmtj.RK4)

        log = junction.getLog()
        m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ["free"]])
        dynamicRx, dynamicRy = calculate_resistance_parallel(
            Rx0=[Rx0],
            Ry0=[Ry0],
            AMR=[AMR],
            AHE=[AHE],
            SMR=[SMR],
            m=m,
            l=[l],
            w=[w],
        )
        dynamicR = dynamicRx if orient == "2p" else dynamicRy
        dynamicI = Irf * np.sin(2 * math.pi * frequency * np.asarray(log["time"]))
        vmix = compute_vsd2(dynamicR, INT_STEP, dynamicI)
        Hsweep[i] = vmix
    return Hspace, Hsweep


def simulate_lorentz_freq_scan(Ms, Ku, Hvalue, orient, alpha=1e-4, Irf=0.5e-3):
    demagTensor = [
        CVector(0.0, 0.0, 0.0),
        CVector(0.0, 0.0, 0.0),
        CVector(0.0, 0.0, 1.0),
    ]
    thickness = 1.45e-9
    l1 = Layer(
        id="free",
        mag=CVector(0.0, 0.0, 1.0),
        anis=CVector(0, 0.0, 1.0),
        Ms=Ms,
        thickness=thickness,
        cellSurface=0,
        demagTensor=demagTensor,
        damping=alpha,
    )
    junction = Junction([l1])

    junction.setLayerAnisotropyDriver("free", constantDriver(Ku))

    HoeAmpl = 5000  # A/m
    theta = np.deg2rad(91)
    if orient == "4p":
        phideg = 0
    elif orient == "2p":
        phideg = 45
    else:
        raise ValueError("Unknown orient")
    phi = np.deg2rad(phideg)
    fspace = np.arange(2, 18, 0.1) * 1e9
    fsweep = np.zeros(fspace.shape[0])
    for i, f in enumerate(fspace):
        junction.clearLog()
        HDriver = AxialDriver(
            constantDriver(Hvalue * math.sin(theta) * math.cos(phi)),
            constantDriver(Hvalue * math.sin(theta) * math.sin(phi)),
            constantDriver(Hvalue * math.cos(theta)),
        )

        HoeDriver = AxialDriver(
            NullDriver(),
            NullDriver(),
            sineDriver(0, -HoeAmpl, f, 0),
        )
        junction.setLayerExternalFieldDriver("all", HDriver)
        junction.setLayerOerstedFieldDriver("all", HoeDriver)
        junction.runSimulation(40e-9, INT_STEP, INT_STEP, solverMode=cmtj.RK4)

        log = junction.getLog()
        m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ["free"]])
        dynamicRx, dynamicRy = calculate_resistance_parallel(
            Rx0=[Rx0],
            Ry0=[Ry0],
            AMR=[AMR],
            AHE=[AHE],
            SMR=[SMR],
            m=m,
            l=[l],
            w=[w],
        )
        dynamicR = dynamicRx if orient == "4p" else dynamicRy
        dynamicI = Irf * np.sin(2 * math.pi * f * np.asarray(log["time"]))
        vmix = compute_vsd2(dynamicR, INT_STEP, dynamicI)
        fsweep[i] = vmix
    return fspace, fsweep


data = defaultdict(list)
hscans = []
vscans = []
fmin = 12
fmax = 17
fscan = np.arange(12, 17)
alpha = 30e-3
Ms = 0.525
Ku = 1.54e5
for orient, irf in zip(("4p", "2p"), (0.75e-3, 0.4e-3)):
    for f in tqdm(fscan):
        hscan, vscan = simulate_lorentz(Ms, Ku, f * 1e9, orient=orient, alpha=alpha, Irf=irf)
        if orient == "2p":
            vscan -= vscan.max()
        data[f"{orient}"].append(vscan)
        data[f"{orient}-field"].append(hscan)


with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots(2, 1, dpi=300, sharex="col")
    colors = plt.get_cmap("rainbow")(np.linspace(0, 1, len(fscan)))
    for j, orient in enumerate(("4p", "2p")):
        for i, f in enumerate(fscan):
            ax[j].plot(
                np.asarray(data[f"{orient}-field"][i]) / 1e3,
                np.asarray(data[orient][i]) * 1e6,
                alpha=1,
                linestyle="--",
                label=f"{f} GHz",
                color=colors[i],
                linewidth=1,
            )
        ax[j].set_xlim([300, 630])
        ax[j].set_ylabel(r"$V_\mathrm{DC} (\mathrm{\mu V})$")

    ax[1].set_xlabel(r"$\mathrm{H}$ ($\mathrm{kA/m}$)")
    ax[1].legend()

    ax[0].legend()
    for label, ax_ in zip(["(a)", "(b)"], ax.flatten()):
        trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
        ax_.text(
            0.0,
            1.0,
            label,
            transform=ax_.transAxes + trans,
            fontsize="medium",
            verticalalignment="top",
            bbox=dict(facecolor="none", edgecolor="none", pad=3.0),
        )
    fig.subplots_adjust(hspace=0)
    fig.align_ylabels()
    fig.tight_layout()

    fig.savefig(
        "./curated-examples/figures/vsd-basic-freq-scan.png",
        dpi=350,
        bbox_inches="tight",
    )
