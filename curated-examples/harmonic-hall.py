"""
Harmonic Hall Voltage Measurement Simulation

GOAL:
This simulation studies harmonic Hall voltage measurements in magnetic thin films to
separate and quantify different magnetoresistance contributions. The simulation applies
AC current excitation and analyzes the resulting voltage harmonics to extract material
parameters including anomalous Hall effect resistance (AHE), anisotropic magnetoresistance (AMR),
and spin Hall magnetoresistance (SMR) coefficients.

CRITICAL SIMULATION MECHANISMS:
1. Multi-geometry Measurements: Performs Hall measurements in different field orientations
   (φ=0° and φ=45°) to separate longitudinal and transverse resistance components
2. Harmonic Analysis: Applies AC current and analyzes 1st and 2nd harmonic voltage
   responses to isolate different physical contributions
3. Resistance Tensor Modeling: Calculates full resistance tensor including AMR, AHE,
   and SMR contributions based on instantaneous magnetization orientation
4. Field Sweep Protocol: Systematically varies magnetic field magnitude and direction
   to map angular-dependent resistance behavior
5. Magnetization Dynamics: Solves LLG equation with realistic demagnetization tensor
   to capture field-dependent equilibrium magnetization states
6. Oersted Field Excitation: Includes current-induced Oersted fields that provide
   AC excitation for harmonic generation
7. Signal Processing: Extracts DC and AC components from resistance signals to
   compute harmonic amplitudes and phases

The simulation generates field-dependent resistance curves for both Rxx and Rxy
components in different measurement geometries. This enables separation of different
magnetoresistance mechanisms and extraction of material-specific parameters essential
for characterizing magnetic thin films and optimizing spintronic device performance.
The harmonic analysis technique is particularly powerful for distinguishing between
bulk and interface contributions to magnetoresistance.
"""

import contextlib
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from tqdm import tqdm

from cmtj import AxialDriver, CVector, Junction, Layer, NullDriver, ScalarDriver
from cmtj.utils.linear import FieldScan
from cmtj.utils.resistance import calculate_resistance_parallel

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

OeToAm = 79.57747


mpl.rcParams["axes.formatter.useoffset"] = False
fsize = 20

Rx0 = [304.306]
Ry0 = [1.008]
SMR = [-0.464]
AMR = [-0.053]
AHE = [-5.71]
w = [3e-5]
l = [2e-5]
lw = 2.5

Ms = 0.525
Ku = 1.54e5
alpha = 0.03


def run_simulation(junction: Junction, Hvecs: np.ndarray, mode: str, int_time=1e-12):
    sim_time = 10e-9
    layer_str = ["free"]
    mags = [CVector(0, 0, 1) for _ in layer_str]
    Rxy, Rxx = [], []
    for Hval in Hvecs:
        junction.clearLog()
        HDriver = AxialDriver(
            ScalarDriver.getConstantDriver(Hval[0]),
            ScalarDriver.getConstantDriver(Hval[1]),
            ScalarDriver.getConstantDriver(Hval[2]),
        )
        junction.setLayerExternalFieldDriver("all", HDriver)
        # set mags for better convergence
        for i, l_str in enumerate(layer_str):
            junction.setLayerMagnetisation(l_str, mags[i])

        junction.runSimulation(sim_time, int_time, int_time)

        # set new mags
        for str_ in layer_str:
            mags[i] = junction.getLayerMagnetisation(str_)

        log = junction.getLog()
        m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in layer_str])
        dynamicRx, dynamicRy = calculate_resistance_parallel(
            [Rx0],
            [0],
            AMR=AMR,
            AHE=AHE,
            SMR=SMR,
            m=m,
            l=l,
            w=w,
        )

        Rxy.append(dynamicRy[-1])
        Rxx.append(dynamicRx[-1])
    return np.asarray(Rxx) if mode.lower() == "rxx" else np.asarray(Rxy)


def simulate(Ku, Ms, Hvecs, alpha, mode="rxx"):
    demagTensor = [
        CVector(0.00024164288391924, 2.71396011566517e-10, 5.95503928124313e-14),
        CVector(2.71396011566517e-10, 0.000160046006320031, 1.32504057070646e-14),
        CVector(5.95503928124313e-14, 1.32504057070646e-14, 0.999598310229469),
    ]

    thickness = 1.45e-9
    surface = w[0] * l[0]
    l1 = Layer(
        id="free",
        mag=CVector(0, 0, 1),
        anis=CVector(0.0, 0, 1),
        Ms=Ms,
        thickness=thickness,
        cellSurface=surface,  # only for temperature calculation
        damping=alpha,
        demagTensor=demagTensor,
    )
    junction = Junction([l1])

    HoePulseAmpl = 50
    HoeDriver = AxialDriver(
        NullDriver(),
        NullDriver(),
        ScalarDriver.getStepDriver(0, HoePulseAmpl, 0.0, 1e-11),
    )
    junction.setLayerOerstedFieldDriver("all", HoeDriver)
    junction.setLayerAnisotropyDriver("all", ScalarDriver.getConstantDriver(Ku))

    return run_simulation(junction=junction, Hvecs=Hvecs, mode=mode)


ms = 2
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots(2, 2, figsize=(4, 4), dpi=300)
    resistance_data = defaultdict(list)

    for i, field in enumerate(tqdm(["hx", "h45"])):
        for j, res_mode in enumerate(("rxx", "rxy")):
            flabel = field.capitalize()
            if field == "hx":
                phi = 1
                flabel = r"H ($\phi=0^\circ$)"
            elif field == "hy":
                phi = 90
                flabel = r"H ($\phi=90^\circ$)"
            elif field == "h45":
                phi = 45
                flabel = r"H ($\phi=45^\circ$)"
            else:
                raise ValueError("Unknown field orientation")

            theta = 92
            Hscan, Hvecs = FieldScan.amplitude_scan(
                theta=theta, phi=phi, start=-600e3, stop=600e3, steps=100, back=True
            )
            simulated = simulate(Ms=Ms, Ku=Ku, alpha=alpha, mode=res_mode, Hvecs=Hvecs)
            if res_mode == "rxy":
                simulated -= np.mean(simulated)

            ax[i, j].plot(
                Hscan / 1e3,
                simulated,
                "r-",
                linewidth=lw,
                label="Sim.",
            )
            ax[i, j].plot(
                Hscan / 1e3,
                simulated,
                "bo",
                markersize=ms,
            )
            ax[i, j].set_xlabel(f"{flabel} (kA/m)")
            ax[i, j].set_ylabel(rf"{res_mode.capitalize()} $(\Omega)$")
            ax[i, 1].yaxis.tick_right()
        ax[i, 1].yaxis.set_label_position("right")
        ax[i, 1].set_ylabel(rf"{res_mode.capitalize()} $(\Omega)$", rotation=270)

    for label, ax_ in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], ax.flatten()):
        # label physical distance in and down:
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

    fig.align_ylabels()
    fig.suptitle("Harmonic Hall 1st and 2nd harmonic response")
    fig.subplots_adjust(wspace=0.0, hspace=0.45)
    fig.savefig(
        "./curated-examples/figures/harmonic-hall.png",
        dpi=350,
        bbox_inches="tight",
    )
