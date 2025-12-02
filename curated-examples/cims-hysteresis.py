"""
Current-Induced Magnetization Switching (CIMS) Hysteresis Simulation

GOAL:
This simulation studies current-induced magnetization switching in spin-orbit torque (SOT)
devices by analyzing hysteresis loops under varying current densities and Oersted field
contributions. The simulation maps the stable magnetization states as a function of
applied current, revealing switching thresholds and stability regions.

CRITICAL SIMULATION MECHANISMS:
1. Spin-Orbit Torque (SOT): Models both damping-like and field-like torque components
   arising from spin-orbit coupling in heavy metal/ferromagnet bilayers
2. Current Sweep Protocol: Applies bidirectional current sweeps with step drivers to
   map complete hysteresis loops and identify switching thresholds
3. Oersted Field Effects: Includes current-induced Oersted fields that provide additional
   torque contributions and can assist or oppose SOT switching
4. Magnetization Dynamics: Solves extended LLG equation including SOT terms using
   adaptive time stepping for accurate switching dynamics
5. Stability Analysis: Monitors final magnetization states after current pulses to
   determine stable configurations and switching probability
6. Resistance Readout: Uses GMR (Giant Magnetoresistance) model to convert magnetization
   states to resistance values for experimental comparison
7. Multi-parameter Sweep: Varies Oersted field scaling to study the interplay between
   SOT and Oersted field contributions

The simulation generates current-resistance hysteresis loops that reveal switching
thresholds, coercive currents, and the influence of Oersted fields on switching
efficiency, providing insights for optimizing SOT-based memory devices.
"""

import contextlib
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from cmtj import AxialDriver, CVector, Junction, Layer, ScalarDriver, constantDriver, stepDriver
from cmtj.utils import TtoAm, compute_gmr

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401


demagTensor = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1)]
"""
Compute the current scan.
The facing area is thickness x width

Below, you can change the torque signs!
Depending on the Hfl sign the Oersted field will assist or oppose the SOT switching.
If sign is changed, or one of the torques is removed, adjust the current density accordingly.
"""
sgnHdl = 1
sgnHfl = -1

Hdl = sgnHdl * 0.000371 * TtoAm
Hfl = sgnHfl * 0.000433 * TtoAm
jden1 = 6e10
HDL = Hdl / jden1
HFL = Hfl / jden1

last_N = 100
dt = 6e-9
tstart = 1e-9

tstep = 5e-13
sim_time = 15e-9


t_hm = 6e-9
t_fm = 1.7e-9
Kdir1 = CVector(0, 0, 1.0)
pdir1 = CVector(0.0, 1.0, 0.0)
Ku1 = 0.3e6
Ms = 1.45
alpha = 0.0064
layer_free = Layer.createSOTLayer(
    id="free",
    mag=CVector(1, 0, 0),
    anis=Kdir1,
    Ms=Ms,
    thickness=t_fm,
    cellSurface=0.0,
    demagTensor=demagTensor,
    damping=alpha,
    dampingLikeTorque=HDL,
    fieldLikeTorque=HFL,
)

j = Junction([layer_free])
j.setLayerAnisotropyDriver("free", ScalarDriver.getConstantDriver(Ku1))
j.setLayerReferenceLayer("free", pdir1)


j.clearLog()
j.runSimulation(14e-9, tstep, tstep)
m_components = defaultdict(list)
hysteresis_scales = []
Hoescales = np.linspace(0, 2500, 3)
Hoescales = np.concatenate((-Hoescales, Hoescales[1:]))

Hoescales = np.asarray([-1000, -500, -100, 0, 100, 500, 1000])
Hoescales = [1, 0.5, 1e-1, 1e-2, 0, -1e-1, -1]
Hoescales = [1, 0, -1]

"""
Below uncomment for specific cases depending on which, damping-like or field-like, torque is used.
"""
if abs(HDL) and abs(HFL):
    i_max = [0.1e12 for _ in range(len(Hoescales))]  # both Hdl and Hfl only
elif not abs(HFL) and not abs(HDL):
    i_max = [1.0e12 for _ in range(len(Hoescales))]  # Hoe only
elif abs(HDL) and not abs(HFL):
    i_max = [3e11, 3e11, 1e12, 3e13, 3e13, 3e13, 3e11]  # Hdl only
elif abs(HFL) and not abs(HDL):
    i_max = [2.0e11 for _ in range(len(Hoescales))]  # Hfl only
else:
    raise ValueError("Invalid torque configuration")

jscans = []
for i, Hoescale in tqdm(enumerate(Hoescales), total=len(Hoescales)):
    v = CVector(0.001, 0.99, 0)
    v.normalize()
    j.setLayerMagnetisation("free", v)
    jscan = np.linspace(-i_max[i], i_max[i], 130)
    jscan = np.concatenate([jscan, jscan[::-1][1:]])

    jscans.append(jscan)
    hysteresis = []
    for current in jscan:
        j.clearLog()
        j.setLayerCurrentDriver("all", stepDriver(0, current, tstart, tstart + dt))
        j.setLayerOerstedFieldDriver(
            "all",
            AxialDriver(
                constantDriver(0),
                stepDriver(0, Hoescale * current * t_hm / 2, tstart, tstart + dt),
                constantDriver(0),
            ),
        )
        j.runSimulation(sim_time, tstep, tstep)
        log = j.getLog()
        m = np.asarray([[log[f"{str_}_mx"], log[f"{str_}_my"], log[f"{str_}_mz"]] for str_ in ("free",)])

        Rstable = compute_gmr(
            Rp=100,
            Rap=200,
            m1=m.squeeze()[:, -last_N:],
            m2=np.asarray([[0, 1, 0] for _ in range(last_N)]).transpose(),
        ).mean()
        str_ = "free"
        m_x_stable = np.mean(log[f"{str_}_mx"][-last_N:])
        m_y_stable = np.mean(log[f"{str_}_my"][-last_N:])
        m_z_stable = np.mean(log[f"{str_}_mz"][-last_N:])
        m_components["x"].append(m_x_stable)
        m_components["y"].append(m_y_stable)
        m_components["z"].append(m_z_stable)
        hysteresis.append(Rstable)

    hysteresis = np.asarray(hysteresis)

    hysteresis_scales.append(hysteresis)


with plt.style.context(["nature"]):
    fig, ax = plt.subplots(1, 1, dpi=250)
    # Set white background
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    # get a color scale
    colors = plt.cm.viridis(np.linspace(0, 1, len(hysteresis_scales)))
    for i, hysteresis in enumerate(hysteresis_scales):
        ax.plot(
            jscans[i] * 1e-12,
            hysteresis,
            label=f"Hoe = {Hoescales[i]}",
            color=colors[i],
        )
    ax.set_xlabel(r"$\mathrm{j}$ ($\mathrm{TA/m^2}$)")
    ax.set_ylabel(r"$\mathrm{R}_\mathrm{stable}$ ($\Omega$)")
    cbar = plt.colorbar(
        plt.cm.ScalarMappable(norm=plt.Normalize(min(Hoescales), max(Hoescales)), cmap=plt.cm.viridis),
        ax=ax,
        label="Hoe",
    )
    cbar.set_ticks(Hoescales)
    ax.set_title(
        rf"$H_\mathrm{{ext}} = 0, \alpha = {alpha}, \mu_0 "
        rf"M_\mathrm{{s}} = {Ms} \, \mathrm{{T}}$"
        "\n"
        rf"$H_\mathrm{{dl}} = {Hdl:.2f} \, \mathrm{{A/m}}, "
        rf"H_\mathrm{{fl}} = {Hfl:.2f} \, \mathrm{{A/m}}$"
        "\n"
        rf"$t_\mathrm{{fm}} = {t_fm * 1e9:.0f}\,  \mathrm{{nm}}, j_\mathrm{{den}} = "
        rf"{jden1 / 1e12:.2f} \, \mathrm{{TA/m^2}}$"
    )
    fig.savefig(
        "./curated-examples/figures/cims-hysteresis.png",
        dpi=350,
        bbox_inches="tight",
    )
