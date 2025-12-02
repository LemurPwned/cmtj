"""
Spin-Transfer Oscillator (STO) Simulation

GOAL:
This simulation studies spin-transfer oscillators (STOs) - magnetic devices that generate
sustained microwave oscillations through spin-transfer torque effects. The simulation
demonstrates how current-driven spin-transfer torque can overcome magnetic damping to
create persistent magnetization oscillations, which are essential for microwave signal
generation and neuromorphic computing applications.

CRITICAL SIMULATION MECHANISMS:
1. Spin-Transfer Torque (STT): Models current-induced torque using Slonczewski formalism
   with both parallel and perpendicular torque components that can overcome damping
2. Single-Layer System: Implements a free layer with fixed reference magnetization
   direction to provide the necessary torque asymmetry for sustained oscillations
3. Nonlinear Dynamics: Captures the transition from damped precession to sustained
   oscillations as current exceeds the critical threshold
4. Magnetization Trajectories: Tracks 3D magnetization dynamics on the unit sphere
   to visualize precession patterns and limit cycle behavior

The simulation generates magnetization trajectories that demonstrate the STO operation
principle. This enables optimization of device parameters for applications in microwave
sources, frequency synthesizers, and neuromorphic computing where controlled oscillatory
behavior is essential.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np

from cmtj import CVector, Junction, Layer, constantDriver

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

r = 1
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0.0:pi:100j, 0.0 : 2.0 * pi : 100j]
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)


def plot_trajectories(log: dict[str, list[float]], title: str):
    with plt.style.context(["science", "no-latex"]):
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(1, 2, 1, projection="3d")
        m = np.asarray([log["free_mx"], log["free_my"], log["free_mz"]])
        ax.set_axis_off()
        ax.plot_surface(
            x,
            y,
            z,
            rstride=2,
            cstride=2,
            alpha=0.3,
            linewidth=0.1,
            color="lavender",
        )
        ax.plot3D(m[0], m[1], m[2], color="crimson", linewidth=3.0)
        ax.scatter([0], [0], [1], color="forestgreen", alpha=1.0, s=50, label="Start")
        ax.legend()
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(
            np.asarray(log["time"]) * 1e9,
            np.asarray(log["R"]),
            color="navy",
            linewidth=1.0,
        )
        ax2.set_xlabel("t (ns)")
        ax2.grid(False)

        for spine in ax2.spines.values():
            spine.set_visible(False)
        ax2.set_ylabel(r"R ($\Omega$)")
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig("./curated-examples/figures/sto.png", dpi=300, bbox_inches="tight")


demagTensor = [CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 0.0), CVector(0.0, 0.0, 1.0)]

damping = 0.3
currentDensity = -6e9
beta = 1
spinPolarisation = 1.0

l1 = Layer.createSTTLayer(
    id="free",
    mag=CVector(0.0, 0.0, 1.0),
    anis=CVector(0, 0.0, 1.0),
    Ms=1.0,
    thickness=1.4e-9,
    cellSurface=7e-10 * 7e-10,
    demagTensor=demagTensor,
    damping=damping,
    SlonczewskiSpacerLayerParameter=1.0,
    spinPolarisation=spinPolarisation,
    beta=beta,
)

l1.setReferenceLayer(CVector(0, 1.0, 1.0))
junction = Junction([l1], 100, 200)

junction.setLayerAnisotropyDriver("free", constantDriver(350e3))
# current driver
junction.setLayerCurrentDriver("free", constantDriver(currentDensity))
junction.runSimulation(50e-9, 1e-13, 1e-13)
log = junction.getLog()
plot_trajectories(log, title="STO with STT on")
