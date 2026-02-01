"""
Domain Wall Motion Simulation

GOAL:
This simulation demonstrates domain wall (DW) motion in magnetic nanowires driven by
spin-orbit torque (SOT) and external magnetic fields. It models the dynamics of a
magnetic domain wall using a collective coordinate approach, tracking the wall's
position, internal angle, and width as it moves through the nanowire.

CRITICAL SIMULATION MECHANISMS:
1. Collective Coordinate Model: Represents domain wall with three degrees of freedom -
   position X, internal angle φ (Néel vs Bloch character), and width δ
2. Spin-Orbit Torque: Models current-induced domain wall motion via spin Hall effect
   with both damping-like and field-like torque components
3. DMI (Dzyaloshinskii-Moriya Interaction): Stabilizes Néel-type domain walls and
   determines the wall chirality through interfacial symmetry breaking
4. Pinning Landscape: Includes artificial pinning potential to model defects and
   structural inhomogeneities that oppose wall motion
5. Dynamic Wall Width: Allows wall width to adapt during motion, capturing
   non-trivial dynamics under strong driving forces
6. Field Contributions: Includes perpendicular anisotropy, in-plane anisotropy,
   external fields, and demagnetization effects

The simulation generates time-series data showing wall position, velocity, and
internal angle evolution. This enables determination of wall mobility, threshold
currents for depinning, and wall structure transformations under various driving
conditions - essential for optimizing racetrack memory and logic devices.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np

from cmtj.models.domain_dynamics import (
    DW,
    DomainWallDynamics,
    MultilayerWallDynamics,
)
from cmtj.utils.general import VectorObj

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Material parameters for Co/Pt nanowire with strong DMI
# NOTE: domain_dynamics module uses Ms in A/m (not T as in core Layer objects)
Ms = 8e5  # Saturation magnetization [A/m]
alpha = 0.03  # Gilbert damping (typical range: 0.01-0.03)
thickness = 1e-9  # Layer thickness [m]
Lx = 2000e-9  # Nanowire length [m]
Ly = 100e-9  # Nanowire width [m]
Lz = thickness  # Nanowire height [m]

# Magnetic interactions
Ku = 8e5  # Perpendicular magnetic anisotropy [J/m^3]
Kp = 0  # In-plane anisotropy [J/m^3]
D = 2e-3  # DMI constant [J/m^2] - stabilizes Néel walls
A = 15e-12  # Exchange stiffness [J/m]

# Spin-orbit torque parameters
SHE_angle = 0.15  # Spin Hall angle
beta = 0.0  # Non-adiabatic STT parameter

# Pinning potential parameters
V0_pin = 1e-20  # Pinning potential depth [J]
pinning = 200e-9  # Pinning period [m]

# External field along z (perpendicular to film)
Hz = 0  # External field [A/m]
H_ext = VectorObj(0, 0, Hz)

# Create domain wall dynamics object
dw = DomainWallDynamics(
    H=H_ext,
    alpha=alpha,
    Ms=Ms,
    thickness=thickness,
    Lx=Lx,
    Ly=Ly,
    Lz=Lz,
    SHE_angle=SHE_angle,
    beta=beta,
    D=D,
    Ku=Ku,
    Kp=Kp,
    A=A,
    V0_pin=V0_pin,
    pinning=pinning,
)

# Define time-dependent current pulse
# Current density ramps up to drive the wall
def current_pulse(t):
    """Current density profile [A/m^2]"""
    t_ns = t * 1e9  # Convert to ns
    if t_ns < 2:
        return 0
    elif t_ns < 4:
        # Ramp up
        return 5e11 * (t_ns - 2) / 2
    elif t_ns < 10:
        # Hold constant
        return 5e11
    else:
        # Ramp down
        return max(0, 5e11 * (12 - t_ns) / 2)


dw.set_current_function(current_pulse)

# Create multilayer wrapper (single wall, no interlayer coupling)
mlw = MultilayerWallDynamics(layers=[dw], J=0)

# Initial conditions: wall at center, Néel-right configuration, natural width
X0 = Lx / 2  # Start at center of nanowire
phi0 = DW.NEEL_RIGHT  # Néel wall pointing right
delta0 = dw.dw0  # Natural DW width
starting_conditions = [X0, phi0, delta0]

# Run simulation
sim_time = 15e-9  # 15 ns simulation
results = mlw.run(sim_time=sim_time, starting_conditions=starting_conditions)

# Extract results
time = np.array(results["t"]) * 1e9  # Convert to ns
position = np.array(results["x_0"]) * 1e9  # Convert to nm
phi = np.array(results["phi_0"])
dw_width = np.array(results["dw_0"]) * 1e9  # Convert to nm
velocity = np.array(results["v_0"])  # m/s
current = np.array(results["je_0"]) * 1e-11  # Scale to 10^11 A/m^2

# Create simplified figure - focus on position and velocity
with plt.style.context(["science", "no-latex"]):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), dpi=300)

    # Plot 1: Wall position vs time
    ax1.plot(time, position, color="crimson", linewidth=2)
    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("Position (nm)")
    ax1.set_title(f"DW Position: $\\Delta x$ = {(position[-1] - position[0]):.1f} nm")
    ax1.grid(True, alpha=0.3)

    # Plot 2: Wall velocity vs time
    ax2.plot(time, velocity, color="navy", linewidth=2)
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel(r"Velocity (m/s)")
    ax2.set_title(f"$M_s$ = {Ms/1e6:.1f} MA/m, $D$ = {D*1e3:.1f} mJ/m², $\\theta_{{SH}}$ = {SHE_angle:.2f}")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(
        "./curated-examples/figures/dw-motion.png",
        dpi=300,
        bbox_inches="tight",
    )

print(f"Domain wall moved {(position[-1] - position[0]):.1f} nm in {sim_time*1e9:.1f} ns")
print(f"Average velocity: {np.mean(np.abs(velocity)):.2f} m/s")
print(f"Maximum velocity: {np.max(np.abs(velocity)):.2f} m/s")
