"""
Interlayer Exchange Coupling (IEC) Dynamics Simulation

GOAL:
This simulation demonstrates interlayer exchange coupling (IEC) between two magnetic
layers in a synthetic antiferromagnet (SAF) structure. It shows how RKKY-type coupling
mediated through a non-magnetic spacer layer creates coupled magnetization dynamics,
enabling synchronized switching and collective behavior essential for spintronic
memory and logic applications.

CRITICAL SIMULATION MECHANISMS:
1. IEC Interaction: Models both linear (J1) and quadratic (J2) coupling terms that
   create effective fields between layers proportional to their relative orientations
2. Dual-Layer System: Simulates two ferromagnetic layers with different magnetic
   properties coupled through IEC, creating coupled precession modes
3. Field-Driven Dynamics: Applies pulsed external field to excite both layers and
   observe their coupled response through magnetization trajectories
4. Antiferromagnetic Coupling: Uses negative J1 to favor antiparallel alignment,
   characteristic of SAF structures for improved thermal stability
5. Coupled Precession: Captures acoustic and optical modes of magnetization dynamics
   arising from symmetric and antisymmetric combinations of layer motions
6. Time-Domain Analysis: Tracks magnetization components to show phase relationships
   and coupling strength effects on dynamic response

The simulation generates magnetization trajectories and phase space plots showing
the coupled dynamics. This enables optimization of IEC strength for applications
requiring synchronized switching (MRAM write), collective oscillations (coupled STOs),
or enhanced stability (SAF-based sensors and memory).
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np

from cmtj import AxialDriver, CVector, Junction, Layer, constantDriver

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Layer parameters following AGENTS.md guidelines
# Ms in Tesla for core Layer objects (not A/m)
Ms1 = 1.0  # Free layer saturation magnetization [T]
Ms2 = 0.8  # Fixed layer saturation magnetization [T]

# Anisotropy values in typical range for PMA
Ku1 = 400e3  # Free layer anisotropy [J/m^3]
Ku2 = 300e3  # Fixed layer anisotropy [J/m^3]

# Anisotropy direction (perpendicular to film plane)
Kdir = CVector(0, 0, 1)

# Damping in recommended range
damping1 = 0.02
damping2 = 0.015

# Demagnetization tensors (perpendicular films)
demag1 = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]
demag2 = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]

# Create two coupled magnetic layers
l1 = Layer(
    "layer1",
    mag=CVector(0, 0, 1.0),  # Initially up
    anis=Kdir,
    Ms=Ms1,
    thickness=1.5e-9,
    damping=damping1,
    demagTensor=demag1,
    cellSurface=np.pi * (50e-9) ** 2,  # Circular cross-section
)

l2 = Layer(
    "layer2",
    mag=CVector(0, 0, -1.0),  # Initially down (antiparallel due to IEC)
    anis=Kdir,
    Ms=Ms2,
    thickness=1.2e-9,
    damping=damping2,
    demagTensor=demag2,
    cellSurface=np.pi * (50e-9) ** 2,
)

# Set anisotropy drivers
l1.setAnisotropyDriver(constantDriver(Ku1))
l2.setAnisotropyDriver(constantDriver(Ku2))

# Create junction with both layers
junction = Junction([l1, l2])

# Set IEC coupling (negative J1 for antiferromagnetic coupling)
# J in mJ/m^2, typical range ±0.001 to ±3.0 according to AGENTS.md
J_linear = -0.5e-3  # -0.5 mJ/m^2, antiferromagnetic coupling
J_quad = 0.1e-3  # 0.1 mJ/m^2, small quadratic term

junction.setIECDriver("layer1", "layer2", constantDriver(J_linear))
junction.setQuadIECDriver("layer1", "layer2", constantDriver(J_quad))

# Apply external field pulse to perturb the system
# Field strength in typical range (±0 – ±500e3 A/m)
field_amplitude = 100e3  # A/m
field_duration = 2e-9  # 2 ns pulse

# Create field pulse function
def field_pulse(t):
    """Pulsed field in z direction"""
    if t < field_duration:
        return field_amplitude
    else:
        return 0


from cmtj import ScalarDriver

# Set external field on both layers
junction.setLayerExternalFieldDriver(
    "layer1",
    AxialDriver(
        ScalarDriver.getConstantDriver(0),
        ScalarDriver.getConstantDriver(0),
        ScalarDriver.getPulseDriver(0, field_amplitude, 0, field_duration),
    ),
)

junction.setLayerExternalFieldDriver(
    "layer2",
    AxialDriver(
        ScalarDriver.getConstantDriver(0),
        ScalarDriver.getConstantDriver(0),
        ScalarDriver.getPulseDriver(0, field_amplitude, 0, field_duration),
    ),
)

# Run simulation with appropriate time step
# Following AGENTS.md: use dt=1e-12 for standard simulations
dt = 1e-12
sim_time = 20e-9  # 20 ns total simulation (within 1-500 ns range)
junction.runSimulation(sim_time, dt, dt)

# Get simulation log
log = junction.getLog()

# Extract data
time = np.array(log["time"]) * 1e9  # Convert to ns
m1x = np.array(log["layer1_mx"])
m1y = np.array(log["layer1_my"])
m1z = np.array(log["layer1_mz"])
m2x = np.array(log["layer2_mx"])
m2y = np.array(log["layer2_my"])
m2z = np.array(log["layer2_mz"])

# Plot results
with plt.style.context(["science", "no-latex"]):
    fig = plt.figure(figsize=(14, 10), dpi=300)

    # Create 3D trajectory plots
    ax1 = fig.add_subplot(2, 3, 1, projection="3d")
    ax1.plot(m1x, m1y, m1z, color="crimson", linewidth=1.5, alpha=0.8)
    ax1.scatter([m1x[0]], [m1y[0]], [m1z[0]], color="green", s=50, label="Start")
    ax1.scatter([m1x[-1]], [m1y[-1]], [m1z[-1]], color="red", s=50, label="End")
    ax1.set_xlabel("$m_x$")
    ax1.set_ylabel("$m_y$")
    ax1.set_zlabel("$m_z$")
    ax1.set_title("Layer 1 Trajectory")
    ax1.legend(fontsize=8)

    ax2 = fig.add_subplot(2, 3, 2, projection="3d")
    ax2.plot(m2x, m2y, m2z, color="navy", linewidth=1.5, alpha=0.8)
    ax2.scatter([m2x[0]], [m2y[0]], [m2z[0]], color="green", s=50, label="Start")
    ax2.scatter([m2x[-1]], [m2y[-1]], [m2z[-1]], color="red", s=50, label="End")
    ax2.set_xlabel("$m_x$")
    ax2.set_ylabel("$m_y$")
    ax2.set_zlabel("$m_z$")
    ax2.set_title("Layer 2 Trajectory")
    ax2.legend(fontsize=8)

    # Phase space plot showing coupling
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.plot(m1z, m2z, color="purple", linewidth=1.5, alpha=0.7)
    ax3.scatter([m1z[0]], [m2z[0]], color="green", s=50, zorder=10, label="Start")
    ax3.scatter([m1z[-1]], [m2z[-1]], color="red", s=50, zorder=10, label="End")
    # Show antiferromagnetic equilibrium
    ax3.axline((1, -1), slope=-1, color="gray", linestyle="--", alpha=0.5, label="AF line")
    ax3.set_xlabel("Layer 1 $m_z$")
    ax3.set_ylabel("Layer 2 $m_z$")
    ax3.set_title("Phase Space (IEC Coupling)")
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.set_aspect("equal")

    # Time domain plots
    ax4 = fig.add_subplot(2, 3, 4)
    ax4.plot(time, m1z, color="crimson", linewidth=1.5, label="Layer 1", alpha=0.8)
    ax4.plot(time, m2z, color="navy", linewidth=1.5, label="Layer 2", alpha=0.8)
    ax4.axvspan(0, field_duration * 1e9, alpha=0.2, color="yellow", label="Field pulse")
    ax4.set_xlabel("Time (ns)")
    ax4.set_ylabel("$m_z$")
    ax4.set_title("Out-of-plane Magnetization")
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    # In-plane components
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.plot(time, m1x, color="crimson", linewidth=1.5, label="Layer 1 $m_x$", alpha=0.8)
    ax5.plot(time, m2x, color="navy", linewidth=1.5, label="Layer 2 $m_x$", alpha=0.8)
    ax5.set_xlabel("Time (ns)")
    ax5.set_ylabel("$m_x$")
    ax5.set_title("In-plane x-component")
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    # Relative angle between layers
    ax6 = fig.add_subplot(2, 3, 6)
    # Compute dot product m1 · m2
    dot_product = m1x * m2x + m1y * m2y + m1z * m2z
    angle = np.arccos(np.clip(dot_product, -1, 1)) * 180 / np.pi  # Convert to degrees
    ax6.plot(time, angle, color="forestgreen", linewidth=2)
    ax6.axhline(y=180, color="red", linestyle="--", alpha=0.5, label="Antiparallel")
    ax6.axhline(y=0, color="blue", linestyle="--", alpha=0.5, label="Parallel")
    ax6.set_xlabel("Time (ns)")
    ax6.set_ylabel("Angle (degrees)")
    ax6.set_title("Interlayer Angle")
    ax6.legend(fontsize=8)
    ax6.grid(True, alpha=0.3)

    # Add parameter information
    param_text = (
        f"IEC Parameters:\n"
        f"$J_1$ = {J_linear*1e3:.2f} mJ/m²\n"
        f"$J_2$ = {J_quad*1e3:.2f} mJ/m²\n"
        f"$M_{{s,1}}$ = {Ms1:.1f} T\n"
        f"$M_{{s,2}}$ = {Ms2:.1f} T\n"
        f"$K_{{u,1}}$ = {Ku1/1e3:.0f} kJ/m³\n"
        f"$K_{{u,2}}$ = {Ku2/1e3:.0f} kJ/m³"
    )
    fig.text(
        0.98,
        0.02,
        param_text,
        transform=fig.transFigure,
        fontsize=8,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.5),
    )

    fig.suptitle("Interlayer Exchange Coupling Dynamics", fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(
        "./curated-examples/figures/iec-coupling.png",
        dpi=300,
        bbox_inches="tight",
    )

print(f"\nSimulation completed: {sim_time*1e9:.1f} ns")
print(f"Final layer 1 mz: {m1z[-1]:.4f}")
print(f"Final layer 2 mz: {m2z[-1]:.4f}")
print(f"Final interlayer angle: {angle[-1]:.2f} degrees")
print(f"Initial interlayer angle: {angle[0]:.2f} degrees")
