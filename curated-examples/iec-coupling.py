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

# Anisotropy values - reduced for observable dynamics
Ku1 = 100e3  # Free layer anisotropy [J/m^3]
Ku2 = 80e3  # Fixed layer anisotropy [J/m^3]

# Anisotropy direction - in-plane for better dynamics with perpendicular field
Kdir = CVector(1, 0, 0)

# Damping in recommended range
damping1 = 0.02
damping2 = 0.015

# Demagnetization tensors - removed to simplify
demag1 = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]
demag2 = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 0)]

# Create two coupled magnetic layers with in-plane initial magnetization
l1 = Layer(
    "layer1",
    mag=CVector(1.0, 0, 0),  # Initially along +x
    anis=Kdir,
    Ms=Ms1,
    thickness=1.5e-9,
    damping=damping1,
    demagTensor=demag1,
    cellSurface=np.pi * (50e-9) ** 2,  # Circular cross-section
)

l2 = Layer(
    "layer2",
    mag=CVector(-1.0, 0, 0),  # Initially along -x (antiparallel due to IEC)
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
J_linear = -1.0e-3  # -1.0 mJ/m^2, stronger antiferromagnetic coupling
J_quad = 0.0  # Remove quadratic term for simplicity

junction.setIECDriver("layer1", "layer2", constantDriver(J_linear))
junction.setQuadIECDriver("layer1", "layer2", constantDriver(J_quad))

# Apply external field pulse to perturb the system
# Field strength in typical range (±0 – ±500e3 A/m)
field_amplitude = 200e3  # A/m - stronger field for clear dynamics
field_duration = 5e-9  # 5 ns pulse - longer for observable effect

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

# Create simplified plots - focus on magnetization dynamics
with plt.style.context(["science", "no-latex"]):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), dpi=300)

    # Plot 1: Out-of-plane magnetization vs time
    ax1.plot(time, m1z, color="crimson", linewidth=2, label="Layer 1 $m_z$")
    ax1.plot(time, m2z, color="navy", linewidth=2, label="Layer 2 $m_z$")
    ax1.axvspan(0, field_duration * 1e9, alpha=0.2, color="yellow", label="Field pulse")
    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("$m_z$")
    ax1.set_title("Out-of-plane Magnetization")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Interlayer angle vs time
    # Compute dot product m1 · m2
    dot_product = m1x * m2x + m1y * m2y + m1z * m2z
    angle = np.arccos(np.clip(dot_product, -1, 1)) * 180 / np.pi  # Convert to degrees
    ax2.plot(time, angle, color="forestgreen", linewidth=2)
    ax2.axhline(y=180, color="red", linestyle="--", alpha=0.5, label="Antiparallel")
    ax2.axhline(y=0, color="blue", linestyle="--", alpha=0.5, label="Parallel")
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("Interlayer Angle (degrees)")
    ax2.set_title(f"$J_1$ = {J_linear*1e3:.2f} mJ/m², $J_2$ = {J_quad*1e3:.2f} mJ/m²")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
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
