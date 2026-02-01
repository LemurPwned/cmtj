"""
Interlayer Exchange Coupling (IEC) Coupled Oscillations Simulation

GOAL:
This simulation demonstrates sustained coupled oscillations in two IEC-coupled magnetic
layers, showing acoustic and optical precession modes. The IEC coupling creates collective
dynamics where the layers oscillate together (acoustic mode) or in opposition (optical mode),
which is essential for coupled spin-torque oscillators and synchronized switching in
spintronic devices.

CRITICAL SIMULATION MECHANISMS:
1. IEC Interaction: Models antiferromagnetic coupling (negative J) that creates 
   effective fields between layers proportional to their relative orientations
2. Dual-Layer System: Two ferromagnetic layers with different properties coupled
   through IEC, enabling acoustic and optical oscillation modes
3. Continuous Drive: Applies oscillating external field to sustain coupled precession
   and demonstrate long-term coupled dynamics
4. Coupled Modes: Captures both in-phase (acoustic) and out-of-phase (optical)
   magnetization oscillations arising from symmetric/antisymmetric layer combinations
5. Phase Relationships: Shows how IEC coupling maintains phase coherence between
   layers during sustained oscillations
6. Low Damping: Uses realistic low damping to enable observable coupled dynamics

The simulation generates magnetization time traces showing sustained coupled oscillations.
This demonstrates IEC-mediated collective behavior crucial for coupled STOs,
synchronized switching, and enhanced stability in SAF-based devices.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np

from cmtj import AxialDriver, CVector, Junction, Layer, constantDriver

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Layer parameters  (docs/tutorials)
# Ms in Tesla for core Layer objects
Ms1 = 1.0  # Free layer - typical from trajectory.ipynb
Ms2 = 1.2  # Reference layer - typical from trajectory.ipynb

# Anisotropy values 
Ku1 = 300e3  # Free layer [J/m^3]
Ku2 = 800e3  # Reference layer [J/m^3]

# Anisotropy direction - perpendicular for clear oscillations
Kdir = CVector(0, 0, 1)

# Low damping for sustained oscillationm
damping1 = 0.011  
damping2 = 0.011

# Standard thin film demagnetization tensor 
demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]

# Create two coupled magnetic layers with perpendicular magnetization
l1 = Layer(
    "layer1",
    mag=CVector(0.0, 0, 1.0),  # Initially up
    anis=Kdir,
    Ms=Ms1,
    thickness=1.4e-9,  # Typical from docs
    damping=damping1,
    demagTensor=demag,
    cellSurface=7e-10 * 7e-10,  # From trajectory.ipynb
)

l2 = Layer(
    "layer2",
    mag=CVector(0.0, 0, -1.0),  # Initially down (antiparallel due to IEC)
    anis=Kdir,
    Ms=Ms2,
    thickness=3e-9,  # Thicker reference layer from trajectory.ipynb
    damping=damping2,
    demagTensor=demag,
    cellSurface=7e-10 * 7e-10,
)

# Set anisotropy drivers
l1.setAnisotropyDriver(constantDriver(Ku1))
l2.setAnisotropyDriver(constantDriver(Ku2))

# Create junction with both layers
junction = Junction([l1, l2])

J_linear = -4e-5  # J/m²
junction.setIECDriver("layer1", "layer2", constantDriver(J_linear))

# Apply continuous oscillating field to sustain coupled oscillations
# Use sinusoidal drive similar to VCMA example in trajectory.ipynb
from cmtj import ScalarDriver

# Oscillation parameters
field_amplitude = 50e3  # A/m - moderate amplitude
oscillation_freq = 5e9  # 5 GHz - typical FMR frequency range
bias_field = 100e3  # A/m - bias field in y direction to set operating point

# Set oscillating field on both layers in y direction
junction.setLayerExternalFieldDriver(
    "layer1",
    AxialDriver(
        ScalarDriver.getConstantDriver(0),
        ScalarDriver.getSineDriver(bias_field, field_amplitude, oscillation_freq, 0),
        ScalarDriver.getConstantDriver(0),
    ),
)

junction.setLayerExternalFieldDriver(
    "layer2",
    AxialDriver(
        ScalarDriver.getConstantDriver(0),
        ScalarDriver.getSineDriver(bias_field, field_amplitude, oscillation_freq, 0),
        ScalarDriver.getConstantDriver(0),
    ),
)

dt = 1e-12
sim_time = 30e-9  # 30 ns to show sustained oscillations (within 1-500 ns range)
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

# Create simplified plots - focus on sustained oscillations
with plt.style.context(["science", "no-latex"]):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), dpi=300)

    # Plot 1: In-plane magnetization vs time showing sustained oscillations
    ax1.plot(time, m1y, color="crimson", linewidth=1.5, label="Layer 1 $m_y$")
    ax1.plot(time, m2y, color="navy", linewidth=1.5, label="Layer 2 $m_y$")
    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("$m_y$")
    ax1.set_title("Sustained Coupled Oscillations")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Interlayer angle vs time
    # Compute dot product m1 · m2
    dot_product = m1x * m2x + m1y * m2y + m1z * m2z
    angle = np.arccos(np.clip(dot_product, -1, 1)) * 180 / np.pi  # Convert to degrees
    ax2.plot(time, angle, color="forestgreen", linewidth=1.5)
    ax2.axhline(y=180, color="red", linestyle="--", alpha=0.5, label="Antiparallel")
    ax2.axhline(y=0, color="blue", linestyle="--", alpha=0.5, label="Parallel")
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("Interlayer Angle (degrees)")
    ax2.set_title(f"$J_1$ = {J_linear*1e3:.2f} mJ/m², f = {oscillation_freq/1e9:.1f} GHz")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(
        "./curated-examples/figures/iec-coupling.png",
        dpi=300,
        bbox_inches="tight",
    )

print(f"\nSimulation completed: {sim_time*1e9:.1f} ns")
print(f"Oscillation frequency: {oscillation_freq/1e9:.1f} GHz")
print(f"Final layer 1 my: {m1y[-1]:.4f}")
print(f"Final layer 2 my: {m2y[-1]:.4f}")
print(f"Final interlayer angle: {angle[-1]:.2f} degrees")
print(f"Mean interlayer angle: {np.mean(angle):.2f} degrees")
print(f"Oscillation amplitude (layer 1): {np.std(m1y):.4f}")
print(f"Oscillation amplitude (layer 2): {np.std(m2y):.4f}")
