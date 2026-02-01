"""
Magnetization Switching and Hysteresis Loop Simulation

GOAL:
This simulation demonstrates basic magnetization switching dynamics in a single
ferromagnetic layer under an applied magnetic field sweep. It captures the hysteretic
behavior characteristic of magnetic memory elements and illustrates the fundamental
relationship between applied field and magnetization state.

CRITICAL SIMULATION MECHANISMS:
1. Single Layer Dynamics: Models a free magnetic layer with perpendicular magnetic
   anisotropy (PMA) that can switch between up and down states
2. Field Sweep Protocol: Applies ramped external magnetic field to drive magnetization
   reversal through coherent rotation or domain nucleation/propagation
3. Hysteresis Behavior: Captures history-dependent switching with distinct coercive
   fields for up-to-down and down-to-up transitions
4. Dynamic Switching: Uses full LLG equation integration to capture realistic
   switching trajectories including precession and damping effects
5. Anisotropy Energy: Includes uniaxial perpendicular anisotropy that creates
   energy barriers and defines the easy axis for magnetization
6. Demagnetization: Accounts for shape anisotropy that favors in-plane orientation
   for thin films, competing with perpendicular anisotropy

The simulation generates magnetization vs field (M-H) hysteresis loops and
time-domain switching trajectories. This enables characterization of coercive
field, switching field distribution, and energy barrier heights - critical
parameters for magnetic memory (MRAM) and sensor applications.
"""

import contextlib

import matplotlib.pyplot as plt
import numpy as np

from cmtj import AxialDriver, CVector, Junction, Layer, ScalarDriver, constantDriver

with contextlib.suppress(ImportError):
    import scienceplots  # noqa: F401

# Layer parameters from documentation (docs/experimental-methods/examples.ipynb)
# Ms in Tesla for core Layer objects
Ms = 1.03  # Saturation magnetization [T] - from examples.ipynb

# In-plane anisotropy - using LOW value from experimental examples
Ku = 0.8e3  # Anisotropy constant [J/m^3] - from examples.ipynb (0.8 kJ/m³)

# Anisotropy direction - IN-PLANE (along x) for in-plane hysteresis
Kdir = CVector(1, 0, 0)

# Damping from documentation
damping = 0.024  # Gilbert damping - from examples.ipynb

# Standard thin film demagnetization tensor from documentation
demag = [CVector(0, 0, 0), CVector(0, 0, 0), CVector(0, 0, 1.0)]

# Create single magnetic layer with in-plane initial magnetization
# Start slightly off-axis to allow switching
layer = Layer(
    "free",
    mag=CVector(1.0, 0.1, 0.1),  # Start near +x with small perturbation
    anis=Kdir,
    Ms=Ms,
    thickness=2.1e-9,  # From examples.ipynb
    damping=damping,
    demagTensor=demag,
    cellSurface=15e-9 * 15e-9 * np.pi,  # From examples.ipynb
)

# Set constant anisotropy
layer.setAnisotropyDriver(constantDriver(Ku))

# Create junction
junction = Junction([layer])

# Field sweep parameters - IN-PLANE field along easy axis (x) for switching
H_max = 100e3  # Maximum field [A/m] in x direction
H_min = -100e3  # Minimum field [A/m]
H_steps = 60  # Number of field steps

# Create field sweep: down sweep first, then up sweep
field_down = np.linspace(H_max, H_min, H_steps)
field_up = np.linspace(H_min, H_max, H_steps)
field_sweep = np.concatenate([field_down, field_up])

# Time parameters
# Following AGENTS.md: use dt=1e-12 for standard simulations
dt = 1e-12
relax_time = 20e-9  # 20 ns relaxation - from examples.ipynb

# Storage for results - track mx since field is in x direction
mx_values = []
applied_fields = []

print("Running hysteresis loop simulation...")
print(f"Ms = {Ms:.2f} T")
print(f"Ku = {Ku/1e3:.1f} kJ/m³")
print(f"Field range: {H_min/1e3:.0f} to {H_max/1e3:.0f} kA/m")
print(f"Damping α = {damping:.3f}\n")

# Run field sweep
for i, H_field in enumerate(field_sweep):
    # Set external field in X direction (along easy axis)
    junction.setLayerExternalFieldDriver(
        "free",
        AxialDriver(
            ScalarDriver.getConstantDriver(H_field),
            ScalarDriver.getConstantDriver(0),
            ScalarDriver.getConstantDriver(0),
        ),
    )

    # Relax magnetization at this field
    junction.clearLog()
    junction.runSimulation(relax_time, dt, dt)

    # Get final magnetization state
    log = junction.getLog()
    mx_final = log["free_mx"][-1]

    mx_values.append(mx_final)
    applied_fields.append(H_field)

    # Progress indicator
    if (i + 1) % 20 == 0:
        print(f"Completed {i+1}/{len(field_sweep)} field points")

# Convert to arrays
mx_values = np.array(mx_values)
applied_fields = np.array(applied_fields)

# Split into down and up sweeps
split_idx = H_steps
H_down = applied_fields[:split_idx]
H_up = applied_fields[split_idx:]
mx_down = mx_values[:split_idx]
mx_up = mx_values[split_idx:]

# Find coercive fields (where mx crosses zero)
# Down sweep: find where mx goes from positive to negative
try:
    idx_down = np.where(np.diff(np.sign(mx_down)))[0][0]
    Hc_down = H_down[idx_down]
except IndexError:
    Hc_down = np.nan

# Up sweep: find where mx goes from negative to positive
try:
    idx_up = np.where(np.diff(np.sign(mx_up)))[0][0]
    Hc_up = H_up[idx_up]
except IndexError:
    Hc_up = np.nan

print(f"\nCoercive field (down sweep): {Hc_down/1e3:.1f} kA/m")
print(f"Coercive field (up sweep): {Hc_up/1e3:.1f} kA/m")
print(f"Coercivity asymmetry: {abs(Hc_down - Hc_up)/1e3:.1f} kA/m")

# Create simplified M(H) plot
with plt.style.context(["science", "no-latex"]):
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    # Main hysteresis loop (mx vs H)
    ax.plot(H_down / 1e3, mx_down, "-", color="crimson", linewidth=2, label="Down sweep")
    ax.plot(H_up / 1e3, mx_up, "-", color="navy", linewidth=2, label="Up sweep")
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
    if not np.isnan(Hc_down):
        ax.axvline(x=Hc_down / 1e3, color="red", linestyle=":", alpha=0.7, label=f"$H_c$↓={Hc_down/1e3:.1f} kA/m")
    if not np.isnan(Hc_up):
        ax.axvline(x=Hc_up / 1e3, color="blue", linestyle=":", alpha=0.7, label=f"$H_c$↑={Hc_up/1e3:.1f} kA/m")
    ax.set_xlabel(r"$H$ (kA/m)")
    ax.set_ylabel(r"$m_x$")
    ax.set_title(f"$M_s$ = {Ms:.2f} T, $K_u$ = {Ku/1e3:.1f} kJ/m³, $\\alpha$ = {damping:.3f}")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-1.1, 1.1])

    fig.tight_layout()
    fig.savefig(
        "./curated-examples/figures/switching-hysteresis.png",
        dpi=300,
        bbox_inches="tight",
    )

print("\nSimulation complete!")
print(f"Hysteresis loop figure saved to curated-examples/figures/switching-hysteresis.png")
